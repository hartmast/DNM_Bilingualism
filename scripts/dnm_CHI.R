# TODO:
# plots for
# - only German
# - only English
# - in comparison with code-mixed data!

# -  more fine-grained periodization options


# This script creates the bigram networks
# for Fion's monolingual and code-mixed data.


# install CRAN packages (if not yet installed)
#sapply(c("tidyverse", "tidytext", "ngram", "igraph", "ggraph", "grid"), function(x) 
#  if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))

# remotes::install_github("david-barnett/microViz")

# load packages
library(tidyverse)
library(tidytext)
library(ngram)
library(igraph)
library(ggraph)
library(patchwork)
library(svglite)
library(geomnet)
library(rgexf)
library(qgraph)
library(genBaRcode)
library(ggraph)
library(ggiraph)
#install.packages(
 # "microViz",
  #repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("phyloseq")
library(microViz)
library(tidygraph)
library(plotly)

# periods for Fion and Silvie
fion_periods <- c("02;03-02;09", "02;10-03;04", "03;05-03;11")
silvie_periods <- c("02;03-02;09", "02;10-03;03", "03;04-03;19")

# helper functions ---------------------------------------------------------

# adapted from https://gk.palem.in/iGraphExport.html
# Converts the given igraph object to GEXF format and saves it at the given filepath location
#     g: input igraph object to be converted to gexf format
#     filepath: file location where the output gexf file should be saved
#
saveAsGEXF = function(g, filepath="converted_graph.gexf") {
  #require(igraph)
  #require(rgexf)
  
  # gexf nodes require two column data frame (id, label)
  # check if the input vertices has label already present
  # if not, just have the ids themselves as the label
  if(is.null(V(g)$label))
    V(g)$label <- as.character(V(g))
  
  # similarily if edges does not have weight, add default 1 weight
  if(is.null(E(g)$weight))
    E(g)$weight <- rep.int(1, ecount(g))
  
  nodes <- data.frame(cbind(V(g), V(g)$label))
  # edges <- ends(g, E(g))
  # get.edge(g, id = 1)
  # ends(g, es = 1, names = FALSE)
  edges <- t(Vectorize(ends, vectorize.args='es')(g, 1:ecount(g), names = FALSE))
  
  # combine all node attributes into a matrix (and take care of & for xml)
  vAttrNames <- setdiff(list.vertex.attributes(g), "label") 
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&", vertex_attr(g, attr))))
  
  # combine all edge attributes into a matrix (and take care of & for xml)
  eAttrNames <- setdiff(list.edge.attributes(g), "weight") 
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&", edge_attr(g, attr))))
  
  # combine all graph attributes into a meta-data
  graphAtt <- sapply(list.graph.attributes(g), function(attr) sub("&", "&", graph_attr(g, attr)))
  
  # generate the gexf object
  output <- write.gexf(nodes, edges, 
                       edgesWeight=E(g)$weight,
                       edgesAtt = edgesAtt,
                       nodesAtt = nodesAtt,
                       meta=c(list(creator="Authors - using a script by Gopalakrishna Palem", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
  
  sink(file = filepath)
  print(output, filepath, replace=T)
  sink()
}


# function for correcting colors in the GEXF file by inserting RGB values
# (calculated in the script below with the help of genBaRcode's .hex2rgb function)
correct_colors <- function(myfile, mynetwork) {
  # read in the file again and correct the colors
  cur_gexf <- readLines(myfile)
  find_cols <- grep("<viz:color", cur_gexf)
  
  for(i in 1:length(find_cols)) {
    # split up current rgb line in gexf file
    cur_splitted_rgb <- unlist(strsplit(cur_gexf[find_cols[i]], "(?<=[rgb]=\").*?\"", perl = T))
    
    # split up rgb value to be inserted
    cur_rgb_values <- gsub("'", "", unlist(strsplit(V(mynetwork)$color[i], ",")))
    
    # paste together
    cur_pasted <- paste0(cur_splitted_rgb[1], cur_rgb_values[1], "\"",
                         cur_splitted_rgb[2], cur_rgb_values[2], "\"",
                         cur_splitted_rgb[3], cur_rgb_values[3], "\"", cur_splitted_rgb[4])
    
    # replace in gexf
    cur_gexf[find_cols[i]] <- cur_pasted
    
    # counter
    #print(i)
    
    
  }
  # output
  return(cur_gexf)
}



# monolingual + mixed data --------------------------------------------------


# read data ---------------------------------------------------------------

d_fion <- read_csv("../../master/fion_CHI.csv")
d_silvie <- read_csv("../../master/silvie_CHI.csv")

# data wrangling ----------------------------------------------------------


# add wordcount
d_fion$wordcount <- sapply(1:nrow(d_fion), 
                           function(i) wordcount(d_fion$Utterance_clean[i]))


d_silvie$wordcount <- sapply(1:nrow(d_silvie), 
                             function(i) wordcount(d_silvie$Utterance_clean[i]))


# bin months in three periods
d_fion$Months <- character(nrow(d_fion))
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[1:7], "02;03-02;09", d_fion$Months)
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[8:14], "02;10-03;04", d_fion$Months)
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[15:21], "03;05-03;11", d_fion$Months)

d_silvie$Months <- character(nrow(d_silvie))
d_silvie$Months <- ifelse(d_silvie$Month %in% levels(factor(d_silvie$Month))[1:6], "02;03-02;09", d_silvie$Months)
d_silvie$Months <- ifelse(d_silvie$Month %in% levels(factor(d_silvie$Month))[7:12], "02;10-03;03", d_silvie$Months)
d_silvie$Months <- ifelse(d_silvie$Month %in% levels(factor(d_silvie$Month))[13:28], "03;04-03;19", d_silvie$Months)


# only multi-word units
mwu_fion <- filter(d_fion, wordcount > 1)
mwu_silvie <- filter(d_silvie, wordcount > 1)



# get bigrams -------------------------------------------------------------

bigrams_fion <- mwu_fion %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)
bigrams_silvie <- mwu_silvie %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)


# add language tags
mwu_fion$Lang_Tags <- gsub("[[:punct:]]", "", mwu_fion$Lang_Tags)
mwu_silvie$Lang_Tags <- gsub("[[:punct:]]", "", mwu_silvie$Lang_Tags)

mwu_fion$Lang_Tags <- sapply(1:nrow(mwu_fion), function(i) ifelse(is.na(mwu_fion$Lang_Tags[i]), ifelse(mwu_fion[i,]$type=="german", paste0(rep("g", mwu_fion[i,]$wordcount), collapse = " "), paste0(rep("e", mwu_fion[i,]$wordcount), collapse = " ")), mwu_fion$Lang_Tags[i]))
mwu_silvie$Lang_Tags <- sapply(1:nrow(mwu_silvie), function(i) ifelse(is.na(mwu_silvie$Lang_Tags[i]), ifelse(mwu_silvie[i,]$type=="german", paste0(rep("g", mwu_silvie[i,]$wordcount), collapse = " "), paste0(rep("e", mwu_silvie[i,]$wordcount), collapse = " ")), mwu_silvie$Lang_Tags[i]))



bigrams_fion <- bind_cols(bigrams_fion,
                          mwu_fion %>% unnest_tokens(bigram_LangTag, Lang_Tags, token = "ngrams", n = 2, drop = FALSE) %>% select(bigram_LangTag))


bigrams_silvie <- bind_cols(bigrams_silvie,
                            mwu_silvie %>% unnest_tokens(bigram_LangTag, Lang_Tags, token = "ngrams", n = 2, drop = FALSE) %>% select(bigram_LangTag))



# one column for each word
bigrams_fion <- bigrams_fion %>% separate(bigram, c("word1", "word2"), sep = " ", remove = F)
bigrams_silvie <- bigrams_silvie %>% separate(bigram, c("word1", "word2"), sep = " ", remove = F)

bigrams_fion <- bigrams_fion %>% separate(bigram_LangTag, c("LangTag1", "LangTag2"), sep = " ", remove = F)
bigrams_silvie <- bigrams_silvie %>% separate(bigram_LangTag, c("LangTag1", "LangTag2"), sep = " ", remove = F)


# add child column
bigrams_fion <- mutate(bigrams_fion, Child = "Fion")
bigrams_silvie <- mutate(bigrams_silvie, Child = "Silvie")

# add together
bigrams <- rbind(bigrams_fion, bigrams_silvie)

# count
bigrams_count <- bigrams %>% group_by(Child, Months, LangTag1, LangTag2) %>% count(word1, word2, sort = T)
#bigrams_count <- bigrams %>% group_by(Child, Months) %>% count(word1, word2, sort = T)




# a) for Fion - monolingual + code-mixed

my_plots <- list()
my_communities <- list()


for(i in 1:length(levels(factor(d_fion$Months)))) {
  
  # code from Silge & Robinson, Tidy Text Mining with R, CC-BY-NC-3.0,
  # https://www.tidytextmining.com/ngrams.html
  
  # check if there are data
  l <- bigrams_count %>%
    filter(n > 5, Child == "Fion",
           Months == levels(factor(d_fion$Months))[i])
  
  if(nrow(l) > 0) {
    # get bigram graph edges and vertices
    bigram_graph <- bigrams_count %>%
      filter(n > 5, Child == "Fion",
             Months == levels(factor(d_fion$Months))[i]) %>%
      ungroup %>% select(word1, word2, n) %>% graph_from_data_frame(directed = FALSE)
    
    # set weight attributes
    bigram_graph <- set_edge_attr(bigram_graph, "weight", value = l$n)
    
    # set labels
    V(bigram_graph)$label <- V(bigram_graph)$name
    
    # Louvain clustering
    lv <- cluster_louvain(bigram_graph)
    
    # save
    my_plots[[i]] <- bigram_graph
    my_communities[[i]] <- lv
    
    
  }
  
}


# node size by freqency (incl. langTag)
unigrams_fion <- bind_cols(mwu_fion %>% unnest_tokens(output = "unigram", input = Utterance_clean, token = "ngrams", n = 1),
                           mwu_fion %>% unnest_tokens(output = "unigram_LangTag", input = Lang_Tags, token = "ngrams", n = 1) %>% select(unigram_LangTag)) %>%
  group_by(Months, unigram, unigram_LangTag) %>%
  summarise(
    n = n()
  )

# add language color column
unigrams_fion$lang_color <- case_when(unigrams_fion$unigram_LangTag == "e" ~ "blue",
                                 unigrams_fion$unigram_LangTag == "g" ~ "orange",
                                 .default = "grey")




# avoid duplicates
unigrams_fion <- unigrams_fion[which(!duplicated(select(unigrams_fion, Months, unigram))),]

# add label color column based on community membership

# color palettes
mycols <- microViz::distinct_palette(n = 41)
mycols <- mycols[which(mycols != "lightgrey")]
mycols <- c(mycols, "#D3D3D3", "#A9A9A9", "#696969")


V(my_plots[[1]])$color <- mycols[1:max(na.omit(membership(my_communities[[1]])))][membership(my_communities[[1]])]
V(my_plots[[2]])$color <- mycols[1:max(na.omit(membership(my_communities[[2]])))][membership(my_communities[[2]])]
V(my_plots[[3]])$color <- mycols[1:max(na.omit(membership(my_communities[[3]])))][membership(my_communities[[3]])]


# convert colors to gephi-usable rgb codes
# V(my_plots[[1]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[1]])$color)
# V(my_plots[[2]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[2]])$color)
# V(my_plots[[3]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[3]])$color)


# add log frequency
unigrams_fion$LogFreq <- log1p(unigrams_fion$n)

# get size
mySizes1 <- left_join(tibble(unigram = V(my_plots[[1]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[1]))
mySizes2 <- left_join(tibble(unigram = V(my_plots[[2]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[2]))
mySizes3 <- left_join(tibble(unigram = V(my_plots[[3]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[3]))

# node size by frequency - squared so that differences become
# actually visible in the plot (but are not as extreme as they
# would be when using the raw frequency values)
V(my_plots[[1]])$size <- mySizes1$LogFreq^2
V(my_plots[[2]])$size <- mySizes2$LogFreq^2
V(my_plots[[3]])$size <- mySizes3$LogFreq^2

# label color according to language
# V(my_plots[[1]])$label.color <- mySizes1$lang_color
# V(my_plots[[2]])$label.color <- mySizes2$lang_color
# V(my_plots[[3]])$label.color <- mySizes3$lang_color

V(my_plots[[1]])$language <- mySizes1$unigram_LangTag
V(my_plots[[2]])$language <- mySizes2$unigram_LangTag
V(my_plots[[3]])$language <- mySizes3$unigram_LangTag


# save as Gephi files
# saveAsGEXF(my_plots[[1]], "fion01_monoandmixed.gexf")
# saveAsGEXF(my_plots[[2]], "fion02_monoandmixed.gexf")
# saveAsGEXF(my_plots[[3]], "fion03_monoandmixed.gexf")

g <- my_plots[[1]]
lc <- cluster_louvain(g)
imc <- cluster_infomap(g)
plot(lc, g)
modularity(lc)

V(g)$community <- membership(lc)
tg <- as_tbl_graph(g)
ggraph(tg, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(fill = as.factor(language), size = log(size)), shape = 21,  color = "black") +
  geom_node_text(aes(label = ifelse(size > 10, name, ""), size = size), repel = TRUE) +
  theme_graph() +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  guides(size = "none") +
  guides(fill = guide_legend(title = "Language")) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 18, face = "bold")) +
  theme(text = element_text(size = 18)) +
  ggtitle("")


get_plot <- function(myplot, myseed = 1999) {
  g <- myplot
  g <- as.undirected(g, mode = "collapse")
  louvain_clusters <- cluster_louvain(g)
  V(g)$community <- membership(louvain_clusters)
  
  # Convert to tidygraph for ggraph plotting
  tg <- as_tbl_graph(g)
  
  set.seed(myseed)
  mygraph <- ggraph(tg, layout = "fr") +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(fill = as.factor(language), size = log(size)), shape = 21,  color = "black") +
    geom_node_text(aes(label = ifelse(size > 10, name, ""), size = size), repel = TRUE) +
    theme_graph() +
    scale_fill_manual(values = c("blue", "orange", "red")) +
    guides(size = "none") +
    guides(fill = guide_legend(title = "Language")) +
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 18, face = "bold")) +
    theme(text = element_text(size = 18)) +
    ggtitle("")
  
  return(mygraph)
  
}

current_plot <- get_plot(my_plots[[1]])
current_plot %>% ggplotly()
girafe(ggobj = current_plot,
       options = list(opts_hover(reactive = TRUE),
                      opts_zoom(min = 1, max = 5)))

g <- my_plots[[1]]
g <- as.undirected(g, mode = "collapse")
louvain_clusters <- cluster_louvain(g)
V(g)$community <- membership(louvain_clusters)

# Convert to tidygraph for ggraph plotting
tg <- as_tbl_graph(g)

set.seed(1999)
ggraph(tg, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(fill = as.factor(language), size = log(size)), shape = 21,  color = "black") +
  geom_node_text(aes(label = ifelse(size > 10, name, ""), size = size), repel = TRUE) +
  theme_graph() +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  guides(size = "none") +
  guides(fill = guide_legend(title = "Language")) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 18, face = "bold")) +
  theme(text = element_text(size = 18)) +
  ggtitle("")
  

cluster_louvain(my_plots[[1]])

g <- as_undirected(V(my_plots[[1]]), mode = "collapse")
louvain_clusters <- cluster_louvain(g)
plot(my_plots[[1]])


# correct_colors("fion01_monoandmixed.gexf", my_plots[[1]]) %>% writeLines("fion01_monoandmixed.gexf")
# correct_colors("fion02_monoandmixed.gexf", my_plots[[2]]) %>% writeLines("fion02_monoandmixed.gexf")
# correct_colors("fion03_monoandmixed.gexf", my_plots[[3]]) %>% writeLines("fion03_monoandmixed.gexf")


# The rest happens in Gephi, where the following steps have to be taken:
# - choose ForceAtlas2 as visualization option
# - adjust node size (for this purpose, duplicate size column specifying type as "float")
# (outdated: -  show labels and set labels to correspond with node size)
# - adjust label color
# - export different SVG files for English and German (with different outline colors)
# plus one with both languages (without node borders) and overlay them in Inkscape.



# b) for Silvie - monolingual + code-mixed


my_plots <- list()
my_communities <- list()


for(i in 1:length(levels(factor(d_silvie$Months)))) {
  
  # code from Silge & Robinson, Tidy Text Mining with R, CC-BY-NC-3.0,
  # https://www.tidytextmining.com/ngrams.html
  
  # check if there are data
  l <- bigrams_count %>%
    filter(n > 5, Child == "Silvie",
           Months == levels(factor(d_silvie$Months))[i])
  
  if(nrow(l) > 0) {
    # get bigram graph edges and vertices
    bigram_graph <- bigrams_count %>%
      filter(n > 5, Child == "Silvie",
             Months == levels(factor(d_silvie$Months))[i]) %>%
      ungroup %>% select(word1, word2, n) %>% graph_from_data_frame(directed = FALSE)
    
    # set weight attributes
    bigram_graph <- set_edge_attr(bigram_graph, "weight", value = l$n)
    
    # set labels
    V(bigram_graph)$label <- V(bigram_graph)$name
    
    # Louvain clustering
    lv <- cluster_louvain(bigram_graph)
    
    # save
    my_plots[[i]] <- bigram_graph
    my_communities[[i]] <- lv
    
    
  }
  
}


# node size by freqency (incl. langTag)
unigrams_silvie <- bind_cols(mwu_silvie %>% unnest_tokens(output = "unigram", input = Utterance_clean, token = "ngrams", n = 1),
                           mwu_silvie %>% unnest_tokens(output = "unigram_LangTag", input = Lang_Tags, token = "ngrams", n = 1) %>% select(unigram_LangTag)) %>%
  group_by(Months, unigram, unigram_LangTag) %>%
  summarise(
    n = n()
  )

# add language color column
unigrams_silvie$lang_color <- case_when(unigrams_silvie$unigram_LangTag == "e" ~ "blue",
                                      unigrams_silvie$unigram_LangTag == "g" ~ "orange",
                                      .default = "grey")




# avoid duplicates
unigrams_silvie <- unigrams_silvie[which(!duplicated(select(unigrams_silvie, Months, unigram))),]

# add label color column based on community membership

# color palettes
mycols <- microViz::distinct_palette(n = 41)
mycols <- mycols[which(mycols != "lightgrey")]
mycols <- c(mycols, "#D3D3D3", "#A9A9A9", "#696969")


V(my_plots[[1]])$color <- mycols[1:max(na.omit(membership(my_communities[[1]])))][membership(my_communities[[1]])]
V(my_plots[[2]])$color <- mycols[1:max(na.omit(membership(my_communities[[2]])))][membership(my_communities[[2]])]
V(my_plots[[3]])$color <- mycols[1:max(na.omit(membership(my_communities[[3]])))][membership(my_communities[[3]])]

# replace NAs
V(my_plots[[1]])$color[which(is.na(V(my_plots[[1]])$color))] <- "#D3D3D3"
V(my_plots[[2]])$color[which(is.na(V(my_plots[[2]])$color))] <- "#D3D3D3"

# convert colors to gephi-usable rgb codes
V(my_plots[[1]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[1]])$color)
V(my_plots[[2]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[2]])$color)
V(my_plots[[3]])$color <- genBaRcode:::.hex2rgbColor(V(my_plots[[3]])$color)


# add log frequency
unigrams_silvie$LogFreq <- log1p(unigrams_silvie$n)

# get size
mySizes1 <- left_join(tibble(unigram = V(my_plots[[1]])$name), filter(unigrams_silvie, Months == levels(factor(mwu_silvie$Months))[1]))
mySizes2 <- left_join(tibble(unigram = V(my_plots[[2]])$name), filter(unigrams_silvie, Months == levels(factor(mwu_silvie$Months))[2]))
mySizes3 <- left_join(tibble(unigram = V(my_plots[[3]])$name), filter(unigrams_silvie, Months == levels(factor(mwu_silvie$Months))[3]))

# node size by frequency - squared so that differences become
# actually visible in the plot (but are not as extreme as they
# would be when using the raw frequency values)
V(my_plots[[1]])$size <- mySizes1$LogFreq^2
V(my_plots[[2]])$size <- mySizes2$LogFreq^2
V(my_plots[[3]])$size <- mySizes3$LogFreq^2

# label color according to language
# V(my_plots[[1]])$label.color <- mySizes1$lang_color
# V(my_plots[[2]])$label.color <- mySizes2$lang_color
# V(my_plots[[3]])$label.color <- mySizes3$lang_color

V(my_plots[[1]])$language <- mySizes1$unigram_LangTag
V(my_plots[[2]])$language <- mySizes2$unigram_LangTag
V(my_plots[[3]])$language <- mySizes3$unigram_LangTag


# save as Gephi files
# saveAsGEXF(my_plots[[1]], "silvie01_monoandmixed.gexf")
# saveAsGEXF(my_plots[[2]], "silvie02_monoandmixed.gexf")
# saveAsGEXF(my_plots[[3]], "silvie03_monoandmixed.gexf")


# use the correct_colors function
correct_colors("silvie01_monoandmixed.gexf", my_plots[[1]]) %>% writeLines("silvie01_monoandmixed.gexf")
correct_colors("silvie02_monoandmixed.gexf", my_plots[[2]]) %>% writeLines("silvie02_monoandmixed.gexf")
correct_colors("silvie03_monoandmixed.gexf", my_plots[[3]]) %>% writeLines("silvie03_monoandmixed.gexf")


# The rest happens in Gephi, where the following steps have to be taken:
# - choose ForceAtlas2 as visualization option
# - adjust node size (for this purpose, duplicate size column specifying type as "float")
# (outdated: -  show labels and set labels to correspond with node size)
# - adjust label color
# - export different SVG files for English and German (with different outline colors)
# plus one with both languages (without node borders) and overlay them in Inkscape.



