# This script creates the bigram networks
# for Fion's code-mixed data.


# install CRAN packages (if not yet installed)
sapply(c("tidyverse", "tidytext", "ngram", "igraph", "ggraph", "grid"), function(x) 
  if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))

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


# helper function ---------------------------------------------------------

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
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))))
  
  # combine all edge attributes into a matrix (and take care of & for xml)
  eAttrNames <- setdiff(list.edge.attributes(g), "weight") 
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))))
  
  # combine all graph attributes into a meta-data
  graphAtt <- sapply(list.graph.attributes(g), function(attr) sub("&", "&",get.graph.attribute(g, attr)))
  
  # generate the gexf object
  output <- write.gexf(nodes, edges, 
                       edgesWeight=E(g)$weight,
                       edgesAtt = edgesAtt,
                       nodesAtt = nodesAtt,
                       meta=c(list(creator="Gopalakrishna Palem", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
  
  sink(file = filepath)
  print(output, filepath, replace=T)
  sink()
}

# read data ---------------------------------------------------------------

d_fion <- read_csv("../../master/fion_mixed.csv")
d_silvie <- read_csv("../../master/silvie_mixed.csv")


# data wrangling ----------------------------------------------------------


# add wordcount
d_fion$wordcount <- sapply(1:nrow(d_fion), 
                           function(i) wordcount(d_fion$Utterance_clean[i]))


d_silvie$wordcount <- sapply(1:nrow(d_silvie), 
                             function(i) wordcount(d_silvie$Utterance_clean[i]))


# add wordcount of LangTags
d_fion$wordcount_LangTags <- sapply(1:nrow(d_fion), 
                                    function(i) wordcount(d_fion$Lang_Tags[i]))


d_silvie$wordcount_LangTags <- sapply(1:nrow(d_silvie), 
                                      function(i) wordcount(d_silvie$Lang_Tags[i]))



# bin months
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

# remove rows without language tags
mwu_fion <- mwu_fion[!is.na(mwu_fion$Lang_Tags),]
mwu_silvie <- mwu_silvie[!is.na(mwu_silvie$Lang_Tags),]

# get bigrams -------------------------------------------------------------

bigrams_fion <- mwu_fion %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)
bigrams_silvie <- mwu_silvie %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)


# add language tags
mwu_fion$Lang_Tags <- gsub("[[:punct:]]", "", mwu_fion$Lang_Tags)
mwu_silvie$Lang_Tags <- gsub("[[:punct:]]", "", mwu_silvie$Lang_Tags)

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


# bind together
bigrams <- rbind(bigrams_fion, bigrams_silvie)

# count
bigrams_count <- bigrams %>% group_by(Child, Months) %>% count(word1, word2, sort = T)





# get bigram plots:


# a) for Fion - code-mixed

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

# add color column for language
unigrams_fion$color <- case_when(unigrams_fion$unigram_LangTag == "e" ~ "blue",
                                 unigrams_fion$unigram_LangTag == "g" ~ "orange",
                                 .default = "grey")


# avoid duplicates
unigrams_fion <- unigrams_fion[which(!duplicated(select(unigrams_fion, Months, unigram))),]


# unigrams_fion1 <- mwu_fion %>% unnest_tokens(output = "unigram", input = Utterance_clean, token = "ngrams", n = 1)  %>%
#   group_by(Months, unigram) %>%
#   summarise(
#     n = n()
#   )

unigrams_fion$LogFreq <- log1p(unigrams_fion$n)


mySizes1 <- left_join(tibble(unigram = V(my_plots[[1]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[1]))
mySizes2 <- left_join(tibble(unigram = V(my_plots[[2]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[2]))
mySizes3 <- left_join(tibble(unigram = V(my_plots[[3]])$name), filter(unigrams_fion, Months == levels(factor(mwu_fion$Months))[3]))


V(my_plots[[1]])$size <- mySizes1$LogFreq
V(my_plots[[2]])$size <- mySizes2$LogFreq
V(my_plots[[3]])$size <- mySizes3$LogFreq


V(my_plots[[1]])$color <- mySizes1$color
V(my_plots[[2]])$color <- mySizes2$color
V(my_plots[[3]])$color <- mySizes3$color


# save as Gephi files
# saveAsGEXF(my_plots[[1]], "fion_mixed01.gexf")
# saveAsGEXF(my_plots[[2]], "fion_mixed02.gexf")
# saveAsGEXF(my_plots[[3]], "fion_mixed03.gexf")
