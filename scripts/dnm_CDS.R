# This script creates the bigram network
# for the child-directed speech data.

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

fion_cds <- read_csv("/Users/stefanhartmann/Library/CloudStorage/Dropbox/Input_Project/Data/master/fion_input_with_language_tags.csv.zip")
silvie_cds <- read_csv("/Users/stefanhartmann/Library/CloudStorage/Dropbox/Input_Project/Data/master/silvie_input_with_language_tags.csv.zip")

# get bigrams -------------------------------------------------------------

bigrams_cds_fion <- fion_cds %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)
bigrams_cds_silvie <- silvie_cds %>% unnest_tokens(bigram, Utterance_clean, token = "ngrams", n = 2)


# one column for each word
bigrams_cds_fion <- bigrams_cds_fion %>% separate(bigram, c("word1", "word2"), sep = " ", remove = F)
bigrams_cds_silvie <- bigrams_cds_silvie %>% separate(bigram, c("word1", "word2"), sep = " ", remove = F)

# count:
# first grouped by Utterance number so that
# the bigrams do not cross utterance boundaries,
# then summing up across utterances.
bigrams_fion_cds_count <- bigrams_cds_fion %>% group_by(Utt_no, word1, word2, lang) %>% summarise(
  n = n()
) %>% na.omit() %>% group_by(word1, word2, lang) %>%
  summarise(
    n = sum(n)
  ) %>% arrange(desc(n))

bigrams_silvie_cds_count <- bigrams_cds_silvie %>% group_by(Utt_no, word1, word2, lang) %>% summarise(
  n = n()
) %>% na.omit() %>% group_by(word1, word2, lang) %>%
  summarise(
    n = sum(n)
  ) %>% arrange(desc(n))

# bigrams_fion_cds_count <- bigrams_cds_fion %>% count(word1, word2, sort = T)
# bigrams_silvie_cds_count <- bigrams_cds_silvie %>% count(word1, word2, sort = T)

# get graphs

# add language tag to the words
bigrams_fion_cds_count$word1 <- paste0(bigrams_fion_cds_count$word1, "_", bigrams_fion_cds_count$lang)
bigrams_fion_cds_count$word2 <-paste0(bigrams_fion_cds_count$word2, "_", bigrams_fion_cds_count$lang)

fion_cds_graph <- bigrams_fion_cds_count %>% 
  na.omit %>% filter(n>5) %>% 
  select(word1, word2, n) %>%
  graph_from_data_frame(directed = F)

# color according to language tag
V(fion_cds_graph)$color <- ifelse(grepl("_en", V(fion_cds_graph)$name), "blue", "green")

# remove language tags from the words
V(fion_cds_graph)$name <- gsub("_en$|_de$", "", V(fion_cds_graph)$name)

# Silvie ------

# add language tag to the words
bigrams_silvie_cds_count$word1 <- paste0(bigrams_silvie_cds_count$word1, "_", bigrams_silvie_cds_count$lang)
bigrams_silvie_cds_count$word2 <-paste0(bigrams_silvie_cds_count$word2, "_", bigrams_silvie_cds_count$lang)


silvie_cds_graph <- bigrams_silvie_cds_count %>% 
  na.omit %>% filter(n>5) %>% 
  select(word1, word2, n) %>%
  graph_from_data_frame(directed = F)

# color according to language tag
V(silvie_cds_graph)$color <- ifelse(grepl("_en", V(silvie_cds_graph)$name), "blue", "green")

# remove language tags from the words
V(silvie_cds_graph)$name <- gsub("_en$|_de$", "", V(silvie_cds_graph)$name)

# add weights
fion_cds_graph <- set.edge.attribute(fion_cds_graph, "weight", value = filter(na.omit(bigrams_fion_cds_count), n > 5)$n)
silvie_cds_graph <- set.edge.attribute(silvie_cds_graph, "weight", value = filter(na.omit(bigrams_silvie_cds_count), n > 5)$n)


# node size by freqency
unigrams_fion_cds <- fion_cds %>% unnest_tokens(output = "unigram", input = Utterance_clean, token = "ngrams", n = 1)  %>%
  group_by(unigram) %>%
  summarise(
    n = n()
  )

unigrams_silvie_cds <- silvie_cds %>% unnest_tokens(output = "unigram", input = Utterance_clean, token = "ngrams", n = 1)  %>%
  group_by(unigram) %>%
  summarise(
    n = n()
  )

# add sizes to nodes
V(fion_cds_graph)$size <- log1p(left_join(tibble(unigram = V(fion_cds_graph)$name), unigrams_fion_cds)$n)
V(silvie_cds_graph)$size <- log1p(left_join(tibble(unigram = V(silvie_cds_graph)$name), unigrams_silvie_cds)$n)

plot(fion_cds_graph)




# save as gexf
# saveAsGEXF(fion_cds_graph, "fion_cds_graph.gexf")
# saveAsGEXF(silvie_cds_graph, "silvie_cds_graph.gexf")



# png("communities_silvie_mono_and_mixed.png", width = 12, height = 5, un = "in", res = 300)
par(mfrow = c(1,3))
plot(my_communities[[1]], my_plots[[1]]) + title(levels(factor(d_silvie$Months))[1])
plot(my_communities[[2]], my_plots[[2]]) + title(levels(factor(d_silvie$Months))[2])
plot(my_communities[[3]], my_plots[[3]]) + title(levels(factor(d_silvie$Months))[3])
# dev.off()
par(mfrow = c(1,1))



# # save as gexf
# saveAsGEXF(fion_cds_graph, "fion_cds_graph_color.gexf")
# saveAsGEXF(silvie_cds_graph, "silvie_cds_graph_color.gexf")
# 
# 
