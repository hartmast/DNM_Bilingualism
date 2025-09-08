# This script identifies the recurrent
# trigrams in the largest communities
# identified by the Louvain method.

# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(tidytext)
library(stringi)
library(ggsankey)
library(patchwork)

# read data
d_fion <- read_csv("../../master/fion_CHI.csv")
d_silvie <- read_csv("../../master/silvie_CHI.csv")

# bin months in three periods
d_fion$Months <- character(nrow(d_fion))
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[1:7], "02;03-02;09", d_fion$Months)
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[8:14], "02;10-03;04", d_fion$Months)
d_fion$Months <- ifelse(d_fion$Month %in% levels(factor(d_fion$Month))[15:21], "03;05-03;11", d_fion$Months)



# read data with community information ------------------------------------

# Fion
f01 <- read_csv("../dnm_graphs/Fion/mono & mixed/fion_mono_mixed_01.csv")
f02 <- read_csv("../dnm_graphs/Fion/mono & mixed/fion_mono_mixed_02.csv")
f03 <- read_csv("../dnm_graphs/Fion/mono & mixed/fion_mono_mixed_03.csv")


# get trigrams from d_fion, one table per time slice
for(j in 1:3) {
  
  d_fion_cur <- filter(d_fion, Months == levels(factor(d_fion$Months))[j])
  
  for(i in 1:length(levels(factor(d_fion_cur$Filename)))) {
    cur <- filter(d_fion_cur, Filename == levels(factor(d_fion_cur$Filename))[i])
    cur <- unnest_tokens(cur, "trigram", 
                         Utterance_clean, token = "ngrams", n = 3)
    
    if(i == 1) {
      all_trigrams_fion <- cur$trigram
    } else {
      all_trigrams_fion <- c(all_trigrams_fion, cur$trigram)
    }
    
    print(i)
  }
  
  # get table of trigrams
  all_tbl_fion <- table(all_trigrams_fion) %>% as_tibble(.name_repair = "unique") %>%
    setNames(c("trigram", "n")) %>% arrange(desc(n))
  
  all_tbl_fion <- separate_wider_delim(data = all_tbl_fion,
                                       cols = trigram,
                                       delim = " ",
                                       #too_few = "align_start",
                                       cols_remove = FALSE,
                                       names = c("gram1", "gram2", "gram3"))
  
  assign(paste0("all_tbl_fion_0", j), all_tbl_fion)
  
}


# find the four largest communities
# in each of the time slices

top_clusters_fion01 <- f01$att2 %>% table %>% as_tibble(.names_repair="unique") %>%
  setNames(c("cluster", "n")) %>% arrange(desc(n)) %>% top_n(4) %>%
  select(cluster) %>% unname() %>% unlist()

top_clusters_fion02 <- f02$att2 %>% table %>% as_tibble(.names_repair="unique") %>%
  setNames(c("cluster", "n")) %>% arrange(desc(n)) %>% top_n(4) %>%
  select(cluster) %>% unname() %>% unlist()

top_clusters_fion03 <- f03$att2 %>% table %>% as_tibble(.names_repair="unique") %>%
  setNames(c("cluster", "n")) %>% arrange(desc(n)) %>% top_n(4) %>%
  select(cluster) %>% unname() %>% unlist()

# PERIOD 1
for(i in 1:4) {
  # get words in top clusters
  cur <- filter(f01, att2 %in% top_clusters_fion01[i])
  
  # which trigrams in the respective period consist 
  # exclusively of words that belong to the cluster?
  cur_tbl <- all_tbl_fion_01[which(all_tbl_fion_01$gram1 %in% cur$Label & 
                          all_tbl_fion_01$gram2 %in% cur$Label &
                          all_tbl_fion_01$gram3 %in% cur$Label),] %>% 
    mutate(cluster = i) %>%
    mutate(months = levels(factor(d_fion$Months))[1])
  
  if(i == 1) {
    fion_top_clusters01 <- cur_tbl
  } else {
    fion_top_clusters01 <- rbind(fion_top_clusters01, cur_tbl)
  }
}

# PERIOD 2
for(i in 1:4) {
  # get words in top clusters
  cur <- filter(f02, att2 %in% top_clusters_fion02[i])
  
  # which trigrams in the respective period consist 
  # exclusively of words that belong to the cluster?
  cur_tbl <- all_tbl_fion_02[which(all_tbl_fion_02$gram1 %in% cur$Label & 
                                     all_tbl_fion_02$gram2 %in% cur$Label &
                                     all_tbl_fion_02$gram3 %in% cur$Label),] %>% 
    mutate(cluster = i) %>%
    mutate(months = levels(factor(d_fion$Months))[2])
  
  if(i == 1) {
    fion_top_clusters02 <- cur_tbl
  } else {
    fion_top_clusters02 <- rbind(fion_top_clusters02, cur_tbl)
  }
}


# PERIOD 3
for(i in 1:4) {
  # get words in top clusters
  cur <- filter(f03, att2 %in% top_clusters_fion03[i])
  
  # which trigrams in the respective period consist 
  # exclusively of words that belong to the cluster?
  cur_tbl <- all_tbl_fion_03[which(all_tbl_fion_03$gram1 %in% cur$Label & 
                                     all_tbl_fion_03$gram2 %in% cur$Label &
                                     all_tbl_fion_03$gram3 %in% cur$Label),] %>% 
    mutate(cluster = i) %>%
    mutate(months = levels(factor(d_fion$Months))[3])
  
  if(i == 1) {
    fion_top_clusters03 <- cur_tbl
  } else {
    fion_top_clusters03 <- rbind(fion_top_clusters03, cur_tbl)
  }
}


# bind together
fion_top_clusters <- rbind(fion_top_clusters01, fion_top_clusters02, fion_top_clusters03)

# export
# writexl::write_xlsx(fion_top_clusters, "fion_top_clusters.xlsx")



# re-import annotated data ------------------------------------------------

fion <- readxl::read_xlsx("fion_top_clusters.xlsx", sheet = 1)
pos_collapse <- readxl::read_xlsx("fion_top_clusters.xlsx", sheet = 2)

# add columns with minimal annotation
fion <- left_join(fion, select(pos_collapse, "POS", "minimal"), by = c("pos1" = "POS"))
colnames(fion)[which(colnames(fion)=="minimal")] <- "POS1"

fion <- left_join(fion, select(pos_collapse, "POS", "minimal"), by = c("pos2" = "POS"))
colnames(fion)[which(colnames(fion)=="minimal")] <- "POS2"

fion <- left_join(fion, select(pos_collapse, "POS", "minimal"), by = c("pos3" = "POS"))
colnames(fion)[which(colnames(fion)=="minimal")] <- "POS3"

# add column with POS trigrams
fion <- mutate(fion, POS_trigrams = paste0(fion$POS1, " ", fion$POS2, " ", fion$POS3))


# tabulate POS trigrams
fion %>% 
  group_by(months, cluster, POS_trigrams) %>% 
  summarise(
    n = n()
  ) %>% filter(n > 10) %>%
  arrange(desc(n)) %>%
  arrange(cluster) %>%
  arrange(months) #%>% 
   #writexl::write_xlsx("pos_trigrams.xlsx")

library(reactable)
library(reactablefmtr)

# also group by language
data <- fion %>% 
  group_by(months, cluster, POS_trigrams) %>% 
  summarise(
    n = n(),
    de = length(which(language == "de")),
    en = length(which(language == "en")),
    mixed = length(which(language == "mixed")),
    de_rel = de/n,
    en_rel = en/n,
    mixed_rel = mixed/n
  ) %>% 
  # group n<10 to "other"
  mutate(POS_trigrams = ifelse(n<10, "other", POS_trigrams)) %>%
  group_by(months, cluster, POS_trigrams) %>% summarise(
    n = sum(n),
    de = sum(de),
    en = sum(en),
    mixed = sum(mixed),
    de_rel = de/n,
    en_rel = en/n,
    mixed_rel = mixed/n
  )

data


# read table with examples
withex <- readxl::read_xlsx("/Users/stefanhartmann/Library/CloudStorage/Dropbox/Input_Project/Data/scripts/dnm_graphs/pos_trigrams_with_examples.xlsx")

# omit whitespace in colnames
colnames(withex) <- gsub(" ", "_", colnames(withex))
withex <- rename(withex, POS_trigrams = POS_Trigram)
withex <- rename(withex, cluster = Cluster)
withex <- rename(withex, months = Age)

d_withex <- left_join(data, withex)


# add percentage
d_withex$d_percent <- paste0(round(d_withex$de_rel*100, digits = 2), "%")
d_withex$e_percent <- paste0(round(d_withex$en_rel*100, digits = 2), "%")
d_withex$m_percent <- paste0(round(d_withex$mixed_rel*100, digits = 2), "%")

d_withex #%>% writexl::write_xlsx("pos_trigrams_with_examples_and_percentages.xlsx")




# adverbs in Fion's utterances --------------------------------------------

d_fion_adv <- d_fion %>% unnest_tokens(output = "word", input = "Utterance_clean") %>% group_by(Month, word) %>% summarise(
  Freq = n()
) %>% ungroup %>% add_count(Month) %>% 
  mutate(rel = Freq/n) %>% 
  filter(word %in% c("jetzt", "now", "da", "there", "hier", "here")) %>%
  mutate(lang = ifelse(word %in% c("jetzt", "da", "hier"), "de", "en")) 

d_fion_adv_col <- tibble(words = c("jetzt", "now", "da", "there", "hier", "here"),
                         cols = c("darkblue","lightblue", "darkred", "orangered", "indianred",  "steelblue3"))

d_fion_adv <- left_join(d_fion_adv, d_fion_adv_col, by = c("word" = "words"))

d_fion_adv %>%
  ggplot(aes(x = Month, y = rel, group = word,  col = cols, label = word)) + 
  geom_line() + 
  geom_label(data = filter(d_fion_adv, Month == "02;07" & word %in% c("da", "jetzt", "there"))) +
  geom_label(data = filter(d_fion_adv, Month == "03;04" & word %in% c("here", "hier", "now"))) +
  geom_label(data = filter(d_fion_adv, Month == "03;11" & word %in% c("here"))) +
  guides(col = "none", lty = "none") +
  scale_color_identity()
  
# Sankey diagram of POS tags
unique(c(fion$POS1, fion$POS2, fion$POS3))


# POS factor levels for right display order
pos_factor_levels <- c("other", "Art", "Interj", "Part", "Prep", 
  "Conj", "Card",
   "Adj", "Adv", "Verb", "Pron", "Noun"
  )

# Sankey diagram

for(i in 1:3) {
  assign(paste0("POS_by_month", i), fion %>% 
           filter(months == levels(factor(fion$months))[i]) %>%
           select("POS1", "POS2", "POS3") %>%
           make_long(POS1, POS2, POS3) %>% 
           mutate(node = factor(node, levels = pos_factor_levels)) %>%
           mutate(next_node = factor(next_node, levels = pos_factor_levels)) %>%
           ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
           geom_sankey() +
           scale_fill_discrete(drop=FALSE) +
           geom_sankey_label(aes(label = node)) +
           labs(x = NULL) +
           theme_bw() +
           theme(legend.position = "none",
                 plot.title = element_text(hjust = .5)) +
           scale_fill_viridis_d(begin = .1, end = .9) +
           ggtitle(levels(factor(fion$months))[i]) +
           theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
           theme(axis.text = element_blank()) +
           theme(axis.ticks = element_blank()) +
           theme(axis.title = element_blank()) +
           theme(strip.text = element_blank()) +
           theme(legend.text = element_text(size = 18)) +
           theme(legend.title = element_text(size = 18, face = "bold")) +
           theme(text = element_text(size = 18))
  )
}


layout <- "
AABB
AABB
CC##
CC##
"

POS_by_month1 + POS_by_month2 + POS_by_month3 + 
  plot_layout(design = layout)
# ggsave("sankey_diagrams.png", width = 14, height = 14)


# most frequent trigrams across periods
d_fion %>% group_by(Months, gram3) %>% summarise(
  n = n()
) %>% arrange(desc(n), Months) %>% na.omit %>%
  top_n(15) %>% arrange(Months) %>% print(n = 46) 


