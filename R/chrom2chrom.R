
library(dplyr)
library(circlize)

load("results/edge_lists.RData")
load("results/noteDF.RData")

df <- tibble(from = left_join(inData, noteDF, by = c("g1" = "id"))$chr, 
             to = left_join(inData, noteDF, by = c("g2" = "id"))$chr)

df <- df %>%
  group_by(from, to) %>%
  summarise(value = n()) %>%
  filter(!is.na(from) & !is.na(to))

chordDiagram(df)

# Apply log10 + 1 to counts
df$value <- log10(df$value) + 1

chordDiagram(df)
