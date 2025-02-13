library(tidyverse)
library(readxl)
cirr_csv <- read_csv("../data/manual/SraRunTable.txt")
cirr_md <- cirr_csv %>%
  mutate(env = ifelse(grepl("^H", `Sample_Name`), "control", "case")) %>%
  separate_wider_delim(`Sample_Name`, delim="_", names=c("subject", "runnum")) %>%
  select(sample=Run, subject, env) %>% mutate(subject=gsub("-","",subject))
write_tsv(cirr_md, "../data/manual/cirr_md.tsv")

