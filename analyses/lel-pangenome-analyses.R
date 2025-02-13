###
### Analyze EggNOG annotations across entire _L. eligens_ pangenome
###

library(tidyverse)
library(polars)
library(tidypolars)

gpal_100060 <- tidypolars::scan_csv_polars(path.expand("../data/raw/100060_genes_presence-absence_locus.csv"))
gmd <- read_tsv("../data/raw/genomes-all_metadata.tsv")
eggnog_annot <- tidypolars::scan_parquet_polars(path.expand("../data/raw/Parquet/UHGP_eggNOG/uhgp_90_eggNOG/"))
uhgp_90_annot <- tidypolars::scan_parquet_polars(path.expand("../data/raw/Parquet/UHGP_Cluster/uhgp_90_cluster/"))

cast_to_int64 <- function(data, col) {
  data <- data$with_columns(
    pl$col(col)$cast(pl$Int64)
  )
}

message("Making cluster IDs")
iso_only <- gmd %>% filter(Genome_type=="Isolate") %>% select(Genome) %>% deframe()

genes_to_cluster_id <- gpal_100060 %>%
  pivot_longer(-c(1:14)) %>%
  select(Gene, name, value) %>%
  mutate(value = str_split(pattern = "\t", value)) %>%
  rename(query_name=value) %>%
  (\(.) .$explode("query_name")) %>%
  mutate(gene_id = as.numeric(str_replace(str_replace(
    query_name, "GUT_GENOME",""), "_", ""))) %>%
  cast_to_int64("gene_id") %>%
  inner_join(., uhgp_90_annot, by="gene_id") %>%
  select(Gene, cluster_id) %>%
  distinct()

genes_to_cluster_id_tbl <- collect(genes_to_cluster_id, streaming=TRUE)
write_csv(genes_to_cluster_id_tbl, file = "genes_to_cluster_id_100060.csv")

message("Mean found (regular)")
genes_to_cluster_id_tbl[1:5,]
gene_clusterID_mean_found <- gpal_100060 %>%
  pivot_longer(-c(1:14)) %>%
  mutate(found = value!="") %>%
  select(Gene, name, found) %>%
  left_join(., genes_to_cluster_id, by="Gene") %>%
  group_by(Gene, cluster_id) %>%
  summarize(m=mean(found))

meanfound_tbl <- collect(gene_clusterID_mean_found, streaming=TRUE)
write_csv(meanfound_tbl, file = "genes_meanfound_tbl.csv")

message("Mean found (iso)")

gene_clusterID_mean_found_iso <- gpal_100060 %>%
  pivot_longer(-c(1:14)) %>%
  mutate(found = value!="") %>%
  select(Gene, name, found) %>%
  filter(name %in% iso_only) %>%
  left_join(., genes_to_cluster_id, by="Gene") %>%
  group_by(Gene, cluster_id) %>%
  summarize(m=mean(found))

meanfound_tbl_iso <- collect(gene_clusterID_mean_found_iso, streaming=TRUE)
write_csv(meanfound_tbl_iso, file = "genes_meanfound_tbl_iso.csv")

message("EggNOG ... ")

eggnog <- eggnog_annot %>% mutate(cluster_id = as.numeric(cluster_id_num)) %>%
  cast_to_int64("cluster_id")

gene_mfound_witheggnog <- gene_clusterID_mean_found %>%
  inner_join(.,
	     eggnog,
             by="cluster_id")

gene_mfound_witheggnog_iso <- gene_clusterID_mean_found_iso %>%
  inner_join(.,
	     eggnog,
             by="cluster_id")

message("Adding EggNOG (regular)")
write_csv(collect(gene_mfound_witheggnog, streaming=TRUE), "gene_mfound_witheggnog.csv")
message("Adding EggNOG (iso)")
write_csv(collect(gene_mfound_witheggnog_iso, streaming=TRUE), "gene_mfound_witheggnog_iso.csv")

gfe_iso <- read_csv("gene_mfound_witheggnog_iso.csv")
gfe <- read_csv("gene_mfound_witheggnog.csv")
prl <- read_csv("possible_ranges_labeled.csv")

core_threshold <- 0.95

gfe_iso_tbl <- gfe_iso %>% filter(m>0) %>% left_join(prl) %>% mutate(core=ifelse(m>=core_threshold, "core", "accessory")) %>% count(core, type)
gfe_all_tbl <- gfe %>% filter(m>0) %>% left_join(prl) %>% mutate(core=ifelse(m>=core_threshold, "core", "accessory")) %>% count(core, type)

write_csv(gfe_iso_tbl, "gfe_iso_tbl.csv")
write_csv(gfe_all_tbl, "gfe_all_tbl.csv")

print(gfe_iso_tbl %>% filter(type != "x") %>% group_by(core) %>% mutate(frac = n / sum(n)))
print(gfe_all_tbl %>% filter(type != "x") %>% group_by(core) %>% mutate(frac = n / sum(n)))
