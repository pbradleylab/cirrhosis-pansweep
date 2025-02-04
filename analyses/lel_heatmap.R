library(tidyverse)
library(ComplexHeatmap)
library(arrow)
library(PanSweep)
library(cluster)
# Load data

cirr_md <- read_tsv("../data/manual/cirr_md.tsv")

species_r <- read_tsv("../data/processed/midas_merged/species/species_relative_abundance.tsv") %>%
  merge_columns_tbl(., cirr_md, mean) %>%
  PanSweep:::tbl_to_mtx()

lel_copynum <- arrow::read_tsv_arrow("../data/processed/midas_snv_merge/genes/100060/100060.genes_copynum.tsv.lz4")
lel_copynum_merged <- merge_columns_tbl(lel_copynum, cirr_md, mean)
lel_copynum_merged_mtx <- PanSweep:::tbl_to_mtx(lel_copynum_merged)

lel_presabs_merged <- merge_columns_tbl(
  arrow::read_tsv_arrow("../data/processed/midas_snv_merge/genes/100060/100060.genes_presabs.tsv.lz4"),
  cirr_md,
  max)
lel_presabs_merged_mtx <- PanSweep:::tbl_to_mtx(lel_presabs_merged)

lel_sig_genes <- scan("../data/manual/lel_sig_genes.txt", "character")
cirr_pres_cases <- cirr_md %>% filter(env=="case") %>% select(subject) %>% filter(subject %in% colnames(lel_presabs_merged_mtx)) %>% deframe() %>% unique()
cirr_pres_controls <- cirr_md %>% filter(env=="control") %>% select(subject) %>% filter(subject %in% colnames(lel_presabs_merged_mtx)) %>% deframe() %>% unique()

lel_sig_presabs_tbl <- sapply(lel_sig_genes %>% unique(), function(x) {
  x_case <- lel_presabs_merged_mtx[x, cirr_pres_cases]
  x_ctrl <- lel_presabs_merged_mtx[x, cirr_pres_controls]
  c(sum_case=sum(x_case), l_case = length(x_case), avg_case=mean(x_case), sum_ctrl=sum(x_ctrl), l_ctrl = length(x_ctrl), avg_ctrl=mean(x_ctrl))
}) %>% t() %>% as_tibble(rownames="Gene_id")

#FliQ
lel_sig_presabs_tbl %>% filter(Gene_id == "UHGG152466_01649") %>% print()
#FliW UHGG148794_01632
lel_sig_presabs_tbl %>% filter(Gene_id == "UHGG148794_01632") %>% print()
#GGDEF
lel_sig_presabs_tbl %>% filter(Gene_id == "UHGG228006_01743") %>% print()
#diguan
lel_sig_presabs_tbl %>% filter(Gene_id == "UHGG239171_01612") %>% print()

lel_eggnog <- read_csv("../data/manual/lel_eggnog.csv.gz")
lel_flag_annot <- read_csv("../data/manual/lel_flag_annot.csv")
additional_genes <- c(
  "UHGG228006_01743",
  "UHGG239171_01612"
)

# add relevant "non-flagellar" genes to annotation, and fix missing short gene names
lel_ccf_genes_tbl <- tibble(gene_id_string = c(lel_flag_annot$gene_id_string,
                                               additional_genes)) %>%
  mutate(gene_id = as.numeric(gsub("UHGG(.*)_(.*)", "\\1\\2", gene_id_string))) %>%
  filter(gene_id_string %in% rownames(lel_copynum_merged_mtx)) %>%
  left_join(., lel_eggnog)

lel_ccf_descs <- lel_ccf_genes_tbl %>%
  mutate(pid = map2_chr(
    Predicted_protein_name,
    eggNOG_free_text_description,
    function(x,y) {
      if (!is.na(x)) { return(x) }
      if (grepl("FliO", y)) return("fliO")
      if (grepl("Flp1-like", y)) return("flp1-like")
      if (grepl("Forms a capping structure",y)) return("fliD")
      if (grepl("FliK", y)) return("fliK")
      if (grepl("VirB11", y)) return("t4ss_virb11")
      if (grepl("Controls the rotational direction of flagella during",
                y)) return("FliL")
      if (grepl("EAL.*GGDEF", y)) return("ealggdef")
      if (grepl("diguanylate cyclase", y)) return("diguancyc")
      return(NA)
    })) %>%
  mutate(ext_pid = paste(pid, gene_id_string, sep="_"))
additional_genes_ext <- lel_ccf_descs %>% filter(gene_id_string %in% additional_genes) %>% .$ext_pid

# Subset and convert row names, and filter to keep only genes that are "present" at least 25% of the time
lel_ccf2 <- lel_copynum_merged_mtx[lel_ccf_descs$gene_id_string,] %>%
  .[union(additional_genes, names(which(rowMeans(. >= 0.5) >= 0.05))),]
rownames(lel_ccf2) <- (lel_ccf_descs %>%
                         select(gene_id_string, ext_pid) %>%
                         deframe)[rownames(lel_ccf2)]

lel_ccf2_casectrl <- (deframe(
  cirr_md[, c("subject", "env")])
)[colnames(lel_ccf2)]
cirr_cases <- cirr_md$subject[cirr_md$env == "case"]
cirr_controls <- cirr_md$subject[cirr_md$env == "control"]
lel_ccf2_signif <- (rownames(lel_ccf2) %>%
                      gsub(".*(UHGG.*)$","\\1",.)) %in% lel_sig_genes

# bootstrap average of t-stat
lel_enrichment_alt <- apply(lel_ccf2, 1, function(x) {
  boot_t_est <- map_dbl(1:500, ~ {
    boot <- sample(names(x), replace=TRUE)
    boot_cases <- boot[boot %in% cirr_cases]
    boot_controls <- boot[boot %in% cirr_controls]
    t.test(x[boot_cases], x[boot_controls])$statistic
  })
  mean(boot_t_est)
})

# This is very similar to the t-stat results
lel_enrichment_nonpara <- apply(lel_ccf2, 1, function(x) {
  w <- wilcox.test(x[cirr_cases], x[cirr_controls], conf.int=TRUE)
  -log10(w$p.value) * sign(w$estimate)
})

# plot results
pdf(width=11,height=8.5,"lel_flagellar.pdf")
lel_ccf2[rownames(lel_ccf2),] %>%
  Heatmap(col=circlize::colorRamp2(breaks=c(0,1,2,3),
                                   colors=c("dodgerblue4","red3",
                                            "yellow4","yellow")),
          clustering_distance_rows="pearson",
          clustering_distance_columns="euclidean",
          column_split = lel_ccf2_casectrl,
          right_annotation = rowAnnotation(
            df=data.frame(significant=lel_ccf2_signif,
                          enrichment=lel_enrichment[rownames(lel_ccf2)]),
            col=list(significant=c("FALSE"="white","TRUE"="black"),
                     enrichment=circlize::colorRamp2(breaks=c(-4,0,4),
                                                     colors=c("dodgerblue",
                                                              "white",
                                                              "indianred")))),
          row_names_gp=gpar(fontsize=5, fontface=2),
          column_names_gp=gpar(fontsize=8, fontface=2))
dev.off()

# cluster genes
corr_dist <- sqrt(0.5*(1-cor(t(lel_ccf2), method='pearson'))) # between 0 and 1
cluster_cors <- lapply(2:25, \(x) pam(corr_dist, diss=TRUE, k=x))
cluster_sil <- map_dbl(cluster_cors, \(x) x$silinfo$avg.width)
best_clustering <- cluster_cors[[which(rank(1-cluster_sil)==1)]]
best_n <- length(best_clustering$medoids)
#plot(cluster_sil)
gene_lists <- best_clustering$clustering %>% enframe(name="gene_id", value="cluster") %>%
  nest(data=gene_id) %>%
  mutate(cluster_name=paste0("C", cluster)) %>%
  mutate(n=map_dbl(data, nrow)) %>%
  filter(n >= 3) %>%
  arrange(cluster) %>%
  select(cluster_name, data) %>%
  deframe

# get and test eigengenes
eigengenes <- lapply(gene_lists, \(x) svd(lel_ccf2[x$gene_id, ])$v[,1])
eig_mtx <- lapply(gene_lists, \(x) svd(lel_ccf2[x$gene_id, ])$v[,1] %>% setNames(., colnames(lel_ccf2))) %>% bind_rows(.id = "gene_id") %>% PanSweep:::tbl_to_mtx()
## note: copy numbers are always strictly positive, so we need to rotate these to get the right interpretation
eig_mtx_pos <- apply(eig_mtx, 1, \(x) {
  avg_sign <- mean(sign(x))
  if (avg_sign==0) return(x)
  return(x * avg_sign)
}) %>% t
eig_wilcox <- apply(eig_mtx_pos, 1, \(x) wilcox.test(x[cirr_cases], x[cirr_controls], conf.int=TRUE))
eig_wilcox_p <- map_dbl(eig_wilcox, ~ .$p.value)
eig_wilcox_est <- map_dbl(eig_wilcox, ~ .$estimate)
eig_wilcox_padj <- p.adjust(eig_wilcox_p, 'BH')

# test whether clusters significant
gene_cluster_signif <- eig_wilcox_padj %>%
  enframe(value="p") %>%
  left_join(enframe(eig_wilcox_est, value="estimate"), by="name") %>%
  left_join(enframe(gene_lists, value="gene_list"), by="name") %>%
  unnest(gene_list) %>%
  mutate(cluster_signif = map2_dbl(p, estimate, ~ ifelse(.x <=0.05, sign(.y), 0))) %>%
  select(gene_id, cluster_signif) %>%
  deframe %>%
  .[rownames(lel_ccf2)] %>%
  setNames(rownames(lel_ccf2))

# plot results
pdf(width=11,height=8,"lel_flagellar_clusters.pdf")
(lel_ccf2[rownames(lel_ccf2),]) %>%
  Heatmap(col=circlize::colorRamp2(breaks=c(0,1,2,3),
                                   colors=c("dodgerblue4","red3",
                                            "yellow4","yellow")),
          clustering_distance_rows="pearson",
          clustering_distance_columns="euclidean",
          column_split = lel_ccf2_casectrl,
          row_split = best_clustering$clustering[rownames(lel_ccf2)],
          row_gap = unit(0.25, "mm"),
          left_annotation = rowAnnotation(
            na_col="#DDDDDD",
            df=data.frame(cluster_signif=gene_cluster_signif),
            col=list(cluster_signif=c(`-1`="dodgerblue", `0`="white", `1`="indianred"))),
          right_annotation = rowAnnotation(
            df=data.frame(pansweep_signif=lel_ccf2_signif,
                          enrichment=lel_enrichment_alt[rownames(lel_ccf2)]),
            col=list(pansweep_signif=c("FALSE"="white","TRUE"="black"),
                     enrichment=circlize::colorRamp2(breaks=c(-4,0,4),
                                                     colors=c("dodgerblue",
                                                              "white",
                                                              "indianred")))),
          row_names_gp=gpar(fontsize=5, fontface=2),
          row_title_gp=gpar(fontsize=7, fontface=2),
          column_names_gp=gpar(fontsize=5, fontface=2))
dev.off()

