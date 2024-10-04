# ROC curves for comparing correlation, EggNOG, and BLAST
library(tidyverse)
require(pROC)

results <- read_csv("SuppTable1.csv")
results_fac <- mutate(results, BLAST_Contaminate = factor(BLAST_Contaminate, levels=c("No","Yes")))
fam_ranks <- sort(unique(results$Family_max_rank))
sp_ranks <- sort(unique(results$Family_max_rank))
fdr_ranks <- sort(unique(results$Fdrs))

make_classif <- function(tbl, col, dec=TRUE) { #dec=TRUE means higher is more likely to be contaminated
  uniq_vals <- sort(unique(deframe(select(tbl, {{col}}))), decreasing=dec)
  overall <- count(tbl, BLAST_Contaminate) %>% deframe
  rstart <- c(No=0, Yes=0, val=uniq_vals[1] - (-1 * dec))
  rend <- c(No=as.numeric(overall["No"]),
            Yes=as.numeric(overall["Yes"]),
            val=uniq_vals[length(uniq_vals)] - (1 * dec))
  classif <- map(uniq_vals, \(x) {
    if (dec) {
      filt <- filter(tbl, {{col}} >= x)
    } else {
      filt <- filter(tbl, {{col}} <= x)
    }
    filt %>%
      count(BLAST_Contaminate, .drop=FALSE) %>%
      pivot_wider(names_from=BLAST_Contaminate, values_from=n) %>%
      mutate(val=x)
    })
  classif <- bind_rows(rstart, classif, rend)
  classif
}

calc_fpr_tpr <- function(tbl) {
  all_no <- max(tbl$No)
  all_yes <- max(tbl$Yes)
  tbl %>%
    mutate(tpr = map_dbl(Yes, ~ (.x / all_yes))) %>%
    mutate(fpr = map_dbl(No, ~ (.x / all_no)))
}

make_roc <- function(tbl, col, dec=TRUE) {
  make_classif(tbl, {{col}}, dec) %>% calc_fpr_tpr %>% arrange(fpr)
}

eggnog_distinct <- unique(results$Predicted_taxonomic_group)
eggnog_translate <- tribble(
  ~Predicted_taxonomic_group, ~EggNOG_contam,
  NA, 0.5,
  "Eubacteriaceae", 0,
  "unclassified Clostridiales", 0,
  "Clostridia", 0,
  "Clostridiaceae", 0,
  "Ruminococcaceae", 0,
  "Paenibacillaceae", 1,
  "Pasteurellales", 1,
  "Negativicutes", 1,
  "Desulfovibrionales", 1,
  "Bacteria", 0.5,
  "Blautia", 0,
  "unclassified Lachnospiraceae", 0,
  "Butyrivibrio", 0,
  "Erysipelotrichia", 1,
  "Oribacterium", 0,
  "Oscillospiraceae", 1
)

results_eggM_fac <- left_join(results_fac, eggnog_translate)

# Make ROCs
roc_eggnogM <- make_roc(results_eggM_fac, EggNOG_contam) %>% mutate(test="EggNOG")
roc_family <- make_roc(results_fac, Family_max_rank, dec=TRUE) %>% mutate(test="family")
roc_species <- make_roc(results_fac, Sp_rank, dec=TRUE) %>% mutate(test="species")
roc_fdrs <- make_roc(results_fac, Fdrs, dec=TRUE) %>% mutate(test="FDR")

pdf("PanSweepROCs.pdf")
bind_rows(roc_eggnogM, roc_family, roc_species, roc_fdrs) %>%
  ggplot(aes(x=fpr, y=tpr, color=test)) +
  geom_abline(slope=1, intercept=0, lty=2, color="#CCCCCC", lwd=2) +
  geom_point(size=2) +
  geom_step(lwd=1.5) +
  ylim(0, 1) +
  xlim(0, 1) +
  ylab("True positive rate (TPR)") +
  xlab("False positive rate (TPR)") +
  scale_color_brewer(type="qual", palette="Dark2") +
  theme_minimal()
dev.off()


# Calculate AUC
blast_contam_vec <- as.numeric(results_eggM_fac$BLAST_Contaminate=="Yes")
auc_eggnog <- pROC::auc(blast_contam_vec, results_eggM_fac$EggNOG_contam)
auc_family <- pROC::auc(blast_contam_vec, results_eggM_fac$Family_max_rank)
auc_species <- pROC::auc(blast_contam_vec, results_eggM_fac$Sp_rank)
auc_fdrs <- pROC::auc(blast_contam_vec, results_eggM_fac$Fdrs, direction="<")

