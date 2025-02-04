# ROC curves for comparing correlation, EggNOG, and BLAST
library(tidyverse)
require(pROC)

results_orig <- read_csv("SuppTable1.csv")
res_blast <- results_orig %>% select(Gene_id=Num, BLAST_Contaminate, Family_max_rank, Sp_rank, Annot=`...44`)

# conditional
cc_results <- read_csv("cond_corr_table.csv")
results <- left_join(cc_results, res_blast, by="Gene_id", suffix=c("_conditional", "_overall")) %>% relocate(Gene_id)
  
results_fac <- mutate(results, BLAST_Contaminate = factor(BLAST_Contaminate, levels=c("No","Yes")))
fam_o_ranks <- sort(unique(results$Family_max_rank_overall))
sp_o_ranks <- sort(unique(results$Sp_rank_overall))
fam_c_ranks <- sort(unique(results$Family_max_rank_conditional))
sp_c_ranks <- sort(unique(results$Sp_rank_conditional))
fdr_ranks <- sort(unique(results$Fdrs))
  
make_classif <- function(tbl, col, col_answers=BLAST_Contaminate, dec=TRUE) {
  tbl <- mutate(tbl, "{{col_answers}}" := as.factor({{col_answers}}))
  uniq_vals <- sort(unique(deframe(select(tbl, {{col}}))), decreasing=dec)
  overall <- count(tbl, {{col_answers}}) %>% deframe
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
      count({{col_answers}}, .drop=FALSE) %>%
      pivot_wider(names_from={{col_answers}}, values_from=n) %>%
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

make_roc <- function(tbl, col, col_correct=BLAST_Contaminate, dec=TRUE) {
  make_classif(tbl, {{col}}, {{col_correct}}, dec) %>% calc_fpr_tpr %>% arrange(fpr)
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
roc_family_o <- make_roc(results_fac, Family_max_rank_overall, dec=TRUE) %>% mutate(test="familyO")
roc_species_o <- make_roc(results_fac, Sp_rank_overall, dec=TRUE) %>% mutate(test="speciesO")
roc_family_c <- make_roc(results_fac, Family_max_rank_conditional, dec=TRUE) %>% mutate(test="familyC")
roc_species_c <- make_roc(results_fac, Sp_rank_conditional, dec=TRUE) %>% mutate(test="speciesC")
roc_fdrs <- make_roc(results_fac, Fdrs, dec=TRUE) %>% mutate(test="FDR")


roc_colors <- c(EggNOG='#1b9e77',
                familyC='#d95f02',
                familyO='#7570b3',
                FDR='#e7298a',
                speciesC='#66a61e',
                speciesO='#e6ab02')
pdf("PanSweepROCs.pdf")
bind_rows(roc_eggnogM, roc_family_o, roc_species_o, roc_family_c, roc_species_c, roc_fdrs) %>%
  mutate(test = fct_relevel(as.factor(test), c("EggNOG","familyO","familyC","speciesO","speciesC","FDR"))) %>%
  ggplot(aes(x=fpr, y=tpr, color=test)) +
  geom_abline(slope=1, intercept=0, lty=2, color="#CCCCCC", lwd=2) +
  geom_point(size=2) +
  geom_step(lwd=1.5) +
  ylim(0, 1) +
  xlim(0, 1) +
  ylab("True positive rate (TPR)") +
  xlab("False positive rate (TPR)") +
  scale_color_manual(values=roc_colors) +
  theme_minimal()
dev.off()


# Calculate AUC
blast_contam_vec <- as.numeric(results_eggM_fac$BLAST_Contaminate=="Yes")
auc_eggnog <- pROC::auc(blast_contam_vec, results_eggM_fac$EggNOG_contam)
auc_family_c <- pROC::auc(blast_contam_vec, results_eggM_fac$Family_max_rank_conditional)
auc_species_c <- pROC::auc(blast_contam_vec, results_eggM_fac$Sp_rank_conditional)
auc_family_o <- pROC::auc(blast_contam_vec, results_eggM_fac$Family_max_rank_overall)
auc_species_o <- pROC::auc(blast_contam_vec, results_eggM_fac$Sp_rank_overall)

auc_egg_fam_either <- pROC::auc(blast_contam_vec, results_eggM_fac$Family_max_rank_overall)

auc_fdrs <- pROC::auc(blast_contam_vec, results_eggM_fac$Fdrs, direction="<")

###

# Comparison to BT
bt <- read_tsv("eggNOG_Genes_to_BLAST.tsv") %>%
  rename(blast_consistent = `BLAST – Bt`) %>% # note, opposite of above
  mutate(blast_contam = blast_consistent=="No") %>%
  mutate(blast_contam_str = ifelse(blast_consistent=="Yes", "No", "Yes"))

# > bt$Predicted_taxonomic_group %>% unique
# [1] NA                     "Bacteroidaceae"       "Clostridiaceae"       "Sutterellaceae"      
# [5] "Eubacteriaceae"       "Oscillospiraceae"     "Streptococcus oralis" "Ruminococcaceae"     
# [9] "Negativicutes"        "Dorea" 

bt_egg_classif <- tribble(~Predicted_taxonomic_group, ~EggNOG_contam,
                          NA, 0.5,
                          "Bacteroidaceae", 0,
                          "Clostridiaceae", 1,
                          "Eubacteriaceae", 1,
                          "Oscillospiraceae", 1,
                          "Streptococcus oralis", 1,
                          "Sutterellaceae", 1,
                          "Ruminococcaceae", 1,
                          "Negativicutes", 1,
                          "Dorea", 1)

# 0.7174
bt_auc_eggnog <- bt %>%
  select(Predicted_taxonomic_group, Family_max_rank, blast_consistent) %>%
  (\(.) pROC::auc(.$blast_consistent,
                  1 * grepl("Bacteroid", .$Predicted_taxonomic_group)))

# 0.9099
bt_auc_fmr <- bt %>%
  select(Predicted_taxonomic_group, Family_max_rank, blast_contam) %>%
  (\(.) pROC::auc(1 * .$blast_contam, .$Family_max_rank))

bt_roc_eggnog <- make_roc(left_join(bt, bt_egg_classif),
                           EggNOG_contam,
                           blast_contam_str) %>% mutate(test="EggNOG")
bt_roc_fmr <- make_roc(bt,
                           Family_max_rank,
                           blast_contam_str) %>% mutate(test="familyO")
bt_roc_sp <- make_roc(bt,
                       Sp_rank,
                       blast_contam_str) %>% mutate(test="speciesO")

pdf("btheta-roc.pdf")
bind_rows(bt_roc_eggnog, bt_roc_fmr, bt_roc_sp) %>%
  mutate(test = fct_relevel(as.factor(test), c("EggNOG","familyO", "speciesO"))) %>%
  ggplot(aes(x=fpr, y=tpr, color=test)) +
  geom_abline(slope=1, intercept=0, lty=2, color="#CCCCCC", lwd=2) +
  geom_point(size=2) +
  geom_step(lwd=1.5) +
  ylim(0, 1) +
  xlim(0, 1) +
  ylab("True positive rate (TPR)") +
  xlab("False positive rate (TPR)") +
  scale_color_manual(values=roc_colors) +
  theme_minimal()
dev.off()
