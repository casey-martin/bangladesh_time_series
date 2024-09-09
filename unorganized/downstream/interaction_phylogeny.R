library(cowplot)
library(data.table)
library(ggsignif)
library(tidyverse)

alldata <- fread('./bangladesh/lag_regress_filtered90.tsv')
colnames(alldata) <- c("OTUi", "OTUj", "pvalue", "r_coef",
                               "rand_r_coef", "skew", "kurt", "sd",
                               "mean_tip_dist_i", "mean_tip_dist_j",
                               "i_j_dist", "anonymized_name")

alldata %>% 
  filter(mean_tip_dist_i == 0 & mean_tip_dist_j == 0 & i_j_dist != 0) -> allinteractions

allinteractions$skew <- as.numeric(allinteractions$skew)

allinteractions %>% 
  arrange(desc(pvalue)) %>% 
  mutate(pvalue = 1 - pvalue,
         pvalue_fdr = p.adjust(pvalue, method = 'fdr'),
         significant = pvalue_fdr <= 0.1,
         interaction_type = case_when(
           pvalue_fdr > 0.1 ~ 'None',
           pvalue_fdr <= 0.1 & r_coef > 0 ~ 'Positive',
           pvalue_fdr <= 0.1 & r_coef < 0 ~ 'Negative')
         ) -> allinteractions

allinteractions %>% 
  ggplot(aes(x = rand_r_coef, y = sd, color = OTUi == OTUj)) +
    geom_point() +
    facet_wrap(~interaction_type) +
    theme_bw()

allinteractions %>% 
  ggplot(aes(x = rand_r_coef)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0) +
    facet_wrap(~interaction_type, scales = "free_y") +
    theme_bw()

allinteractions %>% 
  ggplot(aes(x = interaction_type, y = i_j_dist)) +
    geom_jitter(alpha = 0.1, size=0.9, 
                height = 0, width = 0.2) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_signif(comparisons = list(c(1,2),
                                   c(1,3),
                                   c(2,3)),
                step_increase = 0.07,
                tip_length = 0) +
    theme_bw(base_size = 16) +
    xlab('Interaction Type') +
    ylab('Phylogenetic Distance\nBetween Microbes') -> phylogeny_interaction
    

allinteractions %>% 
  ggplot(aes(x = rand_r_coef, y = r_coef)) +
    geom_point(aes(alpha = significant, fill=anonymized_name), pch=21) +
    xlab('Mean Permuted R Coefficient') +
    ylab('Observed R Coefficient') +
    theme_bw(base_size = 16) +
    labs(alpha = 'Significant', fill = 'ID') +
    theme(legend.position = 'bottom')-> r_coef_corr
  
plot_grid(r_coef_corr, phylogeny_interaction, 
          labels = c('A.', 'B.'),
          label_size = 22,
          rel_widths = c(1, 0.5))


  
