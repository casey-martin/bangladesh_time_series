library(tidyverse)

sig_interact  <- read_table("./shame/downstream/sig_interactions.txt") 
col_names <- colnames(sig_interact)
sig_interact %>% 
  mutate(interact_type = if_else(r_coef > 0, "Positive", "Negative")) -> sig_interact

all_interact <- read_table("./shame/bangladesh/lag_regress_filtered90.tsv",
                           col_names = col_names)

all_interact %>% 
  filter(mean_tip_dist_i == 0 & mean_tip_dist_j ==0) %>% 
  mutate(OTUj = as.double(OTUj)) %>% 
  filter(pvalue > 0.01) -> all_interact

all_interact %>% 
  group_by(anonymized_name) %>% 
  summarise(counts = n(),
            interact_type = "None") -> all_interact_summ

sig_interact %>% 
  group_by(anonymized_name, interact_type) %>% 
  summarise(counts = n()) -> sig_interact_summ

summs <- rbind(all_interact_summ, sig_interact_summ)

#all_interact <- rbind(all_interact, sig_interact)

all_interact %>% 
  ggplot(aes(x = rand_r_coef, y = r_coef)) +
      geom_abline() +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_point(data = sig_interact, alpha = 0.5, pch = 21,
                 aes(x = rand_r_coef, y = r_coef, fill = interact_type)) +

      facet_wrap(~anonymized_name) +
      geom_label(data = all_interact_summ,
                 aes(x = -Inf, y = Inf, 
                     hjust = -0.1, vjust = 1.2,
                     label = paste("N =", counts))) +
    theme_bw(base_size = 16) +
    labs(x = "Mean Permuted Spearman Coefficient",
         y = "Empirical Spearman Coefficient") +
    guides(fill = guide_legend(title = "Association Type"))



            