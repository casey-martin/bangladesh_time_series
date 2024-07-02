library(tidyverse)
library(latex2exp)
library(grid)
library(gridExtra)

mybiom <- read_tsv("./data/qiime_outputs/dada2_otu_table_w_tax_no_pynast_failures_rare7500.tsv", skip = 1)
mymeta <- read_csv("./data/qiime_outputs/full_meta_w_times.csv")
mymeta %>% 
  select(`#SampleID`, Time, ANONYMIZED_NAME) -> min_meta

theme_set(theme_bw(base_size = 16))

mybiom %>% 
  gather(-`#OTU ID`, key = "#SampleID", value = "abundance") %>% 
  left_join(min_meta) -> long_biom

long_biom %>% 
  group_by(ANONYMIZED_NAME, Time, `#OTU ID`) %>% 
  summarise(abundance = mean(abundance)) %>% 
  arrange(ANONYMIZED_NAME, `#OTU ID`, Time) %>% 
  group_by(ANONYMIZED_NAME) %>% 
  mutate(n_time = length(unique(Time))) %>% 
  group_by(`#OTU ID`, ANONYMIZED_NAME) %>% 
  mutate(otu_prev = (1 - sum(abundance == 0) / n_time)) %>% 
  filter(otu_prev >= 0.9) %>% 
  mutate(rel_abund = abundance/7500) %>% 
  group_by(ANONYMIZED_NAME, `#OTU ID`) %>% 
  mutate(first_lag = lead(rel_abund) - rel_abund) -> long_biom


long_biom %>% 
  select(ANONYMIZED_NAME, `#OTU ID`, rel_abund) %>% 
  group_by(ANONYMIZED_NAME, `#OTU ID`) %>%
  slice(sample(1:n())) %>% 
  mutate(first_lag = lead(rel_abund) - rel_abund) -> rand_long_biom
  

rand_long_biom$Time <- long_biom$Time
rand_long_biom$data_type <- "Randomized"
long_biom$data_type <- "Observed"

long_biom <- rbind(long_biom, rand_long_biom)

#long_biom$rand_rel_abund <- rand_long_biom$rel_abund
#long_biom$rand_first_lag <- rand_long_biom$first_lag


long_biom %>% 
  ggplot(aes(x = Time, y = rel_abund)) +
  geom_vline(xintercept = 0, color="red") +
  geom_line(aes(group=`#OTU ID`), alpha=0.25) +
  facet_grid(data_type~ANONYMIZED_NAME) +
  scale_y_log10() +
  xlab("") +
  ylab("") -> plot_1

long_biom %>% 
  ggplot(aes(x = Time, y = first_lag)) +
    geom_vline(xintercept = 0, color="red") +
    geom_line(aes(group=`#OTU ID`), alpha=0.25) +
    facet_grid(data_type~ANONYMIZED_NAME) +
  ylab("") +
  xlab("") -> plot_2

plot_grid(plot_1 +
            ylab("Log10 Relative Abundance"),
          plot_2 +
            ylab("1st Derivative\n(Lag Difference)"),
          ncol = 1,
          align="hv")


long_biom %>% 
  ggplot(aes(x = first_lag)) +
    geom_histogram() +
    scale_y_log10() +
    facet_wrap(~ANONYMIZED_NAME)


long_biom %>% 
  ggplot(aes(x = rel_abund, y = first_lag, group = `#OTU ID`)) +
    geom_point(size=0.1) +
    geom_smooth(method = 'lm', se=FALSE, alpha=0.1, size=0.2) +
    geom_hline(yintercept=0) +
    facet_grid(data_type~ANONYMIZED_NAME, scales = "free")


long_biom %>% 
  filter(rel_abund > 0, data_type == "Observed") %>% 
  select(ANONYMIZED_NAME, `#OTU ID`) %>% 
  distinct() %>% 
  group_by(`#OTU ID`) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))



  filter(`#OTU ID` == "125624") %>% 
  ggplot(aes(x = rel_abund, y = first_lag, group = `#OTU ID`)) +
  geom_point(size=0.1, color = "") +
  geom_smooth(method = 'lm', se=FALSE, alpha=0.1, size=0.2) +
  geom_hline(yintercept=0) +
  facet_grid(data_type~ANONYMIZED_NAME, scales = "free")



