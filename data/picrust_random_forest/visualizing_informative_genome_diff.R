library(data.table)
library(tidyverse)

mydiffs <- fread('~/Documents/lag_regress/genome_diffs.txt', data.table = F)
myinteractions <- fread('~/Documents/lag_regress/sig_interactions.txt', data.table = F)
featureImportances <- fread('~/Documents/lag_regress/feature_importances.csv', data.table = F)
colnames(featureImportances) <- c('KO', 'importance')

mytopKos <- fread('~/Documents/lag_regress/predictive_kos.csv', data.table = F)
colnames(mytopKos) <- c('V1', 'KO', 'ko_name')

mytopKos <- merge(mytopKos, featureImportances)

mytopKos <- mytopKos[rev(order(mytopKos$importance)),]


koDict <- sapply(mytopKos[,3], function(x) str_wrap(x, 20))

names(koDict) <- mytopKos[,1]

topKoDiffs <- mydiffs[,colnames(mydiffs) %in% mytopKos[,1]]

myinteractions <- cbind(myinteractions[c('OTUi', 'OTUj', 'r_coef', 'i_j_dist')], topKoDiffs)

myinteractions %>% 
  group_by(OTUi, OTUj, r_coef, i_j_dist) %>% 
  gather(-OTUi, -OTUj, -r_coef, -i_j_dist, 
         key = KO, value = Difference) -> long_interactions


long_interactions <- merge(long_interactions, featureImportances)
long_interactions$ko_name <- koDict[long_interactions$KO]

long_interactions %>% 
  ggplot(aes(x = Difference, y = i_j_dist)) +
    geom_point() + 
    geom_boxplot(alpha=0) +
    facet_wrap(~ko_name)

long_interactions %>% 
  filter(importance > 0.0033) %>% 
  ggplot(aes(x = if_else(r_coef > 0, 'Positive', 'Negative'),
                         fill = Difference)) +
    theme_bw(base_size = 20) +
    geom_bar() +
    xlab('Interaction Type') +
    facet_wrap(~ko_name, nrow = 1) +
    theme(legend.position = 'bottom')

long_interactions %>% 
  group_by(OTUi, OTUj, r_coef, i_j_dist, importance, ko_name) %>% 
  spread(KO, Difference)

head(long_interactions, n=20)
