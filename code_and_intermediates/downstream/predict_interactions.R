library(data.table)
library(tidyverse)

mydata <- fread('/home/operon/Documents/lag_regress/bangladesh/lag_regress_filtered90.tsv')
mygenomes <- fread('/home/operon/Documents/lag_regress/OTU_KO_table.txt')

mydata$V3 = 1 - mydata$V3

mydata %>% 
    filter(!str_detect(V1, 'node'), !str_detect(V2, 'node'), V1 != V2, V3 < 0.01) -> mydata



mydata$sig <- abs(mydata$V4) > abs(mydata$V5)

iGenomes <- t(sapply(mydata$V1, function(x) mygenomes[mygenomes$OTU == x,]))
jGenomes <- t(sapply(mydata$V2, function(x) mygenomes[mygenomes$OTU == x,]))
