library(tidyverse)
library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

df_1_can = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.14.lvl-3.TCGA-E9-A1RF-01A-11D-A161-05.gdc_hg38.txt")
df_1_con = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.14.lvl-3.TCGA-E9-A1RF-11A-32D-A161-05.gdc_hg38.txt")
df_2_can = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0BZ-01A-31D-A12R-05.gdc_hg38.txt" )
df_2_con = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0BZ-11A-61D-A12R-05.gdc_hg38.txt" )
df_3_can = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0HA-01A-11D-A12R-05.gdc_hg38.txt" )
df_3_con = read_tsv("Data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0HA-11A-31D-A12R-05.gdc_hg38.txt" )



colnames(df_1_can)[1:2] = c("cpg_site","beta_sample1_cancer")
colnames(df_2_can)[1:2] = c("cpg_site","beta_sample2_cancer")
colnames(df_3_can)[1:2] = c("cpg_site","beta_sample3_cancer")

colnames(df_1_con)[1:2] = c("cpg_site","beta_sample1_control")
colnames(df_2_con)[1:2] = c("cpg_site","beta_sample2_control")
colnames(df_3_con)[1:2] = c("cpg_site","beta_sample3_control")



#merge(x = df_1_can[1:2], y = df_2_can[1:2], by = "cpg_site", all.x = TRUE) 



merged_df = Reduce(function(x,y) merge(x,y,all=TRUE), list(df_1_can[1:2],df_2_can[1:2],df_3_can[1:2],df_1_con[1:2],df_2_con[1:2],df_3_con[1:2]))
merged_df_cleaned  = merged_df[complete.cases(merged_df),]

# merged_df[4,][2:4] 
# merged_df[4,][5:7] 
dim(merged_df_cleaned)

# t.test(merged_df[4,][2:4]  ,merged_df[4,][5:7])$p.value 
# 
# p_val = tibble(cpg_site = "", p_val = "")
# 
# for (i in 1:394564) {
#   site_name = merged_df_cleaned$cpg_site[i]
#   site_p = t.test(merged_df_cleaned[i,][2:4]  ,merged_df_cleaned[i,][5:7])$p.value 
#   p_val[i,] = c(site_name,site_p)
#   print(i)
# }

# p_val_adj = p.adjust(p_val, method = "BH")
# 
# p_val_adj[p_val_adj < 0.05]
# 
# hist(p_val_adj)


p_vec = vector()
for (i in 1:394564) {
  p_vec[i] = t.test(merged_df_cleaned[i,][2:4]  ,merged_df_cleaned[i,][5:7])$p.value 
  print(i)
}

cpt_result = tibble(cpg_site = merged_df_cleaned$cpg_site, p_value = p_vec)
cpt_result = merge(x = cpt_result, y = df_1_can[,c(1,3)], by = "cpg_site", all.x = TRUE)
# add start
cpt_result = merge(x = cpt_result, y = df_1_can[,c(1,4)], by = "cpg_site", all.x = TRUE)
# add end 
cpt_result = merge(x = cpt_result, y = df_1_can[,c(1,5)], by = "cpg_site", all.x = TRUE)
# bp
cpt_result$bp = cpt_result$End - cpt_result$Start

cpt_result$Chromosome = gsub("chr", "", cpt_result$Chromosome)
cpt_result$Chromosome = gsub("[a-zA-Z ]", 23, cpt_result$Chromosome)
#cpt_result$Chromosome = gsub("*", 24, cpt_result$Chromosome)
p_val_adj = p.adjust(cpt_result$p_value, method = "BH")
p_val_adj[p_val_adj < 0.05] %>% length()
hist(p_val_adj)

cpt_result = cpt_result[cpt_result$Chromosome != "*", ]


library(qqman)

manhattan(cpt_result,chr = "Chromosome", p ="p_value",bp = "bp")

cpt_result$Chromosome %>% table()
cpt_result$bp %>% table()
