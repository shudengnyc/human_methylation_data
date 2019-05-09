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

 

# check data distritbuion 
histo_raw = merged_df_cleaned$beta_sample1_cancer %>% hist(main = "Raw Beta Value")
# try log2 transformation 
merged_df_cleaned$beta_sample1_cancer %>% log2() %>% hist()

# use log2ratio transformation 
foo = sapply(merged_df_cleaned$beta_sample1_cancer, function(x) log(x[1]/(1-x[1])))
histo_transfomed = foo %>% hist(main = "Log 2 Ratio Transformed Beta Value")

# plot 
par(mfrow=c(1,2))  
merged_df_cleaned$beta_sample1_cancer %>% hist(main = "Raw Beta Value")
histo_transfomed = foo %>% hist(main = "Log 2 Ratio Transformed Beta Value")
# #### some other automatic Tranformation 
# 
library(rcompanion)
# merged_df_cleaned$beta_sample1_cancer %>% hist()
# plotNormalHistogram(merged_df_cleaned$beta_sample1_cancer)

plotNormalHistogram(merged_df_cleaned$beta_sample1_cancer,linecol = "claer",main = "Raw Beta Value",xlab = "")
plotNormalHistogram(foo,main = "Log 2 Ratio Transformed Beta Value",xlab = "")

# qq plot
qqnorm(foo, ylab = "Sample Quantiles for beta")
qqline(foo, col = "red")
# qqnorm(merged_df_cleaned$beta_sample1_cancer, ylab = "Sample Quantiles for Turbidity")
# qqline(merged_df_cleaned$beta_sample1_cancer, col = "red")
# ## box - cox tranfrom 
# 
# library(MASS)
# 
# Box = boxcox(merged_df_cleaned$beta_sample1_cancer ~ 1,              # Transform Turbidity as a single vector
#              lambda = seq(-6,6,0.1)      # Try values -6 to 6 by 0.1
# )
# Cox = data.frame(Box$x, Box$y) 
# Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
# Cox2[1,] 
# lambda = Cox2[1, "Box.x"]
# T_box = (merged_df_cleaned$beta_sample1_cancer ^ lambda - 1)/lambda 
# plotNormalHistogram(T_box)
# 
# #merged_df_cleaned$beta_sample1_cancer %>% log2() %>% hist() 
# library(EnvStats)
# boxcox.list <- boxcox(merged_df_cleaned$beta_sample1_cancer)
# plot(boxcox.list)
# plot(boxcox.list, plot.type = "Q-Q Plots", same.window = FALSE) 
# #### tranformation 


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

# conduct t-test 
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

# Manhattan Plot
# library(qqman)
# manhattan(cpt_result,chr = "Chromosome", p ="p_value",bp = "bp")

cpt_result$Chromosome %>% table()
cpt_result$bp %>% table()
