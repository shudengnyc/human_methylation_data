library(tidyverse)

df_1_can = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.14.lvl-3.TCGA-E9-A1RF-01A-11D-A161-05.gdc_hg38.txt")
df_1_con = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.14.lvl-3.TCGA-E9-A1RF-11A-32D-A161-05.gdc_hg38.txt")
df_2_can = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0BZ-01A-31D-A12R-05.gdc_hg38.txt" )
df_2_con = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0BZ-11A-61D-A12R-05.gdc_hg38.txt" )
df_3_can = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0HA-01A-11D-A12R-05.gdc_hg38.txt" )
df_3_con = read_tsv("data/jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0HA-11A-31D-A12R-05.gdc_hg38.txt" )



colnames(df_1_can)[1:2] = c("cpg_site","beta_sample1_cancer")
colnames(df_2_can)[1:2] = c("cpg_site","beta_sample2_cancer")
colnames(df_3_can)[1:2] = c("cpg_site","beta_sample3_cancer")

colnames(df_1_con)[1:2] = c("cpg_site","beta_sample1_control")
colnames(df_2_con)[1:2] = c("cpg_site","beta_sample2_control")
colnames(df_3_con)[1:2] = c("cpg_site","beta_sample3_control")



merge(x = df_1_can[1:2], y = df_2_can[1:2], by = "cpg_site", all.x = TRUE) 



merged_df = Reduce(function(x,y) merge(x,y,all=TRUE), list(df_1_can[1:2],df_2_can[1:2],df_3_can[1:2],df_1_con[1:2],df_2_con[1:2],df_3_con[1:2]))
merged_df_cleaned  = merged_df[complete.cases(merged_df),]

# merged_df[4,][2:4] 
# merged_df[4,][5:7] 


# t.test(merged_df[4,][2:4]  ,merged_df[4,][5:7])$p.value 
# 
p_val = vector()

for (i in 1:400000) {
  p_val[i] = t.test(merged_df_cleaned[i,][2:4]  ,merged_df_cleaned[i,][5:7])$p.value 
}

p_val_adj = p.adjust(p_val, method = "BH")

p_val_adj[p_val_adj < 0.05]

hist(p_val_adj)
