library(methyAnalysis)
library(tidyverse)
library(CpGassoc)

df=  merged_df_cleaned
rownames(df) = merged_df_cleaned$cpg_site
df = df[-1]



df_result = cpg.assoc(df,indep = c(0,0,0,1,1,1) )
df_result$coeff
plot(df_result)



data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
results




library(gtools)
#dev.off()
df$beta_sample1_cancer %>% logit() %>% hist()
df$beta_sample1_cancer%>% hist()
