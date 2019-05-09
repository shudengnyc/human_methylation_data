# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 23 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 23) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}
# tranfor data 
plot_df = cpt_result
colnames(plot_df) = c("SNP","P","CHR","BP","end","BP1")
plot_df$CHR = as.numeric(plot_df$CHR)
# get FDR cut off line 
sort(plot_df$P)[25]
# Variables ====
mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette

mysnps <- c("rs11801961","rs116558464","rs61703161") # snps to highlight
sig = 0.05/394564 # significant threshold line
sugg = sort(plot_df$P)[25] # suggestive threshold line

# Define Function ====
# see below
# Run Function 
gg.manhattan(plot_df, threshold=1e-6, hlight=mysnps, col=mypalette, ylims=c(0,10), title="My Manhattan Plot")
gg.manhattan(plot_df, threshold=1e-6, hlight=mysnps, col=sample(color_list), ylims=c(0,10), title="Manhattan Plot")


# use qqman to plot 
#manhattan(plot_df)
