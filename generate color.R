library(RColorBrewer)
brewer.pal(6,"Set1")

library(grDevices)
#Heat map color 
heat.colors(23, alpha=1)
#color_list = dput(as.character(heat.colors(23, alpha=1)))
cat(paste(shQuote(heat.colors(23, alpha=1), type="cmd"), collapse=", "))
#rainbow palette 
n = 23
rainbow_list = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
color_list = dput(rainbow_list)



library(colorspace)
rainbow_hcl(23)


# diverge_hcl 
# diverge_hsl
# terrain_hcl
# sequential_hcl
# rainbow_hcl
#However, all palettes are fully customizable:
#diverge_hcl(7, h = c(246, 40), c = 96, l = c(65, 90))

# color choosing tool 
pal <- choose_palette()


library(colorRamps)
colorRamps::primary.colors(23)
