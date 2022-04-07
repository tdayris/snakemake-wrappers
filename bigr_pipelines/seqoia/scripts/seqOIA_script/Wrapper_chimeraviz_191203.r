
library(chimeraviz) 

# https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/


args <- commandArgs(trailingOnly=TRUE)
#print(args)

input = args[1]
output= args[2]

fusions <- import_starfusion(input, "hg38")

svg(output)
plot_circle(fusions)
dev.off()


#print(input)
#print(output)

