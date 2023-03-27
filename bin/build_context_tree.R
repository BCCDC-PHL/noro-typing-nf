library(ggtree)
library(ggplot2)
library(ape)
library(stringr)
library(grid)

args = commandArgs(trailingOnly = T)

if (length(args) != 2 && length(args) != 3){
    print("USAGE: Rscript build_context_tree.R INPUT_TREE OUTPUT_PDF [HEADER_DELIM]")
    quit()
}

if (length(args) == 2){
  delim = '|'
}else{
  delim = args[3]
}

delim = '|'
# 1: input tree
# 2: metadata file indicating target / background 
# 3: output PDF plot 

# setwd("/Users/jpalmer/work/norovirus/references/G9-test")

tree = read.tree(args[1])
# tree = read.tree('~/work/norovirus/results-jan-31-9/phylo/capsid/tree/jan-31-9-gtype.nwk')

target_mask = str_count(tree$tip.label, '\\|') == 3
print(target_mask)
target_labels = tree$tip.label[target_mask]
background_labels = tree$tip.label[!target_mask]
print(target_labels)
print(background_labels)

branches <- list(Target=target_labels, Background=background_labels)
tree <- groupOTU(tree,branches)

print("ggtree")
p <- ggtree(tree, aes(color=group)) +
  geom_tiplab(align=T, linesize=0.5) + xlim(NA,2.5) +   
  scale_color_manual(name="Sequence Type", values=c(Background = "darkblue",Target= "red")) + 
  guides(color = guide_legend(override.aes = list(size = 5))) 

print("ggsave")
# pdf(file = args[2], width = 12, height = 8)
# plot(p)
ggsave(args[2], device="pdf")
# dev.off()


#tree$tip.label = sapply(tree$tip.label, function(x){ if (grepl("\\|", x)){ paste(strsplit(x, '\\|')[[1]][2:3], collapse="_")} else {x}})

