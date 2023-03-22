require(ggtree)
require(ggplot2)
require(ape)
require(stringr)

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

setwd("/Users/jpalmer/work/norovirus/results/mar-8/tree/capsid")
tree = read.tree('mar-8-4_gtype.nwk')
# tree = read.tree(args[1])

target_mask = str_count(tree$tip.label, paste0('\\',delim)) == 3
print(target_mask)
target_labels = tree$tip.label[target_mask]
background_labels = tree$tip.label[!target_mask]
print(target_labels)
print(background_labels)

branches <- list(Target=target_labels, Background=background_labels)
tree <- groupOTU(tree,branches)

d = fortify(tree)
d = subset(d, isTip)
tip.order = with(d, label[order(y, decreasing=F)])

g <- ggtree(tree, aes(color=group)) + geom_tiplab(align=T, linesize=0.5, size=2) + xlim(NA,2.5) +   scale_color_manual(name="Sequence Type", values=c(Background = "darkblue",Target= "red")) + guides(color = guide_legend(override.aes = list(size = 5))) + theme(legend.position = "none")


alleles <- read.table('alleles.tsv', header=T)

# order SNPs according to tip
alleles$name = factor(alleles$name, levels=tip.order)
alleles$pos = factor(alleles$pos)

bases = c("A", "C", "G", "T", "N", "X")
alleles$alt_allele = factor(alleles$alt_allele, bases)
ambiguous_bases = !(alleles$alt_allele %in% bases)
alleles$alt_allele[ambiguous_bases] = "X"

# draw panel with variants
cols <- c("blue", "red", "green3", "purple3", "lightgrey", "black")
panel.snps <- ggplot(alleles, aes(x=pos, y=name)) + 
  geom_tile(aes(fill=alt_allele), color="white") +
  ylim(tip.order) +
  theme_bw() + 
  theme(axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "top") +
  scale_fill_manual(name="Variant", values=cols, drop=FALSE) + 
  scale_x_discrete(breaks = levels(alleles$pos)[c(T, rep(F, 49))])

plot <- g + panel.snps 

ggsave('test.pdf', plot, device="pdf", height=9)
