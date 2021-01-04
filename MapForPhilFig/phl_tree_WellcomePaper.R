### Trees

#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
# install_github("YuLab-SMU/ggtree")
library(ape)
library(stringr)
library(data.table)
library(tidytree)
library(phylobase)
library(ggrepel)
library(ggtree)
library(treeio)
library(devtools)
library(dplyr)
library(rgdal)
library(ggplot2)
library(broom)
library(rgeos)
library(prettymapr)
library(viridis)
library(RColorBrewer)
library(wesanderson)

## simple tree visualisation
rawtree = read.tree("AL_WG_ASIAN_SEA4_og_rerooted.tree")
tree=drop.tip(rawtree,grep("GU358653|GU647092",rawtree$tip.label))

# extract tip info
# meta=tree$tip.label
# meta2=as.data.table(tstrsplit(meta, "/"))
# first column must be taxa labels (or node) to allow mapping back to tree
# meta2=cbind(tree$tip.label,meta2)
# names(meta2)=c("taxa","aln", "clade", "member","source","sample_id")

#extracted tip info combined with additional data:
meta = read.csv("AL_WG_ASIAN_SEA4_metadata.csv")
meta$label2 = paste("Reg", meta$region,"|",meta$year,"|",meta$sample_id.1, sep="" )
# meta$region
region_levels = sort(unique(meta$region))
#removing outgroup labels (don't fit in plot)
#meta$label2[3]=""
#meta$label2[4]=""
meta$label2=gsub("Regunknown","unknown", meta$label2)
meta$label2=gsub("Regoutgroup","outgroup", meta$label2)
meta$label2[meta$label2=="outgroup|na|GU647092"]=""
#colours
##ggplot colours:
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# n = 6
# reg_colours = c(gg_color_hue(n),"transparent")

##wes colours:
##continuous if more than 5
pal=wes_palette(7, name = "FantasticFox1", type = "continuous")
reg_colours = c(pal[c(1:4,6:7)],"black","burlywood4")

##attach metadata to tree
p <- ggtree(tree) %<+% meta

##plot tree
pdf("phl_phylogeny1_lab.pdf", height=11, width=8)
p+geom_nodepoint(size=2, shape=15,alpha=.3,aes(label=node,subset = !is.na(as.numeric(label)) & as.numeric(label) > 80))+geom_tippoint(size=4,aes(shape = source, color = region, fill=region))  + theme(legend.text=element_text(size=14), legend.title=element_text(size=14),legend.position = "bottom")+geom_treescale(linesize=0.5, offset=-2)+scale_fill_manual(values=reg_colours)+scale_shape_manual(values=c(24,21))+scale_colour_manual(values=rep("black",8))+theme(legend.position = "bottom",legend.text=element_text(size=14), legend.title=element_text(size=14))+
  guides(fill = guide_legend(override.aes = list(colour = reg_colours))) +theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
#add for tip labels:
geom_tiplab(aes(label=label2), size=1.5, linesize=0, align=F, offset=0.001)

dev.off()

