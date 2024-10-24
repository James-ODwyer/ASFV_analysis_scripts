

library(ggtree)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(phangorn)
library(tidyverse)
library(ggnewscale)
library(ggtreeExtra)
library(ggplot2)
library(adegenet)
library(seqinr)
library(treeio)
library(ggrepel)
library(gtable)
library(scales)
library(grid)
library(TDbook)
library(ape)
library(adegenet)
# Whole genomes 

# denovo ref decided
#tree <- read.iqtree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/trees/PNG_denovo_alignments_23s/no_outliers/2/all_PNG_retree_20_multi-csar_realigned23s_homopolymers_removed_no_ITRs_removed_outliers_2.fasta.treefile")


# Whole genome no ITR 
tree <- read.tree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Genomes_no_ITR_for_tree/all_genomes_multifasta_20240116_mafft_no_ITR.fasta.treefile")


# genes

tree <- read.tree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/results/aligned.treefile")

# ASFV Marker genes

tree <- read.tree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Marker_genes/alignments.treefile")


#tree <- read.tree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/20000boot/Full_genomes/final_fastas_combined_seqs_ASFV_Mafft_align.mafft.treefile")

#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/Asia_pacific_refs/Full_genomes")

#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/Asia_pacific_refs/Genes/MAFFT/FNA")

# Whole genome
#metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/redo_without_viet_error/trees/wgs/metadataASFVPNGs.csv")

metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Genomes_no_ITR_for_tree/metadataASFVPNGs_wholegenomes_trees.csv")

# Genes
#metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/redo_without_viet_error/trees/genes/metadataASFVPNGsgenes.csv")

# Whole genome
#metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/metadataASFV52samples.csv", header = TRUE)

mid_root <- midpoint(tree)


p_mid <- ggtree(mid_root) +
  geom_tiplab(aes(label=label), hjust=.05, size =3.3, align=F) +
  xlim(values=c(0, 0.001)) +
  theme_tree2()
p_mid

# save a tree so that I can see what the genotypes are called

#ggsave("JEV_genome_names.png", height=49, width=49)

# color the branches according to JeV genotype
# first check clade numbers

p_clade <- ggtree(mid_root) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
p_clade


# Look at specific clades as the base
viewClade(p_clade, node=38)

#metadata and trees aren't in order
# create second metadata file, populate as empty
# Then use grep to cycle through and reorder rows




#
metadata$tip_lab <-paste0(metadata$Country,sep="/",metadata$Year,sep="/",metadata$SAN_number)


q <- ggtree(tree)


d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

tree
to_drop_excess_europe <- c("OX376272.1",
             "OX376271.1",
             "OX376268.1",
             "OX376267.1",
             "OX376266.1",
             "OX376265.1",
             "OX376264.1",
             "OX376263.1",
             "OX376262.1",
             "OX376261.1",
             "OX376260.1",
             "OX376259.1",
             "OX376258.1",
             "OX376257.1",
             "OX376255.1",
             "OX376254.1",
             "OX376253.1",
             "OX376252.1",
             "OX376250.1",
             "OM966721.1",
             "OM966720.1",
             "OM966720.1",
             "OM966719.1",
             "OM966714.1",
             "MT847623.2",
             "MT847622.1",
             "OM799941.1",
             "OP823269.1",
             "LC659088.1",
             "MK543947.1",
             "MK333180.1",
             "OM161110.1",
             "OL692744.1",
             "OM481276.1",
             "MZ614662.1",
             "OM966718.1",
             "LC659089.1",
             "OX376251.1",
             "OP628183.1",
             "OM966716.1",
             "MK128995.1")

to_drop_excess_europe_and_africa <- c("OX376272.1",
                           "OX376271.1",
                           "OX376268.1",
                           "OX376267.1",
                           "OX376266.1",
                           "OX376265.1",
                           "OX376264.1",
                           "OX376263.1",
                           "OX376262.1",
                           "OX376261.1",
                           "OX376260.1",
                           "OX376259.1",
                           "OX376258.1",
                           "OX376257.1",
                           "OX376255.1",
                           "OX376254.1",
                           "OX376253.1",
                           "OX376252.1",
                           "OX376250.1",
                           "OM966721.1",
                           "OM966720.1",
                           "OM966720.1",
                           "OM966719.1",
                           "OM966714.1",
                           "MT847623.2",
                           "ON263123.1",
                           "ON456300.2",
                           "MT847622.1",
                           "OM799941.1",
                           "OP823269.1",
                           "LC659088.1",
                           "LR899193.1",
                           "OM966715.1",
                           "MK543947.1",
                           "MK333180.1",
                           #"MT882025.1",
                           #"MT872723.1",
                           "OM161110.1",
                           "OL692744.1",
                           "OM481276.1",
                           "MZ614662.1",
                           "OP628183.1",
                           "OM966718.1",
                           "OX376251.1",
                           "OM966716.1",
                           "MK128995.1",
                           "LC659089.1",
                           "ON409983.1",
                           "MW856068.1",
                           "ON409979.1",
                           "LR813622.1",
                           "1904-02-8314ABCD_S1014_L001_consensus_rd1",
                           "20-01023-02_S4_L001_consensus_rd1",
                           "20-01023-03_S8_L001_consensus_rd1",
                           "LC659086.1",
                           "LC659087.1",
                           "LC659088.1",
                           "LC659089.1"
                           )

p1 <- ggtree(tree) + geom_tiplab(aes(color = label %in% to_drop_excess_europe_and_africa)) +
  scale_color_manual(values=c("black", "red")) + xlim(0, 0.002)

tree_reduced <- drop.tip(tree, to_drop_excess_europe_and_africa)

q <- ggtree(tree_reduced)


d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

ggtree <- ggtree(tree_reduced) %<+% metadata +
  #xlim(values=c(0.0000, 0.0012)) +
  xlim(values=c(0.0000, 0.0004)) +
  geom_tippoint(aes(color=Country, shape=Source), size=3.5, stroke=1) +
  scale_shape_manual(values=c(16, 1)) +
  geom_text_repel(data=d, aes(x = x, y = y,label=label),box.padding = 0.5) + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.2,nudge_y = 0.2), hjust=(-0.1), align=F,size =3.3) +
  #geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 75)) +
  theme_tree2()

#ggtree <- ggtree + geom_label_repel(data=d, aes(label=label))
ggtree

ggsave("PNG_ASF_tree_denovo_reference_aligned_all_genomes_no_ITR.png", height=18, width=14)
ggsave("PNG_ASF_tree_original_reference_aligned.png", height=13, width=14)

#geom_label2(aes(label=label,   subset = !is.na(as.numeric(label)) & as.numeric(label) > 50)) +

#geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
#xlim(values=c(0.016, 0.0217)) +
#geom_nodepoint() +
#scale_fill_manual(values = cols, na.value="grey30") +



# Now to add SNP data to tree. Use Alignment to do this. For full genome
# Next will be gene specific

library(ape)
library(adegenet)


# FULL genome

# analyse SNPs using Investigating mutational differences R script


snp_data$name2 <- NA

for (i in c(1:nrow(snp_data))) {
  
  
  grep(snp_data$name[i],x = metadata$id) -> idx
  snp_data$name2[i] <-  metadata$tip_lab[idx]

  
}
snp_data$name <- snp_data$name2

snp_data <- snp_data[,-3]

# Full SNPs

# It appears that sample names need to be columns and SNPs on rows
library(ggtree)
library(TDbook)

# building off the prior ggtree generated

#p <- ggtree(tree)

## attach the sampling information data set 
## and add symbols colored by location
#p <- p %<+% metadata + geom_tippoint(aes(color=Country))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure

# Whole genome no ITR_all SNPS
#
#
#
#
#snp_data works



# Repeat above for reduced size tree

#Final tree made here including 95 and 80% conf boots and resizing


split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))


# for only when both measures are used
#dots95 = c(rep(FALSE,Ntip(tree_reduced)),as.numeric(tree_reduced$node.label)>95)
#dots80 = c(rep(FALSE,Ntip(tree_reduced)),as.numeric(tree_reduced$node.label)<=95 & as.numeric(tree_reduced$node.label) > 80)

#SingleSNPs_figure_wgs_all_mutation_num <-read.table("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Genomes_no_ITR_for_tree/name_plus_position_all_mutations_for_tree_whole_genome.txt",header = TRUE)
#SingleSNPs_figure_wgs_all_mutation_num <-read.table("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Genomes_no_ITR_for_tree/name_plus_position_common_mutations_for_tree_whole_genome.txt",header = TRUE)
#SingleSNPs_figure_wgs_all_mutation_num <-read.table("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Genomes_no_ITR_for_tree/name_plus_position_all_mutations_min_2_inds_for_tree_whole_genome.txt",header = TRUE)

SingleSNPs
SingleSNPs_figure_wgs <- SingleSNPs[,1:2]
# Make initial ggtree and SNP facet box
ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.1,nudge_y = 0.2), hjust=(-0.04), align=F,size =8.6) + 
  geom_treescale(x=0, y=-2.2, width=0.00096, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=5.6, stroke=1) +  
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=5.3) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=5.3) +
  scale_shape_manual(values=c(16, 1)) + theme_tree2() +
  geom_facet(panel = "SNP position", data = SingleSNPs_figure_wgs, geom = geom_point, mapping=aes(x = pos, color = Country), shape = '|',size=6.5)
ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(3, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 34), # Adjust the size of x-axis labels
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size=26)
  ) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(legend.position = "none")



gt = ggplot_gtable(ggplot_build(ggtree))
#gt$widths[7] = 0.4*gt$widths[7]

gt$widths[9] = 0.25*gt$widths[9]


# redraw ggplot with revised sizing 
grid.draw(gt)




ggsave("ASFV_no_ITR_whole_genome_common.pdf", height=22, width=22,gt)
ggsave("ASFV_no_ITR_whole_genome_common.png", height=22, width=22,gt)



ggsave("ASFV_no_ITR_whole_genome_common_no_leg.pdf", height=22, width=22,gt)
ggsave("ASFV_no_ITR_whole_genome_common_no_leg.png", height=22, width=22,gt)




# Whole genomes. For GARA figs


# GARA remove extras
metadata$tip_lab <- paste0(metadata$Country,sep="-",metadata$Year)

to_drop_excess_all_gara <- c("LR722599.1",
                             "ON456300.2",
                             "OP612151.1",
                             "MK645909.1",
                             "ON380539.1",
                             "MW049116.1",
                             "ON263123.1",
                             "ON380540.1",
                             "OM105586.1",
                             "LR899193.1",
                             "OX376256.1",
                             "OP605386.1",
                             "NC_044948.1",
                             "MH681419.1",
                             "LC659087.1",
                             "LC659086.1",
                             "MW306191.1",
                             "OX376272.1",
                             "OX376271.1",
                             "OX376268.1",
                             "OX376267.1",
                             "OX376266.1",
                             "OX376265.1",
                             "OX376264.1",
                             "OX376263.1",
                             "OX376262.1",
                             "OX376261.1",
                             "OX376260.1",
                             "OX376259.1",
                             "OX376258.1",
                             "OX376257.1",
                             "OX376255.1",
                             "OX376254.1",
                             "OX376253.1",
                             "OX376252.1",
                             "OX376250.1",
                             "OM966721.1",
                             "OM966720.1",
                             "OM966720.1",
                             "OM966719.1",
                             "OM966714.1",
                             "MT847623.2",
                             "ON263123.1",
                             "ON456300.2",
                             "MT847622.1",
                             "OM799941.1",
                             "OP823269.1",
                             "LC659088.1",
                             "LR899193.1",
                             "OM966715.1",
                             "MK543947.1",
                             "MK333180.1",
                             "MT882025.1",
                             "MT872723.1",
                             "OM161110.1",
                             "OL692744.1",
                             "OM481276.1",
                             "MZ614662.1",
                             "OP628183.1",
                             "OM966718.1",
                             "OX376251.1",
                             "OM966716.1",
                             "MK128995.1",
                             "LC659089.1",
                             "ON409983.1",
                             "MW856068.1",
                             "ON409979.1",
                             "LR813622.1")


tree_reduced <- drop.tip(tree, to_drop_excess_all_gara)



split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))

#alphabet_pallete <- c( "#AA0DFE", "#3283FE","#85660D","#782AB6","#565656","#1C8356","#16FF32","#F7E1A0", "#E2E2E2" , "#1CBE4F" , "#C4451C",  "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16","#F8A19F"  ,"#90AD1C"  )

# Want to change the whole genome SNP table to a SNP position out of 600 table. 


#SingleSNPs
#SingleSNPs_figure_wgs <- SingleSNPs[,1:2]
# Make initial ggtree and SNP facet box
ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.2,nudge_y = 0.2), hjust=(-0.1), align=F,size =12.6) + 
  geom_treescale(x=0, y=-2.2, width=0.0010, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=7.3, stroke=1) +  
  #scale_color_manual(values = alphabet_pallete) +
  scale_shape_manual(values = c(17, 16)) +
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=3.3) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=3.3) +
  #scale_shape_manual(values=c(1,16)) + theme_tree2() +
  geom_facet(panel = "SNP", data = SingleSNPs_figure_wgs_all_mutation_num, geom = geom_point, mapping=aes(x = Mutation_number, color = Country), shape = '|',size=3.5)
  ggtree <- facet_labeller(ggtree, c(Tree = "Whole genome phylogeny", SNP = "SNP position"))
  ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(3, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 34), # Adjust the size of x-axis labels
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size=26)
  ) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(strip.text.x = element_text(size = 30))



gt = ggplot_gtable(ggplot_build(ggtree))
gt$widths[7] = 0.44*gt$widths[7]
# redraw ggplot with revised sizing 
grid.draw(gt)




ggsave("ASFV_no_ITR_whole_genome_common_and_fixed_SNPs_GARA.pdf", height=22, width=22,gt)
ggsave("ASFV_no_ITR_whole_genome_common_and_fixed_SNPs_GARA.png", height=22, width=22,gt)





##############


##################################






##################################################



# Repeat above for gene specific analysis

#load("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/redo_without_viet_error/trees/360_fixed_genes/results/Results_SNP_tables_and_Trees_gene_regions_ASFV.Rdata")
tree <- read.tree("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/tree/aligned.treefile")

setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/results")

p1 <- ggtree(tree) + geom_tiplab(aes(color = label %in% to_drop_excess_europe_and_africa)) +
  scale_color_manual(values=c("black", "red")) + xlim(0, 0.002)

tree_reduced <- drop.tip(tree, to_drop_excess_europe_and_africa)


# Repeat above for reduced size tree

#Final tree made here including 95 and 80% conf boots and resizing


#tree_reduced <- drop.tip(tree, to_drop_excess_all_gara)

split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))


# for only when both measures are used
#dots95 = c(rep(FALSE,Ntip(tree_reduced)),as.numeric(tree_reduced$node.label)>95)
#dots80 = c(rep(FALSE,Ntip(tree_reduced)),as.numeric(tree_reduced$node.label)<=95 & as.numeric(tree_reduced$node.label) > 80)

SingleSNPS_allseqs_commmon_regional
SingleSNPs
SingleSNPs_figure_wgs <- SingleSNPS_allseqs_commmon_regional[,2:3]
# Make initial ggtree and SNP facet box
ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.2,nudge_y = 0.2), hjust=(-0.04), align=F,size =8.6) + 
  geom_treescale(x=0, y=-2.2, width=0.00066, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=4.3, stroke=1) +  
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=4.3) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=4.3) +
  scale_shape_manual(values=c(16, 1)) + theme_tree2() +
  geom_facet(panel = "SNP position", data = SingleSNPs_figure_wgs, geom = geom_point, mapping=aes(x = pos, color = Country), shape = '|',size=8.5)
ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(3, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 34), # Adjust the size of x-axis labels
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size=26)
  ) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(legend.position = "none")



gt = ggplot_gtable(ggplot_build(ggtree))
gt$widths[9] = 0.25*gt$widths[7]
# redraw ggplot with revised sizing 
grid.draw(gt)




ggsave("ASFV_genes_common_and_fixed_SNPs.pdf", height=22, width=22,gt)
ggsave("ASFV_genes_common_and_fixed_SNPs.png", height=22, width=22,gt)

ggsave("ASFV_genes_common_and_fixed_SNPs_no_legend.pdf", height=22, width=22,gt)
ggsave("ASFV_genes_common_and_fixed_SNPs_no_legend.png", height=22, width=22,gt)


# Whole genomes. For GARA figs


# GARA remove extras
metadata$tip_lab <- paste0(metadata$Country,sep="-",metadata$Year)

to_drop_excess_all_gara <- c("LR722599.1",
                             "ON456300.2",
                             "OP612151.1",
                             "MK645909.1",
                             "ON380539.1",
                             "MW049116.1",
                             "ON263123.1",
                             "ON380540.1",
                             "OM105586.1",
                             "LR899193.1",
                             "OX376256.1",
                             "OP605386.1",
                             "NC_044948.1",
                             "MH681419.1",
                             "LC659087.1",
                             "LC659086.1",
                             "MW306191.1",
                             "OX376272.1",
                             "OX376271.1",
                             "OX376268.1",
                             "OX376267.1",
                             "OX376266.1",
                             "OX376265.1",
                             "OX376264.1",
                             "OX376263.1",
                             "OX376262.1",
                             "OX376261.1",
                             "OX376260.1",
                             "OX376259.1",
                             "OX376258.1",
                             "OX376257.1",
                             "OX376255.1",
                             "OX376254.1",
                             "OX376253.1",
                             "OX376252.1",
                             "OX376250.1",
                             "OM966721.1",
                             "OM966720.1",
                             "OM966720.1",
                             "OM966719.1",
                             "OM966714.1",
                             "MT847623.2",
                             "ON263123.1",
                             "ON456300.2",
                             "MT847622.1",
                             "OM799941.1",
                             "OP823269.1",
                             "LC659088.1",
                             "LR899193.1",
                             "OM966715.1",
                             "MK543947.1",
                             "MK333180.1",
                             "MT882025.1",
                             "MT872723.1",
                             "OM161110.1",
                             "OL692744.1",
                             "OM481276.1",
                             "MZ614662.1",
                             "OP628183.1",
                             "OM966718.1",
                             "OX376251.1",
                             "OM966716.1",
                             "MK128995.1",
                             "LC659089.1",
                             "ON409983.1",
                             "MW856068.1",
                             "ON409979.1",
                             "LR813622.1")


tree_reduced <- drop.tip(tree, to_drop_excess_all_gara)



split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))


###

#SingleSNPs
#SingleSNPs_figure_wgs <- SingleSNPs[,1:2]
# Make initial ggtree and SNP facet box
ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.2,nudge_y = 0.2), hjust=(-0.1), align=F,size =12.3) + 
  geom_treescale(x=0, y=-2.2, width=0.0006, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=7.3, stroke=1) +  
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=3.3) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=3.3) +
  #scale_shape_manual(values=c(16, 1)) + theme_tree2() +
  scale_shape_manual(values = c(17, 16)) +theme_tree2() +
  geom_facet(panel = "SNP position", data = SingleSNPS2col, geom = geom_point, mapping=aes(x = pos, color = Country), shape = '|',size=11.5)
ggtree <- facet_labeller(ggtree, c(Tree = "Core gene phylogeny", SNP = "SNP position"))
ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(2, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 26.4), # Adjust the size of x-axis labels
    legend.text = element_text(size = 26),
    plot.title = element_text(size = 21),
    legend.title = element_text(size=28)
  ) +
  theme(strip.text.x = element_text(size = 30))



gt = ggplot_gtable(ggplot_build(ggtree))
gt$widths[7] = 0.36*gt$widths[7]
# redraw ggplot with revised sizing 
grid.draw(gt)




ggsave("ASFV_no_ITR_genes_common_and_fixed_SNPs_GARA.pdf", height=22, width=22,gt)
ggsave("ASFV_no_ITR_genes_common_and_fixed_SNPs_GARA.png", height=22, width=22,gt)








##############################################




###############################################





##############################################







































































####






###########

############

# Whole genome restricted SNPS to common within 1-2 countries. 

############
metadata$tip_lab <-paste0(metadata$Country,sep="/",metadata$Year,sep="/",metadata$SAN_number)

for (i in c(1:nrow(SingleSNPs))) {
  
  
  grep(SingleSNPs$name[i],x = metadata$id) -> idx
  SingleSNPs$name2[i] <-  metadata$tip_lab[idx]
  
  
}

SingleSNPs$name <- SingleSNPs$name2

SingleSNPs <- SingleSNPs[,1:2]



ggtree <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label=tip_lab), hjust=(-0.1), align=F,size =3.3) +
  geom_treescale(x=0, y=-2.2, width=0.0012, color='white') +
  #xlim(values=c(0.0000, 0.0012)) +
  geom_tippoint(aes(color=Country, shape=Source), size=3.5, stroke=1) +
  scale_shape_manual(values=c(16, 1)) +
  #geom_text_repel(data=d, aes(label=label),box.padding = 0.001) + 
  #geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 75)) +
  theme_tree2()

ggtree <- ggtree + geom_facet(panel = "SNP", data = SingleSNPs, geom = geom_point, 
                              mapping=aes(x = pos, color = Country), shape = '|')
#+
#  theme_tree2(legend.position=c(.05, .85))

ggsave("ASFV_no_ITR_whole_genome_common_SNPs_only.pdf", height=26, width=24)





# Gene specific analysis




metadata$tip_lab <-paste0(metadata$Country,sep="/",metadata$Year,sep="/",metadata$SAN_number)


q <- ggtree(tree)

d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

ggtree <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label=tip_lab), hjust=(-0.1), align=F,size =3.3) +
  xlim(values=c(0.0000, 0.0005)) +
  geom_tippoint(aes(color=Country, shape=Source), size=3.5, stroke=1) +
  scale_shape_manual(values=c(16, 1)) +
  geom_text_repel(data=d, aes(label=label),box.padding = 0.3) + 
  #geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 75)) +
  theme_tree2()

ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label=tip_lab), hjust=(-0.1), align=F,size =3.5) +
  geom_text(aes(label=node), hjust=-.3)

ggtree

ggsave("PNG_ASF_tree_denovo_reference_aligned_all_genomes_gene_regions.png", height=28, width=14)


# Now to add the SNP panel to the genes 


metadata$tip_lab <-paste0(metadata$Country,sep="/",metadata$Year,sep="/",metadata$SAN_number)


#SingleSNPs$name <- SingleSNPs$name2



SingleSNPS2col <- as.data.frame(matrix(nrow=nrow(SingleSNPs),ncol=2))

SingleSNPs$name -> SingleSNPS2col[,1]
SingleSNPs$pos -> SingleSNPS2col[,2]

colnames(SingleSNPS2col) <- c("id","pos")



ggtree <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label=tip_lab), hjust=(-0.1), align=F,size =3.3) +
  geom_treescale(x=0, y=-2.2, width=0.0012, color='white') +
  #xlim(values=c(0.0000, 0.0012)) +
  geom_tippoint(aes(color=Country, shape=Source), size=3.5, stroke=1) +
  scale_shape_manual(values=c(16, 1)) +
  #geom_text_repel(data=d, aes(label=label),box.padding = 0.001) + 
  #geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 95)) +
  theme_tree2()

ggtree <- ggtree + geom_facet(panel = "SNP", data = SingleSNPS2col, geom = geom_point, 
                              mapping=aes(x = pos, color = Country), shape = '|', size=5)

ggtree <- ggtree + scale_x_continuous(labels = label_comma())

ggtree <- ggtree +
  theme(
    legend.key.size = unit(2, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 15), # Adjust the size of x-axis labels
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16),
  ) +
  theme(strip.text.x = element_text(size = 30))


ggsave("ASFV_core_genes_common_SNPs_only.pdf", height=22, width=22)
ggsave("ASFV_core_genes_common_SNPs_only.png", height=22, width=22)


# SNP 2 column works as of 20/06/23 so use this one


# I have added both forms of bootstrapping (Left is branch specific, right is UFboots)

split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))

# Make initial ggtree and SNP facet box
ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.2,nudge_y = 0.2), hjust=(-0.1), align=F,size =6.5) + 
  geom_treescale(x=0, y=-2.2, width=0.00061, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=4.3, stroke=1) +  
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=4.5) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=4.5) +
  scale_shape_manual(values=c(16, 1)) + theme_tree2() +
  geom_facet(panel = "SNP", data = SingleSNPS2col, geom = geom_point, mapping=aes(x = pos, color = Country), shape = '|',size=5.5)

ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(2, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 22), # Adjust the size of x-axis labels
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 20),
    legend.title = element_text(size=22)
  ) +
  theme(strip.text.x = element_text(size = 30))



gt = ggplot_gtable(ggplot_build(ggtree))
gt$widths[7] = 0.4*gt$widths[7]
# redraw ggplot with revised sizing 
grid.draw(gt)




ggsave("ASFV_no_ITR_whole_genome_common_and_fixed_SNPs_test.pdf", height=22, width=22,gt)
ggsave("ASFV_no_ITR_whole_genome_common_and_fixed_SNPs.png", height=22, width=22,gt)








SingleSNPs
SingleSNPs_figure_genes <- SingleSNPs[,2:3]
ggtree <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label=tip_lab), hjust=(-0.1), align=F,size =3.8) +
  geom_treescale(x=0, y=-2.2, width=0.00115, color='white') +
  #xlim(values=c(0.0000, 0.0012)) +
  geom_tippoint(aes(color=Country, shape=Source), size=3.5, stroke=1) +
  scale_shape_manual(values=c(16, 1)) +
  #geom_text_repel(data=d, aes(label=label),box.padding = 0.001) + 
  #geom_text2(aes(subset = !isTip, label=label),vjust=0.0001) +
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 75)) +
  theme_tree2()

ggtree <- ggtree + geom_facet(panel = "SNP", data = SingleSNPs_figure_genes, geom = geom_point, 
                              mapping=aes(x = pos, color = Country), shape = '|',size=5)



ggsave("ASFV_genes_common_SNPs_only.pdf", height=26, width=24)










SingleSNPS2col <- as.data.frame(matrix(nrow=nrow(SingleSNPs),ncol=2))

SingleSNPs$name -> SingleSNPS2col[,1]
SingleSNPs$pos -> SingleSNPS2col[,2]

colnames(SingleSNPS2col) <- c("id","pos")








# Only fix for the full genome by changing the id names slightly to fit the whole genome format
# Only required for full genome 
#metadata$id <- metadata$id.1
Study <- data.frame("Study" = metadata[,c("Study")])
rownames(Study) <- metadata$id
# Need to specify a new column for metadata for full genome (MAFFT fomatted the names differently to the gene alignments)

metadata$tip_lab <-paste0(metadata$SAN_number,sep="/",metadata$Country,sep="/",metadata$Year)

# Sub in the following into the SNP section of after filtering out additional non informative snps 
# SingleSNPsalllocals
#SingleSNPsalllocals <- SingleSNPsalllocals[,1:2]

p <- ggtree(tree)
p <- p + geom_treescale(x=0, y=-2.2, width=0.0010, color='black')
p <- p %<+% metadata + geom_tiplab(aes(label=tip_lab), hjust=0.0003, align=F,size =6.0)
p <- p + geom_tippoint(aes(color=Country),size=3.8)
## attach the sampling information data set 
## and add symbols colored by location
#p <- p %<+% metadata + geom_tippoint(aes(color=Country))
#p <- p %<+% metadata + geom_tiplab(aes(label=tip_lab), hjust=0, align=F,size =3.5)
#offset = 0.00006
h1 <- p + new_scale_fill()
h1 <- h1 + scale_x_continuous(labels = scales::comma_format())
h1 <-  gheatmap(p, Study,                                 # we add a heatmap layer of the study dataframe to our tree plot
                offset = 0.00028, width = 0.010 ,                             # offset shifts the heatmap to the right,                         # width defines the width of the heatmap column,                              # color defines the boarder of the heatmap columns
                colnames = FALSE) +
  scale_fill_manual(name = "Study",                       # define the coloring scheme and legend for gender
                    values = c("#00d1b1", "purple"),
                    breaks = c("This study", "Reference"),
                    labels = c("This study", "Genbank reference")) +
  theme_tree2(legend.position=c(.05, .85)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        title =element_text(size=22, face='bold'),
        legend.box = "horizontal", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 2,byrow = TRUE))

# Picking up from complete snp data from SNP analysis, just need to delete genes 
#snp_data <- snp_data[,-1]

h2 <- h1 + geom_facet(panel = "SNP", data = SingleSNPs_figure_wgs, geom = geom_point,
                      mapping=aes(x = pos, color = Country), shape = '|', size = 5) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

h2

ggsave("ASFV_no_ITR_whole_genome_common_and_fixed_SNPs.png", height=26, width=28)



###############################


# Repeat with the revised down only single/double country hits


###########################
SingleSNPsoneortwo_short <- SingleSNPsoneortwo[,1:2]

# To do all SNPs for all locations
SingleSNPssubset <- snp_data_meta[,2:3]

#to do SNPs for all locations but restricted to common/fixed/prevelant

# Only this study SNPs that are country specific
SingleSNPssubset <- SingleSNPs[,2:3]

# All SNPs in specific countries
SingleSNPssubset <- SingleSNPsalllocals[,2:3]


#metadata$id <- metadata$id.1
Study <- data.frame("Study" = metadata[,c("Study")])
rownames(Study) <- metadata$id


p <- ggtree(tree) + theme_tree2()
#p <- p + geom_treescale(x=0, y=-2.2, color='black') + theme_tree2()
p <- p %<+% metadata + geom_tiplab(aes(label=tip_lab), hjust=0.0000002, align=F,size =5.8)
p <- p + geom_tippoint(aes(color=Country),size=4.4) + theme_tree2()
p <- p + geom_treescale(x=0, y=-2.2, color='black') 
## attach the sampling information data set 
## and add symbols colored by location
#p <- p %<+% metadata + geom_tippoint(aes(color=Country))
#p <- p %<+% metadata + geom_tiplab(aes(label=tip_lab), hjust=0, align=F,size =3.5)
#offset = 0.00006
h1 <- p + new_scale_fill()
h1 <-  gheatmap(p, Study,                                 # we add a heatmap layer of the gender dataframe to our tree plot
                offset = 0.0003, width = 0.020 ,                             # offset shifts the heatmap to the right,                         # width defines the width of the heatmap column,                              # color defines the boarder of the heatmap columns
                colnames = FALSE) +
  scale_fill_manual(name = "Study",                       # define the coloring scheme and legend for gender
                    values = c("#00d1b1", "purple"),
                    breaks = c("This study", "Reference"),
                    labels = c("This study", "Genbank reference")) +
  theme_tree2(legend.position=c(.05, .85)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        title =element_text(size=22, face='bold'),
        legend.box = "horizontal", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 2,byrow = TRUE))

# Picking up from complete snp data from SNP analysis, just need to delete genes 
#snp_data <- snp_data[,-1]

h2 <- h1 + geom_facet(panel = "SNP", data = SingleSNPS2col, geom = geom_point, size=10,
                      mapping=aes(x = pos, color = Country), shape = '|') +
  theme(strip.text.x = element_text(size = 38)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(plot.margin = margin(1.5,1.5,1.5,1.5, "cm")) 

h2

ggsave("ASFV_whole_genome_no_ITR_phylogenetic_tree_with_local_SNPs.pdf", height=26, width=28)
ggsave("ASFV_whole_genome_no_ITR_phylogenetic_tree_with_local_SNPs.png", height=26, width=28)

ggsave("ASFV_gene_gene_comparison_phylogenetic_tree_with_local_SNPs.pdf", height=26, width=28)
ggsave("ASFV_gene_gene_comparison_phylogenetic_tree_with_local_SNPs.png", height=26, width=28)

# also works for 2 columns SNPs as of 20/06/23







# Generate marker gene trees


to_drop_excess_europe_and_africa <- c("OX376272.1",
                                      "OX376271.1",
                                      "OX376268.1",
                                      "OX376267.1",
                                      "OX376266.1",
                                      "OX376265.1",
                                      "OX376264.1",
                                      "OX376263.1",
                                      "OX376262.1",
                                      "OX376261.1",
                                      "OX376260.1",
                                      "OX376259.1",
                                      "OX376258.1",
                                      "OX376257.1",
                                      "OX376255.1",
                                      "OX376254.1",
                                      "OX376253.1",
                                      "OX376252.1",
                                      "OX376250.1",
                                      "OM966721.1",
                                      "OM966720.1",
                                      "OM966720.1",
                                      "OM966719.1",
                                      "OM966714.1",
                                      "MT847623.2",
                                      "ON263123.1",
                                      "ON456300.2",
                                      "MT847622.1",
                                      "OM799941.1",
                                      "OP823269.1",
                                      "LC659088.1",
                                      "LR899193.1",
                                      "OM966715.1",
                                      "MK543947.1",
                                      "MK333180.1",
                                      #"MT882025.1",
                                      #"MT872723.1",
                                      "OM161110.1",
                                      "OL692744.1",
                                      "OM481276.1",
                                      "MZ614662.1",
                                      "OP628183.1",
                                      "OM966718.1",
                                      "OX376251.1",
                                      "OM966716.1",
                                      "MK128995.1",
                                      "LC659089.1",
                                      "ON409983.1",
                                      "MW856068.1",
                                      "ON409979.1",
                                      "LR813622.1",
                                      "1904-02-8314ABCD_S1014_L001_consensus_rd1",
                                      "20-01023-02_S4_L001_consensus_rd1",
                                      "20-01023-03_S8_L001_consensus_rd1",
                                      "LC659086.1",
                                      "LC659087.1",
                                      "LC659088.1",
                                      "LC659089.1",
                                      "MN194591.1",
                                      "NC_044948.1"
)

tree_reduced <- drop.tip(tree, to_drop_excess_europe_and_africa)


split_values <- strsplit(tree_reduced$node.label, "/")

# Extracting X and Y into separate vectors
alrt_values <- sapply(split_values, function(x) x[1])
UFboots <- sapply(split_values, function(x) x[2])

boots= as.data.frame(matrix(nrow= length(tree_reduced$node.label),ncol=2))
boots$V1 <- as.numeric(alrt_values)
boots$V2 <- as.numeric(UFboots)


dots95 = c(rep(FALSE,Ntip(tree_reduced)),(boots$V1>95 & boots$V2>95))
dots80 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1<95 & boots$V2<95 & boots$V1>80 & boots$V2>80) | (boots$V1<95 & boots$V2>95 & boots$V1>80) | (boots$V1>95 & boots$V2<95 & boots$V2>80)))

dots20 = c(rep(FALSE,Ntip(tree_reduced)),((boots$V1>10) | (boots$V2>10) ))




ggtree <- ggtree(tree_reduced) %<+% metadata + 
  geom_tiplab(aes(label=tip_lab,geom = "label_repel",direction = "both",method = "lastb",nudge_x = 0.1,nudge_y = 0.2), hjust=(-0.04), align=F,size =8.6) + 
  geom_treescale(x=0, y=-2.2, width=0.00096, color='white') + geom_tippoint(aes(color=Country, shape=Source), size=5.6, stroke=1) +  
  geom_nodepoint(aes(subset=dots95),shape=22,fill='black',size=5.3) + 
  geom_nodepoint(aes(subset=dots80),shape=22,fill='white',size=5.3) +
  scale_shape_manual(values=c(16, 1)) + theme_tree2() +
  xlim(values=c(0.0000, 0.0003)) +
  theme(
    legend.key.size = unit(3, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 34), # Adjust the size of x-axis labels
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size=26)
  ) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(legend.position = "none")






ggtree <- ggtree + scale_x_continuous(labels = label_comma())
ggtree <- ggtree +
  theme(
    legend.key.size = unit(3, 'lines'),  # Adjust the legend key size
    axis.text.x = element_text(size = 34), # Adjust the size of x-axis labels
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 18),
    legend.title = element_text(size=26)
  ) +
  theme(strip.text.x = element_text(size = 38)) +
  theme(legend.position = "none") +
  xlim(values=c(0.0000, 0.0001)) +
  scale_x_continuous(labels = label_comma())




ggsave("ASFV_marker_genes.pdf", height=26, width=28)
ggsave("ASFV_marker_genes.png", height=26, width=28)








