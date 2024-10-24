# Analysis of gene mutations and positions of mutations in ASFV 
library(ape)
library(adegenet)
library(seqinr)
library(dplyr)
library(ComplexHeatmap)
library(Biostrings)


#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/trees/All_genomes_rd2_consensus_full_genome_no_ITR/genes")
setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/aligned")
fastas <- list.files(path =getwd(),pattern = "fasta")
snplist <- list()
for (i in c( 1:length(fastas))) {
  
  
  FNAdata <- seqinr::read.fasta(file= fastas[i],seqtype = "DNA", as.string = FALSE)
  
  FNAdata[] <- lapply(FNAdata, c)
  
  FNAdata2 <- bind_rows(FNAdata,.id = "Sampleid")
  
  
  snps <- as.data.frame(FNAdata2)
  names <- colnames(snps)
  gapChar <- "n"
  snps %>% relocate("NC_044959.2") -> snps
  snp <- t(snps)
  # Need rownames to be samples
  lsnp <- apply(snp, 1, function(x) {
    x != snp[1,] & x != gapChar & snp[1,] != gapChar
  })
  lsnp <- as.data.frame(lsnp)
  
  
  if (!is.null(nrow(lsnp))) {
    if (isTRUE(grep("lsnp",x = colnames(lsnp)) >=1)) {
      lsnp <- t(lsnp)
    }
    positions <- rownames(snps)
    lsnp <- as.data.frame(lsnp)
    lsnp$pos <- as.numeric(positions)
    lsnp <- tidyr::gather(lsnp, name, value, -pos)
    snp_data <- lsnp[lsnp$value, c("name", "pos")]
    
    snp_data <- as.data.frame(snp_data)

      snplist[[i]] <- snp_data
    }
    
    if (nrow(snp_data) ==0) {
      snplist[[i]] <- NULL
    }
  }









heatmapslist <- list()
generanges <- list()

for (i in c( 1:length(fastas))) {

FNA <-read.alignment(fastas[i],format = "fasta")
aa <- dist.alignment(FNA, matrix = "similarity")
aa1 <- as.matrix(aa)


generanges[[i]] <- range(aa1)

heatmapslist[[i]] <- Heatmap(aa1)

}

heatmapslist[[11]]
# Combined genes looks a lot cooler. Can clearly see the demarkation between regions (Vietnam very dissimilar, PNG very similar to each other. German samples similar etc)




#1.  Define all SNPs 
#2. Add metadata to the SNP file
#3. Sort/filter to see SNPs found in only 1 country 


# 1.

# FNA files but with Georgia root defining SNPs  
#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/Asia_pacific_refs/Genes/MAFFT/FNA/alignments")



names(snplist) <- fastas


snplist

snplist2 <-bind_rows(snplist, .id = "gene")

snplist2$genepos <- snplist2$pos


snplist2$gene <- sub(pattern = ".mafft",replacement = "",x = snplist2$gene)

snplist2$gene <- sub(pattern = "_fixed",replacement = "",x = snplist2$gene)
snplist2$gene <- sub(pattern = ".fasta",replacement = "",x = snplist2$gene)
snplist2$gene <- sub(pattern = "_aligned",replacement = "",x = snplist2$gene)

snplist2$gene <- gsub("^[0-9]*_",replacement = "",x = snplist2$gene)

geneandpos <-read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/attempt_multicsar/trees/All_genomes_rd2_consensus_full_genome_no_ITR/genes/geneposfile.csv",header=FALSE)
snplist2$gene <- gsub(x = snplist2$gene,pattern = "_CDS",replacement = "")

table(snplist2$name)
snplist2$pos <- as.numeric(snplist2$pos)
snplist3 <- snplist2

#for ( j in c(1:190)) {
  
 # lowerlim=as.numeric(geneandpos[j,1])
  #upperlim=as.numeric(geneandpos[j,2])
  #gene=geneandpos[j,3]
  
#for (i in c(1:nrow(snplist2))) {
  
#  if (snplist2$pos[i] <= upperlim && snplist2$pos[i] >= lowerlim) {
    
#    snplist2$gene[i] <- gene
    
#  }

  
 #   }
#}

for ( i in c(1:nrow(snplist2))) {
  
  grep(snplist2$gene[i],x = geneandpos$V3) -> idx
  idx <- idx[1]
  snplist3$pos[i] <-  (snplist2$pos[i] + geneandpos$V1[idx])
  
}


snp_data <- snplist3 
#snp_data <- snp_data[,-1]

#

#
##
#

#2.

metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/metadataASFVPNGs.csv", header = TRUE)


snp_databothpos <- snp_data
snp_databothpos$posgene <- snplist2$pos
snp_data_meta <- snp_databothpos
snp_data_meta$Country <- NA
snp_data_meta$Study <- NA



for (i in c(1:nrow(snp_data_meta))) {
  
  
  grep(snp_data_meta$name[i],x = metadata$id) -> idx 
  
  snp_data_meta$Country[i] <- metadata$Country[idx]
  snp_data_meta$Study[i] <- metadata$Study[idx]
  
}

uniqpositions<- unique(snp_data_meta$pos)
snp_data_meta$SNPtype <- NA


for (i in c(1:length(uniqpositions))) {
  
  
  snp_subset <- subset(snp_data_meta,snp_data_meta$pos==uniqpositions[i]) 
  
  countrynum <- length(unique(snp_subset$Country))
  samplenum <- length(unique(snp_subset$name))
  
  if(countrynum<2){
    totalsamplesfromcountry <- sum(grepl(snp_subset$Country[1],metadata$Country))
    
    
    
    if(samplenum>=(0.5*totalsamplesfromcountry)) {
      
      information <- "Common single country"
      
    }
    
    if(samplenum>=2 && samplenum<(0.5*totalsamplesfromcountry)) {
      
      information <- "Present single country"
      
    }
    
    if(totalsamplesfromcountry==samplenum) {
      
      information <- "fixed single country"
      
    }
    
  }
  
  
  if(countrynum>=2 && countrynum<=4){
    
    information <- "Multi country"
  }
  
  if(countrynum==2) {
    countries <- unique(snp_subset$Country)
    totalsamplesfromcountry1 <- sum(grepl(countries[1],metadata$Country))
    totalsamplesfromcountry2 <- sum(grepl(countries[2],metadata$Country)) 
    totalsamplesfromcountry <- (totalsamplesfromcountry1 + totalsamplesfromcountry2)
    if(samplenum>=(0.5*totalsamplesfromcountry)) {
      
      information <- "common two countries"
      
    }
    
    if(samplenum>=2 && samplenum<(0.5*totalsamplesfromcountry)) {
      
      information <- "present two countries"
      
    }
    
    if(totalsamplesfromcountry==samplenum) {
      
      information <- "fixed two countries"
      
    }
    
  }
  
  
  if(countrynum>=5){
    
    information <- "Widespread"
  }
  if(samplenum==1){
    
    information <- "non_informative"
  }
  
  grep(pattern = paste0("^",uniqpositions[i],"$"),x = snp_data_meta$pos) -> idx
  
  snp_data_meta$SNPtype[idx] <- information
  
}
table(snp_data_meta$SNPtype)

# So 35 substitutions (not unique positions but positions multiplied by samples) were common (found in over 50% of samples from that country)
# 17,059 substitutions (98.1%) were found across more than 5 countries and are widespread 
# 5 substitutions were found across multiple but fewer than 5 countries.
# 84 substitutions were found in a single country 
# 36 were found fixed
# 166 were found to be sample specific


aggregatedstats <-aggregate(as.factor(snp_data_meta$Country),by=list(snp_data_meta$SNPtype,snp_data_meta$Country),length)

SingleSNPs <- snp_data_meta
SingleSNPs <- subset(snp_data_meta,snp_data_meta$SNPtype %in% c("Common single country","fixed single country","common two countries", "fixed two countries") )
#SingleSNPs <- subset(snp_data_meta,snp_data_meta$SNPtype %in% c("Common single country","fixed single country") )
SingleSNPs <- subset(SingleSNPs,SingleSNPs$Study== "This study" )

# removing ambig code Indian ref seq 




PNGsingles <- subset(SingleSNPs,SingleSNPs$Country=="PNG")

unique(PNGsingles$pos)
unique(PNGsingles$gene)
# 19 substitution sites across 6 genes


Timor_leste <- subset(SingleSNPs,SingleSNPs$Country=="Timor-Leste")

unique(Timor_leste$pos)
unique(Timor_leste$gene)
# Only two genomes for East Timor so can a statement about fixation be made 
# 1 substitution sites in one genes was found in both East Timor samples. 


Vietnamsingles <- subset(SingleSNPs,SingleSNPs$Country=="Viet nam")

unique(Vietnamsingles$pos)
unique(Vietnamsingles$gene)



Germanysingles <- subset(SingleSNPs,SingleSNPs$Country=="Germany")

unique(Germanysingles$pos)
unique(Germanysingles$gene)


Indiasingles <- subset(SingleSNPs,SingleSNPs$Country=="India")

unique(Indiasingles$pos)
unique(Indiasingles$gene)


Russiasingles <- subset(SingleSNPs,SingleSNPs$Country=="Russia")

unique(Russiasingles$pos)
unique(Russiasingles$gene)

Southkoreasingles <- subset(SingleSNPs,SingleSNPs$Country=="South Korea")

unique(Southkoreasingles$pos)
unique(Southkoreasingles$gene)


# 16 sites within the one gene MGF_360_21R were found in two Vietnam samples sequenced.

# Have done the three main countries where samples for this study were taken from 
# Can do all but discuss with David and Matt
# 

# Next, need to get the gene specific SNP position for each important substitution and translate to AA to see whether they are synonymous vs nonsynonymous
# Added gene position into above code. Now will re extract relevant base from fastas
# 

SingleSNPs$refbase <- NA
SingleSNPs$refbasesequence <- NA
SingleSNPs$changedbase <- NA
SingleSNPs$sequence <- NA
setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/aligned")


fastas <- list.files(path =getwd(),pattern = "fasta")

for (i in c( 1:length(fastas))) {
  
  FNAdata <- seqinr::read.fasta(file= fastas[i],seqtype = "DNA", as.string = FALSE)
  
  FNAdata[] <- lapply(FNAdata, c)
  
  FNAdata2 <- bind_rows(FNAdata,.id = "Sampleid")
  
  snps <- as.data.frame(FNAdata2)
  
  for (l in c(1:nrow(snps))) {
    
    if (sum(grep("n",snps[l,])) >=1) {
      grep("n",snps[l,]) -> idx
      baseposns <- snps[l,]
      baserefns <- apply(baseposns[], 1, function(x) {names(which.max(table(factor(x,unique(x)))))})
      
      snps[l,idx] <- baserefns
      
    }
    
  }
  # There is a single ambiguous code from the reference sequence from India. It crashes the AA. Changing to standard for that position to
  for (l in c(1:nrow(snps))) {
    
    
    if (sum(grep("m|r|y|k|b|d|v|h|w|s",snps[l,])) >=1) {
      grep("m|r|y|k|b|d|v|h|w|s",snps[l,]) -> idx
      baseposns <- snps[l,]
      baserefns <- apply(baseposns[], 1, function(x) {names(which.max(table(factor(x,unique(x)))))})
      
      snps[l,idx] <- baserefns
      
    }
    
    
  }
  
  
  gene <- sub(pattern = ".mafft",replacement = "",x = fastas[i])
  
  gene <- sub(pattern = "_fixed",replacement = "",x = gene)
  
  
  gene <- sub(pattern = ".fasta",replacement = "",x = gene)
  gene <- sub(pattern = "_aligned",replacement = "",x = gene)
  
  gene <- gsub("^[0-9]*_",replacement = "",x = gene)
  
  gene <- gsub("_CDS",replacement = "",x = gene)
  
  #gene <- gsub("_cds",replacement = "",x = gene)
  
  grep(gene,x = SingleSNPs$gene) -> idx
  
  if (length(idx)>=1) {
    
    for (j in c(1:length(idx))) {
      
      target<- SingleSNPs[idx[j],]   
      sampname <- target$name 
      basepos <- target$posgene  
      grep(sampname,colnames(snps)) -> samploc
      
      sequence <- paste(snps[,samploc],collapse = "")
      allseqs <- as.data.frame(matrix(nrow=ncol(snps),ncol=1))
      for ( k in c(1:ncol(snps))) {
        
        allseqs[k,1] <- paste(snps[,k],collapse = "")
      }
      allbasesposition <- snps[basepos,]
      baseref <- apply(allbasesposition[], 1, function(x) {names(which.max(table(factor(x,unique(x)))))})
      baserefsequence <- apply(allseqs[], 2, function(x) {names(which.max(table(factor(x,unique(x)))))})
      altbase <- substring(sequence,basepos,basepos)
      names(which.max(table(factor(allseqs$V1,unique(allseqs$V1)))))
      
      #grep
      # Need to change n's to standard now!
      
      SingleSNPs[idx[j],9] <- baseref
      SingleSNPs[idx[j],10] <- baserefsequence
      SingleSNPs[idx[j],11] <- altbase
      SingleSNPs[idx[j],12] <- sequence
      
      
    }
    
  }  
  
  
  
  
  
}

rownames(SingleSNPs) <- c(1:nrow(SingleSNPs))

SingleSNPs <-subset(SingleSNPs,SingleSNPs$refbase != SingleSNPs$changedbase)

SingleSNPs$changedbase






#####################





####################






#######################






# Ran into an issue, The above works for substitutions but not indels. These can throw the whole AA chain out of whack (or more likely 
#there are sets of 3 each when they occur)
#
# Maybe a better method is to start from scratch with the fasta alignments, convert all n's to their average base and run each through translation to compare
# the whole sequence. 



# I inserted this code into the above SNPs analysis for SingleSNPs 

##############


# First. Identify how much of an issue - and ambig codes are 
SingleSNPs$refbase
SingleSNPs$changedbase

SingleSNPsnoambig <- subset(SingleSNPs,SingleSNPs$changedbase == "a"| SingleSNPs$changedbase == "t" | SingleSNPs$changedbase ==  "c" | SingleSNPs$changedbase ==  "g" | SingleSNPs$changedbase ==  "-")


SingleSNPsnoambig$changedbase

# Now picking up from the end of generating single SNPs to now generate complete AA codes 


SingleSNPsnoambig$refseqshort <- NA 
SingleSNPsnoambig$altseqshort <- NA 
SingleSNPsnoambig$refAA <- NA 
SingleSNPsnoambig$altAA <- NA 
SingleSNPsnoambig$identical <- NA 
SingleSNPsnoambig$refgenelength <- NA
SingleSNPsnoambig$altgenelength <- NA
SingleSNPsnoambig$refAAlength <- NA
SingleSNPsnoambig$altAAlength <- NA

for ( i in c(1:nrow(SingleSNPsnoambig))) {
  
  # to extract the relevant codon will need to identify what bp position the mutation is and extract around it
  # divide by 3 and whole number = end base of codon (add prior two bases)
  # divide by 3 and x.67 middle base of codon (add one prior and one after base)
  # divide by three and x.33 start of codon (add two after)
  # 
  SingleSNPsnoambig$refseqshort[i] <- gsub(pattern = "-",replacement = "",x = SingleSNPsnoambig$refbasesequence[i])
  SingleSNPsnoambig$altseqshort[i] <- gsub(pattern = "-",replacement = "",x = SingleSNPsnoambig$sequence[i])
  
  AAstringref <- Biostrings::translate(DNAStringSet(SingleSNPsnoambig$refseqshort[i]))
  AAstringalt <- Biostrings::translate(DNAStringSet(SingleSNPsnoambig$altseqshort[i]))
  
  SingleSNPsnoambig$refAA[i] <- as.character(unlist(AAstringref))
  SingleSNPsnoambig$altAA[i] <- as.character(unlist(AAstringalt)) 
  
  SingleSNPsnoambig$identical[i] <- identical(x = SingleSNPsnoambig$refAA[i],SingleSNPsnoambig$altAA[i])
  SingleSNPsnoambig$refgenelength[i] <- nchar(SingleSNPsnoambig$refseqshort[i])
  SingleSNPsnoambig$altgenelength[i] <- nchar(SingleSNPsnoambig$altseqshort[i])
  SingleSNPsnoambig$refAAlength[i] <- nchar(SingleSNPsnoambig$refAA[i])
  SingleSNPsnoambig$altAAlength[i] <- nchar(SingleSNPsnoambig$altAA[i])
  
}

proteininspect <- SingleSNPsnoambig[,c(1:2,15:16)]

# There are a small number of nonesense mutations but appear to be only in refs seq


SingleSNPssynonymous <- subset(SingleSNPsnoambig,SingleSNPsnoambig$identical==TRUE)
SingleSNPsnonsynonymous <- subset(SingleSNPsnoambig,SingleSNPsnoambig$identical==FALSE)

# Some weired no change inds that reported identical bases.
# it must be the ambiguious seqs. I changed the ambigs to the average which would come up as a mutation still! They
# should be deleted as they are not actually there 
SingleSNPssynonymous <-subset(SingleSNPssynonymous,SingleSNPssynonymous$refbase != SingleSNPssynonymous$changedbase)

aggregatedstatsnonsynom <-aggregate(as.factor(SingleSNPssynonymous$Country),by=list(SingleSNPssynonymous$refbase,SingleSNPssynonymous$changedbase,SingleSNPssynonymous$name),length)

colnames(aggregatedstatsnonsynom) <- c("Georgia reference base","changed base","Sample ID", "Frequency")
aggregatedstatsnonsynom


aggregatedstats <-aggregate(as.factor(SingleSNPsnonsynonymous$Country),by=list(SingleSNPsnonsynonymous$refbase,SingleSNPsnonsynonymous$changedbase,SingleSNPsnonsynonymous$name),length)

colnames(aggregatedstats) <- c("Georgia reference base","changed base","Sample ID", "Frequency")
aggregatedstats


#SingleSNPS_allseqs_commmon_regional <- SingleSNPsnoambig
#SingleSNPS_allseqs_allSNPS <- SingleSNPsnoambig
#SingleSNPS_this_study_commmon_regional <- SingleSNPsnoambig
#SingleSNPS_this_study_all <- SingleSNPsnoambig

setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Final_analyses_Jan2024/Core_genes_for_tree/results")


#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Results/SNP analysis/Gene regions")

write.table(aggregatedstats,file = "Frequency of synonymous base pair changes per sample genes ASFV_mutations_all_seqs_all_SNPs.tsv",sep = "\t",row.names = rownames(aggregatedstats),col.names =(colnames(aggregatedstats)))

write.table(aggregatedstatsnonsynom,file = "Frequency of nonsynonymous base pair changes per sample genes ASFV_mutations_all_seqs_all_SNPs.tsv",sep = "\t",row.names = rownames(aggregatedstatsnonsynom),col.names=(colnames(aggregatedstatsnonsynom)))

write.table(SingleSNPsnoambig,file = "All_SNPS_table_ASFV_gene_regions_all_mutations_all_seqs_all_SNPs.tsv",sep = "\t",row.names = FALSE,col.names=(colnames(SingleSNPsnoambig)))

write.table(SingleSNPsnonsynonymous,file = "All_SNPS_table_ASFV_gene_regions_non_synonymous_mutations_all_seqs_all_SNPs.tsv",sep = "\t",row.names =FALSE,col.names=(colnames(SingleSNPsnoambig)))

write.table(SingleSNPssynonymous,file = "All_SNPS_table_ASFV_gene_regions_synonymous_mutations_all_seqs_all_SNPs.tsv",sep = "\t",row.names = FALSE,col.names=(colnames(SingleSNPsnoambig)))

# At this point I should just count them manually and confirm what the mutations are. There are 106 in total but Many will be different inds for same gene
# and same position


