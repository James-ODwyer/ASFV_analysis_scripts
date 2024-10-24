# Analysis of gene mutations and positions of mutations in ASFV 

library(seqinr)
library(ggtree)
library(plyr)
library(dplyr)
library(ape)
library(Biostrings)



#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Results/denovo_multistep/Rago_assembled/final_trees/2ndassembly/genes/FNA/2")

setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/Asia_pacific_refs/Genes/MAFFT/FNA/alignments")

fastas <- list.files(path =getwd(),pattern = "fas")
snplist <- list()
for (i in c( 1:length(fastas))) {
  
  
  FNAdata <- read.fasta(file= fastas[i],seqtype = "DNA", as.string = FALSE)
  
  FNAdata[] <- lapply(FNAdata, c)
  
  FNAdata2 <- bind_rows(FNAdata,.id = "Sampleid")
  
  
  snps <- as.data.frame(FNAdata2)
  names <- colnames(snps)
  gapChar <- "n"
  snps %>% relocate("FR682468.2") -> snps
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
# Didn't really accomplish anything but there are now 164 heatmaps




#1.  Define all SNPs 
#2. Add metadata to the SNP file
#3. Sort/filter to see SNPs found in only 1 country 


# 1.

# FNA files but with Georgia root defining SNPs  
setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/Final runs/Asia_pacific_refs/Genes/MAFFT/FNA/2")

fastas <- list.files(path =getwd(),pattern = "fas")
snplist <- list()
for (i in c( 1:length(fastas))) {
  
  
  FNAdata <- seqinr::read.fasta(file= fastas[i],seqtype = "DNA", as.string = FALSE)
  
  FNAdata[] <- lapply(FNAdata, c)
  
  FNAdata2 <- bind_rows(FNAdata,.id = "Sampleid")
  
  
  snps <- as.data.frame(FNAdata2)
  names <- colnames(snps)
  gapChar <- "n"
  snps %>% relocate("FR682468.2") -> snps
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
    #names(snplist[[i]]) <- fastas[i]
  }
  
  if (nrow(snp_data) ==0) {
    
    snplist[[i]] <- NULL
  }
  
  if (i == length(fastas)) {
    
    if (nrow(snp_data) ==0) {
      
      snplist[[i]] <- "NULL"
      
      
    }
    
  }
  
}



if(snplist[[i]]!="NULL") {
  
  names(snplist[1:184]) <- fastas[1:184]
}


if (snplist[[i]]=="NULL") {
  names(snplist) <- fastas
  snplist[[i]] <- NULL
  
}


snplist



snplist2 <-bind_rows(snplist, .id = "gene")



snplist2$gene <- sub(pattern = ".fas",replacement = "",x = snplist2$gene)

snplist2$gene <- sub(pattern = "_fixed",replacement = "",x = snplist2$gene)

snplist2$gene <- gsub("^[0-9]*_",replacement = "",x = snplist2$gene)
snplist2$gene <- gsub("_CDS",replacement = "",x = snplist2$gene)
snplist2$gene <- gsub(" CDS",replacement = "",x = snplist2$gene)

geneandpos <-read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Results/denovo_multistep/Rago_assembled/final_trees/FR682468 Annotations.csv",header=FALSE)

geneandpos$V3 <- gsub(" ",replacement = "_",x = geneandpos$V3)

snplist2$pos <- as.numeric(snplist2$pos)
snplist3 <- snplist2
for (i in c(1:nrow(snplist2))) {
  
  grep(paste0("^",snplist2$gene[i]),x = geneandpos$V3) -> idx
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

metadata <- read.csv("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/phylogenetic trees/metadataASFV52samples.csv", header = TRUE)


snp_databothpos <- snp_data
snp_databothpos$posgene <- snplist2$pos
# For full genome do snp_data$pos
#snp_databothpos$posgene <- snp_data$pos
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

# So 71 substitutions (not unique positions but positions multiplied by samples) were common (found in over 50% of samples from that country)
# 17,059 substitutions (98.1%) were found across more than 5 countries and are widespread 
# 5 substitutions were found across multiple but fewer than 5 countries.
# 84 substitutions were found in a single country 

#country_specific_SNPs <- subset(snp_data_meta,snp_data_meta$SNPtype=="Country specific")

aggregatedstats <-aggregate(as.factor(snp_data_meta$Country),by=list(snp_data_meta$SNPtype,snp_data_meta$Country),length)


SingleSNPs <- subset(snp_data_meta,snp_data_meta$SNPtype %in% c("Common single country", "Present single country","fixed single country") )
SingleSNPsalllocals <- subset(snp_data_meta,snp_data_meta$SNPtype %in% c("Common single country", "Present single country","fixed single country") )
SingleSNPs <- subset(SingleSNPs,SingleSNPs$Study== "This study" )

# removing ambig code Indian ref seq 

countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="PNG")
countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="Vietnam")


aggregatedstats <-aggregate(as.factor(countrysingles$Country),by=list(countrysingles$SNPtype,countrysingles$pos),length)

aggregatedstats

PNGsingles <- subset(SingleSNPs,SingleSNPs$Country=="PNG")

unique(PNGsingles$pos)
unique(PNGsingles$gene)
# 19 substitution sites across 6 genes


EastTimorsingles <- subset(SingleSNPs,SingleSNPs$Country=="East Timor")

unique(EastTimorsingles$pos)
unique(EastTimorsingles$gene)
# Only two genomes for East Timor so can a statement about fixation be made 
# 1 substitution sites in one genes was found in both East Timor samples. 


Vietnamsingles <- subset(SingleSNPs,SingleSNPs$Country=="Vietnam")

unique(Vietnamsingles$pos)
unique(Vietnamsingles$gene)

# 16 sites within the one gene MGF_360_21R were found in two Vietnam samples sequenced.

# Have done the three main countries where samples for this study were taken from 
# Can do all but discuss with David and MAtt
# 

# Next, need to get the gene specific SNP position for each important substitution and translate to AA to see whether they are synonymous vs nonsynonymous
# Added gene position into above code. Now will re extract relevant base from fastas
# 

SingleSNPs$refbase <- NA
SingleSNPs$refbasesequence <- NA
SingleSNPs$changedbase <- NA
SingleSNPs$sequence <- NA

#setwd("C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/ASFV/Results/denovo_multistep/Rago_assembled/final_trees/2ndassembly/genes/FNA/2")

fastas <- list.files(path =getwd(),pattern = "mafft")

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
    
      
      if (sum(grep("m|r|y|k|b|d|v|h|w",snps[l,])) >=1) {
        grep("m|r|y|k|b|d|v|h|w",snps[l,]) -> idx
        baseposns <- snps[l,]
        baserefns <- apply(baseposns[], 1, function(x) {names(which.max(table(factor(x,unique(x)))))})
        
        snps[l,idx] <- baserefns
        
      }
      
    
  }
  
  gene <- sub(pattern = ".mafft",replacement = "",x = fastas[i])
  
  gene <- sub(pattern = "_fixed",replacement = "",x = gene )
  
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
    
    SingleSNPs[idx[j],8] <- baseref
    SingleSNPs[idx[j],9] <- baserefsequence
    SingleSNPs[idx[j],10] <- altbase
    SingleSNPs[idx[j],11] <- sequence
    
    
    }
  
  }  
    
  
  

  
}

rownames(SingleSNPs) <- c(1:nrow(SingleSNPs))







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

# Now picking up from the end of generating single SNPs to now generate complete AA codes 

SingleSNPs$refseqshort <- NA 
SingleSNPs$altseqshort <- NA 
SingleSNPs$refAA <- NA 
SingleSNPs$altAA <- NA 
SingleSNPs$identical <- NA 
SingleSNPs$SNPchainlengthref <- NA 
SingleSNPs$SNPchainlengthalt <- NA 
for ( i in c(1:nrow(SingleSNPs))) {
  
  # to extract the relevant codon will need to identify what bp position the mutation is and extract around it
  # divide by 3 and whole number = end base of codon (add prior two bases)
  # divide by 3 and x.67 middle base of codon (add one prior and one after base)
  # divide by three and x.33 start of codon (add two after)
  # 
  SingleSNPs$refseqshort[i] <- gsub(pattern = "-",replacement = "",x = SingleSNPs$refbasesequence[i])
  SingleSNPs$altseqshort[i] <- gsub(pattern = "-",replacement = "",x = SingleSNPs$sequence[i])
  
  AAstringref <- translate(DNAStringSet(SingleSNPs$refseqshort[i]))
  AAstringalt <- translate(DNAStringSet(SingleSNPs$altseqshort[i]))
  
  SingleSNPs$refAA[i] <- as.character(unlist(AAstringref))
  SingleSNPs$altAA[i] <- as.character(unlist(AAstringalt)) 
  
  SingleSNPs$identical[i] <- identical(x = SingleSNPs$refAA[i],SingleSNPs$altAA[i])
  
  
  
}


SingleSNPssynonymous <- subset(SingleSNPs,SingleSNPs$identical==TRUE)
SingleSNPsnonsynonymous <- subset(SingleSNPs,SingleSNPs$identical==FALSE)

SingleSNPssynonymousthisstudy <- subset(SingleSNPssynonymous,SingleSNPssynonymous$Study=="This study")
SingleSNPsnonsynonymousthisstudy <- subset(SingleSNPsnonsynonymous,SingleSNPsnonsynonymous$Study=="This study")
# At this point I should just count them manually and confirm what the mutations are. There are 106 in total but Many will be different inds for same gene
# and same position


countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="PNG")
countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="Vietnam")
countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="China")
countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="South Korea")
countrysingles <- subset(SingleSNPs,SingleSNPs$Country=="India")


aggregatedstats <-aggregate(as.factor(countrysingles$Country),by=list(countrysingles$SNPtype,countrysingles$pos,countrysingles$identical),length)

aggregatedstats

