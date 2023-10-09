suppressPackageStartupMessages(library(data.table))
options(warn = -1)

args <- commandArgs(trailingOnly = T)
tardir <- args[1]
basedir <- args[2]
basepop <-args[3] 


# tardir <- '/media/tem/ssd/skin-code/C+T/db_report/bfile'
# basedir <- '/media/tem/ssd/skin-code/C+T/gwas_base/'
# basepop <- 'skin_freckles.hg38'
setwd(tardir)

pheno <- read.table(gzfile(paste(basedir,'/',basepop,".QC.gz",sep='')),header = T,stringsAsFactors = F, sep="\t")
colnames(pheno) <- c('CHR','BP','SNP','A1','A2','N','SE','P','OR','INFO','MAF')
pheno$beta <-log(pheno$OR)


# Change all alleles to upper case for easy comparison
pheno$A1 <- toupper(pheno$A1)
pheno$A2 <- toupper(pheno$A2)
pheno$SNP <- paste(paste0('chr',pheno$CHR),pheno$BP,pheno$A2,pheno$A1,sep = ':')
write.table(pheno,
            paste(basedir,'/',basepop,".QC.Transformed",sep=""),
            col.names = T, row.names = F, quote = F,sep ='\t')
tarpop.list <- list.files()[grepl('bim',list.files())]
tarpop.list <- gsub('.bim','',tarpop.list)
tarpop.list

# Read in bim file
lapply(tarpop.list, function(tarpop){
    bim <- fread(paste(tarpop,".bim",sep=''))
    colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
    
    # Read in QCed SNPs
    # qc <- fread(paste(tarpop,".QC.snplist",sep=''), header = F, stringsAsFactors = F)
    
    # Read in the GWAS data
    bim$B.A1 <- toupper(bim$B.A1)
    bim$B.A2 <- toupper(bim$B.A2)
    pheno$CHR <- as.numeric(pheno$CHR)
    # Merge summary statistic with target
    info <- merge(bim, pheno, by = c("SNP", "CHR", "BP"))
    # Filter QCed SNPs
    # info <- info[info$SNP %in% qc$V1,]
    # Function for finding the complementary allele
    complement <- function(x) {
        switch (
            x,
            "A" = "T",
            "C" = "G",
            "T" = "A",
            "G" = "C",
            return(NA)
        )
    }
    # Get SNPs that have the same alleles across base and target
    info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
    # Identify SNPs that are complementary between base and target
    info$C.A1 <- sapply(info$B.A1, complement)
    info$C.A2 <- sapply(info$B.A2, complement)
    info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
    # Update the complementary alleles in the bim file
    # This allow us to match the allele in subsequent analysis
    complement.snps <- bim$SNP %in% info.complement$SNP
    bim[complement.snps,]$B.A1 <-
        sapply(bim[complement.snps,]$B.A1, complement)
    bim[complement.snps,]$B.A2 <-
        sapply(bim[complement.snps,]$B.A2, complement)
    
    # identify SNPs that need recoding
    info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
    # Update the recode SNPs
    recode.snps <- bim$SNP %in% info.recode$SNP
    tmp <- bim[recode.snps,]$B.A1
    bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
    bim[recode.snps,]$B.A2 <- tmp
    
    # identify SNPs that need recoding & complement
    info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
    # Update the recode + strand flip SNPs
    com.snps <- bim$SNP %in% info.crecode$SNP
    tmp <- bim[com.snps,]$B.A1
    bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
    
    bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
    
    # Output updated bim file
    write.table(
        bim[,c("SNP", "B.A1")],
        paste(tarpop,".tmp.a1",sep=''),
        quote = F,
        row.names = F,
        col.names = F,
        sep="\t"
    )
    
    mismatch <- bim$SNP[!(bim$SNP %in% info.match$SNP | bim$SNP %in% info.complement$SNP | bim$SNP %in% info.recode$SNP | bim$SNP %in% info.crecode$SNP)]
    write.table(
        mismatch,
        paste(tarpop,".tmp.mismatch",sep=''),
        quote = F,
        row.names = F,
        col.names = F
    )
    print(paste0(tarpop," DONE!"))
})




