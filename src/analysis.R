#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args[1] <- '/home/tem/Desktop/code/skin/new_branch/vgt-report-skin/output/PRS_ASA_part1/'
setwd(args[1])
# args[7] <- '/home/tem/Desktop/code/skin/database/trait_ratio.xlsx'
##install packages
# list.of.packages <- c("ggplot2","stringr","readxl","magrittr","data.table","doParallel")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
norm <- function(x){(x-min(x))/(max(x)-min(x))}
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
filenames <- list.files(args[1], pattern = "profile")

personFIDs <- as.character(read.table(filenames[1], header=T,colClasses = (c(rep('character',2),rep('numeric',4))))$FID)
# personIDs <- unlist(lapply(personFIDs, function(x) { strsplit(as.character(x), "_", fixed = TRUE)[[1]][2] }))
personIDs <- personFIDs
score.all <- lapply(personFIDs, function(x) {
  unlist(lapply(filenames, function(f){
    prs <- read.table(f, header=T,colClasses = (c(rep('character',2),rep('numeric',4))))
    prs.norm <- norm(prs$SCORE)
    names(prs.norm) <- prs$FID
    prs.norm[x]
  }))
})

names(score.all) <- personIDs
score.all.mean <- unlist(lapply(score.all, function(x) mean(x)))

##raw report generation
library(readxl)
#args[2] <- '/home/tem/private_NAS/giang/PRS_205540230122.imputed_old/'
skin_folder <- args[2]
setwd(skin_folder)

##summarized by traits
library(stringr)
skin_traits <- c("UV protection","Moisturizing",
                 "Inflamatory Cytokine",
                 "Antioxidation response",
                 "Elasticity","Skin sensitivity","Glycation",
                 "Skin aging","Sagging eyelid",
                 "Acne","Wrinkle","Collagen degradation","Stretch mark",
                 "Tanning response","Freckle","Vitamin B9 need",
                 "Vitamin A need","Vitamin B2 need","Vitamin B6 need",
                 "Vitamin B12 need","Vitamin C need",
                 "Vitamin D need","Vitamin E need","Vitamin K need",
                 "Omega need")
# skin_traits <- str_to_title(skin_traits)

#args[3] <-'/home/tem/Desktop/code/skin/new_branch/vgt-report-skin/db_report/gwas.skin.catalog.rda'
load(args[3])

context.sig <- names(table(as.character(gwas.skin.catalog$CONTEXT)))
context.sig <- context.sig[context.sig!=""]
context.sig.map <- c("stop_gained",
                     "missense_variant",
                     "synonymous_variant",
                     "5_prime_UTR_variant",
                     "3_prime_UTR_variant",
                     "intron_variant",
                     "non_coding_transcript_exon_variant",
                     "regulatory_region_variant",
                     "TF_binding_site_variant",
                     "intergenic_variant")
context.sig.map <- cbind(context.sig.map, rev(seq(1:length(context.sig.map))))
colnames(context.sig.map) <- c("name","rank")
top1 <- c("stop_gained")
top2 <- c(top1, "missense_variant")
top3 <- c(top2, "synonymous_variant")
top4 <- c(top3, "5_prime_UTR_variant")
top5 <- c(top4, "3_prime_UTR_variant")
topk_thr <- top5

##curated DB
#args[4] <- '/home/tem/Desktop/code/skin/new_branch/vgt-report-skin/db_report/db_v3.xlsx'
db <- read_excel(args[4],sheet = 1)
db$strongest_snp <- unlist(lapply(db$SNP,function(x)strsplit(x,"-")[[1]][1]))
topTraits <- sort(table(db$trait))
names(topTraits) <- str_to_title(names(topTraits))
#args[5] <- 'skin'
diseaseID <- as.character(args[5])
#args[6] <- '/home/tem/private_NAS/giang/each_sample_snp_folder/'
snp_folder <- args[6]
ratio <- read_xlsx(args[7])
# library(doParallel)
# cl <- makeCluster(8)
# registerDoParallel(cl)
# 
# print("done")
# foreach(i = 1:length(personIDs)) %dopar% {

for (i in 1:length(personIDs)) {  
  personID <- personIDs[i]
  dir.create("raw_report", showWarnings = FALSE)
  dir.create(paste0("raw_report/",personID), showWarnings = FALSE)
  mean.all <- mean(score.all.mean)
  sd.all <- sd(score.all.mean)
  score <- mean(score.all[[personID]])
  sd <- sd(score.all[[personID]])
  
  if (score > (mean.all+2*sd.all)) {
    main_text <- "high risk"
    main_text_vn <- "nguy cơ cao"
    color <- "red"
  } else if (score > (mean.all+sd.all)){
    main_text <- "increased risk"
    main_text_vn <- "tăng nguy cơ"
    color <- "orange"
  } else if (score < (mean.all-sd.all)){
    main_text <- "low risk"
    main_text_vn <- "nguy cơ thấp"
    color <- "blue"
  } else {
    main_text <- "medium risk"
    main_text_vn <- "nguy cơ trung bình"
    color <- "orange"
  }
  dens <- density(score.all.mean, kernel = "gaussian")
  
  x1 <- min(which(dens$x >= (score-sd)))
  x2 <- max(which(dens$x <  (score+sd)))
  plot_points <- c(score.all.mean, score, sd)
  write.csv(round(plot_points,digits = 2), file = paste0(skin_folder,"/raw_report/",personID,"/histogram_points_",personID,".csv"), row.names = F)
  
  pdf(paste0(skin_folder,"/raw_report/",personID,"/histogram_",personID,".pdf"), height = 5, width = 5)
  plot(dens,
       main = paste0("n = ", length(score.all)),
       xlab = "Polygenic risk score", xlim = c(0,1))
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col = color))
  dev.off()
  
  ##truncated two sides
  if(score-sd < 0) {
    lowerbound <- 0
  } else {
    lowerbound <- score-sd
  }
  if(score+sd.all > 1) {
    upperbound <- 1
  } else {
    upperbound <- score+sd
  }
  
  text <- cbind(
    personID,
    diseaseID,
    main_text,
    paste0(round(score*100, digits = 1)," %"),
    paste0(round(lowerbound*100, digits = 1)," %"),
    paste0(round(upperbound*100, digits=1)," %"),
    paste0("The person with ID of ", personID,
           " has ", main_text,
           " (averge ", round(score*100, digits = 1),
           "%, from ", round(lowerbound*100, digits = 1),
           "% to ", round(upperbound*100, digits=1),"%)"),
    paste0("Phân tích kiểu gen của ", personID,
           "cho thấy nguy cơ mắc bệnh về da thẩm mỹ ở mức " ,
           main_text_vn," (trung bình " ,
           round(score*100, digits = 1)," % từ " ,round(lowerbound*100, digits = 1),
           " % đến ",round(upperbound*100, digits=1),"%)" ))
  colnames(text)[3:8] <- c("association_risk","average","lowerbound","upperbound","explanation_vn","explanation_en")
  write.csv(text, file = paste0(skin_folder,"/raw_report/",personID,"/summary_",personID,".csv"), row.names = F)
  
  
  filenames <- list.files(args[6])
  file_personID <- personFIDs[i]
  snps_sample <- data.table::fread(paste0(snp_folder,'/',filenames[grepl(file_personID, filenames)]),header = F)$V1

  skin.snp.tbl <- db[grepl("^skin", db$trait, ignore.case = T),]
  skin.snp.rsIDs <- skin.snp.tbl$strongest_snp
  nskin <- length(skin.snp.rsIDs)
  traitw <- 0.7
  skinw <- 1-traitw

  traitM <- do.call("rbind",lapply(skin_traits,function(x){
    trait.snp.tbl <- db[grepl(x, db$trait, ignore.case = T),]
    #print(trait.snp.tbl)
    trait.snp.rsIDs <- trait.snp.tbl$strongest_snp
    ntrait <- length(trait.snp.rsIDs)
    nt <- 0
    ns <- 0
    if(ntrait>0&nskin>0){
      # tt <- intersect(snps_sample, trait.snp.rsIDs)
      tt <- snps_sample[snps_sample %in% trait.snp.rsIDs]
      # ss <- intersect(snps_sample, skin.snp.rsIDs)
      ss <- snps_sample[snps_sample %in% skin.snp.rsIDs]
      nt <- ifelse(length(tt)>0,length(tt),0)
      ns <- ifelse(length(ss)>0,length(ss),0)
      sc <- round(nt/(2*ntrait)*traitw+ns/(2*nskin)*skinw, digits = 2)*100
      r <- ifelse((sc >=50 & ratio[ratio$Trait == x,]$type == 0) | (sc < 50 & ratio[ratio$Trait == x,]$type == 1),
                  0,ratio[ratio$Trait == x,]$ratio)
      paste(sc,r,sep=',')
    } else if (nskin>0) {
      sc <- round(ns/nskin*skinw, digits = 2)*100
      r <- ifelse((sc >=50 & ratio[ratio$Trait == x,]$type == 0) | (sc < 50 & ratio[ratio$Trait == x,]$type == 1),
                  0,ratio[ratio$Trait == x,]$ratio)
      paste(sc,r,sep=',')
    } else {
      print("Not detected!")
    }
  }))
  colnames(traitM) <- "trait_skin,ratio"
  rownames(traitM) <- skin_traits
  write.csv(traitM, file = paste0(skin_folder,"/raw_report/",personID,"/traitM_",personID,".csv"),quote = F)

  traitDetailM <- do.call("rbind",lapply(skin_traits,function(x){
    trait.snp.tbl <- db[grepl(x, db$trait, ignore.case = T),]
    trait.snp.rsIDs <- unlist(lapply(trait.snp.tbl$SNP,function(s) strsplit(s,"-")[[1]][1]))
    ntrait <- length(trait.snp.rsIDs)
    tt <- intersect(snps_sample, trait.snp.rsIDs)
    ss <- intersect(snps_sample, skin.snp.rsIDs)
    if(length(tt)>0&length(ss)>0){
      rest <- trait.snp.tbl[trait.snp.tbl$strongest_snp%in%tt,]
      return(rest)
    } else if(length(ss)>0){
      ress <- skin.snp.tbl[skin.snp.tbl$strongest_snp%in%ss,]
      ress <- ress[sample(rownames(ress),2),]
      ress$trait <- x
      return(ress)
    } else {
      print("Not detected!")
    }
  }))
  write.csv(traitDetailM, file = paste0(skin_folder,"/raw_report/",personID,"/traitDetailM_",personID,".csv"),quote = F)
}
# stopCluster(cl)



