suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = T)
batch.dosage <- args[1]
mapping.file <- args[2]
out_path <- args[3]

df_dosage <- fread(batch.dosage,header =T)
# print(head(df_dosage,5))
sample_name <- colnames(df_dosage)[!(colnames(df_dosage) %in% 'SNP')]
map <- fread(mapping.file)
# print(head(map))
df_dosage <- merge(df_dosage,map,by='SNP',all.x = T)
sample <- lapply(sample_name,function(x){
  sample_snp <- ifelse(round(df_dosage[[x]]) == 0,"",
                       ifelse(round(df_dosage[[x]]) == 1,
                              df_dosage$rsID,
                              paste(df_dosage$rsID,df_dosage$rsID,sep=',')))
  sample_snp <- sample_snp[sample_snp != ""]
  sample_snp <- unlist(strsplit(sample_snp,','))
  write.table(sample_snp,paste0(out_path,'/',x,'.rsID.txt'),sep = "\n",col.names = F,
              row.names = F, quote = F)
})
