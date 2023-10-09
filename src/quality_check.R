suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = T)

raw_report <- args[1]
print(paste("raw_report on :",raw_report,sep = ''))
setwd(raw_report)
sampleL <- list.files(raw_report)
print(paste("Number of sample: ",length(sampleL)))
check_sample <- function(sample){
 out_file <- list.files(sample)
 histogram_points__colname <- c('x')
 traitDetail_colname <- c("V1","trait","GENE","SNP","REF","ALT","CHR","AUTHOR","PUBMEDID","CONTEXT","strongest_snp")
 files <- c(paste0('histogram_',sample,'.pdf'),
            paste0('histogram_points_',sample,'.csv'),
            paste0('summary_',sample,'.csv'),
            paste0('traitDetailM_',sample,'.csv'),
            paste0('traitM_',sample,'.csv'))
 if (sum(!(files %in% out_file)) == 0 ){
  list_table <- c('histogram_points_','summary_','traitM_')
  check_list <- lapply(list_table, function(x){
    file <- fread(paste0(sample,'/',x,sample,'.csv'),header = T)
    if (sum(is.na(file)) == 0) {
      status <- 'OK'
    }
  })
  file <- fread(paste0(sample,'/','traitDetailM_',sample,'.csv'),header = T)
  if (sum(!(traitDetail_colname %in% colnames(file)))==0){
    status <-'OK'
  } else {status <- 'error'}

  check_list[[4]] <- status
  if (all(check_list=='OK')){
    out <- T
  } else {out <- F}
 } else { out <- F}
 return(out)
}

result_check <- lapply(sampleL, function(x){
  check_sample(x)
})
names(result_check) <- sampleL
if (all(result_check == T)) {
  print("Testing output's format: PASS")
} else {
  print(sampleL[result_check == F])
}



