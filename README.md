This script is to generate a raw report folder ready for System group.


Requirements:

* R version: R version 4.1.1
* Packages for R: 
	* data.table 
	* dplyr
* Plink 1.9 : https://www.cog-genomics.org/plink/

INPUT:
	
	* Imputed vcf file

OUTPUT:
	
	* Raw reports for all samples in vcf files : 

		* histogram*.pdf
		
		* histogram_points*.pdf
	
		* summary*.csv
	
		* traitDetailM*.csv

		* traitM*.csv   

CHECKING Quality of Raw report:

	* Rscript quality_check.R path_to_raw_report
	

RUNNING STEPS:

1. Edit parameters in src/parameters.sh:
	
	* MAF_thres : MAF threshold
	
	* IMP_thres : Imputation threshold

	* hwe : Hardy-weiberg

2. Run pipeline for generating report :
	
	./run_pipline.sh basedir_path basepop_prefix tardir_path tarpop_prefix outdir_path db_report	
	
	* basedir_path : GWAS summary statistic, Clumped SNPs from vn1008
	
	* basepop : prefix Gwas sumstat
	
	* tardir_path : vcf_file directory contains imputed VCF file
	
	* tarpop_prefix : prefix of vcf files ($tarpop_prefix.vcf.gz)

	* db_report : db for generating report

	

	

