WD=/media/tem/ssd/skin_v1_v2/vgt-report-skin/

basedir=$WD/gwas_base/
basepop=skin_freckles.hg38
tardir=/home/tem/private_NAS/giang/ASA_imputed/
tarpop=$1
out_dir=$WD/output
db_report=${WD}/db_report

#basedir=$1
#basepop=$2
#tardir=$3
#tarpop=$4
#out_dir=$5
#db_report=$6

#############################################################################
#./clumping
start=$(date +'%s')
./QC_target.sh $basedir/vn1008.valid.snp  $tardir $tarpop $basedir $basepop
#generate_PRS
./generate_PRS.sh $basedir $basepop $tardir $tarpop $out_dir

# Extract SNP for each sample
./sample.extract_2.sh $tardir $tarpop ${db_report} $out_dir

# # Generate_report
Rscript --slave analysis.R $out_dir/PRS_${tarpop} $out_dir ${db_report}/gwas.skin.catalog.rda ${db_report}/db_v3.xlsx "skin" $out_dir/each_sample_snp_folder ${db_report}/trait_ratio.xlsx
echo "It took $(($(date +'%s') - $start)) seconds"
