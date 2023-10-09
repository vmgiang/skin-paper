tardir=$1
tarpop=$2
db_report=$3
out_dir=$4
vcf_file=${tardir}/${tarpop}.vcf.gz
snp_list=${db_report}/db_snp_report.txt
mapping_rsID=${db_report}/mapping_rsID.txt

echo "extract snp for each sample from $vcf_file"

#Extract sample name 
bcftools query -l ${vcf_file} > ${out_dir}/sample.txt
(echo "SNP" & cat ${out_dir}/sample.txt) |\
    tr '\n' ',' |\
    sed 's/.$/\n/' > ${out_dir}/sample2.txt
bcftools filter -i ID==@$snp_list $vcf_file |\
    bcftools query -f '%ID[,%DS]\n' > ${vcf_file/.vcf.gz/.dosage.tmp.txt}

cat ${out_dir}/sample2.txt ${vcf_file/.vcf.gz/.dosage.tmp.txt} > ${out_dir}/${tarpop}.dosage.txt

rm ${vcf_file/.vcf.gz/.dosage.tmp.txt} ${out_dir}/sample*.txt

mkdir -p $out_dir/each_sample_snp_folder
Rscript sample.extract.snp.R ${out_dir}/${tarpop}.dosage.txt ${mapping_rsID} $out_dir/each_sample_snp_folder