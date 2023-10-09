
valid_snp=$1
tardir=$2
tarpop=$3
basedir=$4
basepop=$5

echo "query clumped SNP from imputed files $tarpop.vcf.gz!"
#bcftools filter -e 'R2 > 0.7 & TYPED=1' ${tardir}/${tarpop}.vcf.gz -Oz -o ${tardir}/${tarpop}.clumped1.vcf.gz
bcftools filter -i ID==@${valid_snp} ${tardir}/${tarpop}.vcf.gz -Oz -o ${tardir}/${tarpop}.clumped.vcf.gz

plink --vcf ${tardir}/${tarpop}.clumped.vcf.gz \
    --make-bed --out ${tardir}/${tarpop} \
    --silent \
    --double-id

echo "correct ambigious SNPS"
Rscript filter_mismatching_snps.R  $tardir $basedir $basepop

plink --bfile ${tardir}/${tarpop} \
    --a1-allele ${tardir}/${tarpop}.tmp.a1 \
    --make-bed \
    --out ${tardir}/${tarpop}.QC \
    --write-snplist \
    --silent

