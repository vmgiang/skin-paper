basedir=$1
basepop=$2
tardir=$3
tarpop=$4
out_dir=$5

awk '{print $3,$8}' ${basedir}/${basepop}.QC.Transformed > ${tardir}/SNP.pvalue

touch ${tardir}/range_list
for range_list in "0.001 0 0.001" "0.05 0 0.05" "0.1 0 0.1" "0.2 0 0.2" "0.3 0 0.3" "0.4 0 0.4" "0.5 0 0.5"
do
  echo $range_list >> ${tardir}/range_list
done


mkdir ${out_dir}/PRS_${tarpop}

plink \
    --bfile ${tardir}/${tarpop}.QC \
    --score ${basedir}/${basepop}.QC.Transformed 3 4 12 header \
    --q-score-range ${tardir}/range_list ${tardir}/SNP.pvalue \
    --extract ${tardir}/${tarpop}.QC.snplist \
    --out ${out_dir}/PRS_${tarpop}/${tarpop} \
    --silent
rm ${tardir}/${tarpop}.QC* ${tardir}/range_list ${tardir}/SNP.pvalue 
rm ${tardir}/${tarpop}.{bim,bed,log,nosex,fam}
rm ${tardir}/${tarpop}.clumped*.vcf.gz 
rm ${tardir}/*.tmp*
