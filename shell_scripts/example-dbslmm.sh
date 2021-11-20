### change file permission for dbslmm
#DIR=~/tmp/DBSLMM
DIR=~/research/DBSLMM
dbslmm=../bin/dbslmminterval
#dbslmm=${DIR}/scr/dbslmm
### Parameters for DBSLMM
let chr=1
#for chr in $(seq 1 22)
#do 
echo ${chr}
DBSLMM=${DIR}/software/DBSLMM.R
summf=${DIR}/test_dat/summary_gemma_chr
outPath=${DIR}/test_dat/out/
plink=/usr/cluster/bin/plink-1.9
ref=${DIR}/test_dat/ref_chr
outfile=armafieldChr${chr}.dat

blockf=${DIR}/block_data/EUR/chr
m=`cat ${summf}${chr}.assoc.txt | wc -l` 
h2=0.5 # need to replace with input from ldscore output
nobs=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
# write a row to control file
#declare -A my_array # -A means associative array: https://linuxconfig.org/how-to-use-arrays-in-bash-script
declare -a my_array2
#my_array=([chr] = "${chr}" [block] = ${blockf}${chr}.bed [n]=${n} [nsnp]=${m} \
#    [outfile] = ${outfile} [model] = DBSLMM [type] = auto [ref] = ${ref}${chr} [dbslmm] = ${dbslmm} \
#    [h2] = ${h2} [summary] = ${summf}${chr}.assoc.txt [outpath] = ${outPath})
my_array2=(${chr} \
            # b 
            ${blockf}${chr}.bed \
            # eff 
            ${outPath}${chr}.dbslmm \
            # heritability
            ${h2} \
            # l, large effects summary file path
            ${outPath}${chr}l.txt \
            # mafmax default 0.2
            0.2 \
            # n sample size
            ${n} \
            # nsnp
            ${m} \
            # r, bfile for ref data
            ${DIR}/test_dat/ref_chr${chr} \
            # s small effects summary data
            ${outPath}${chr}s.txt \
            # type
            auto 
            # 
)
join_arr() {
  local IFS="$1"
  shift
  echo "$*"
}

join_arr , "${my_array2[@]}" >> control.csv # https://linuxize.com/post/bash-append-to-file/
#done

## execute Rscript
Rscript ${DBSLMM} --summary ${summf}${chr}.assoc.txt --outPath ${outPath} \
  --plink ${plink} --dbslmm ${dbslmm} --outfile ${outfile} --ref ${ref}${chr} --n ${n} \
  --nsnp ${m} --type auto --model DBSLMM --block ${blockf}${chr}.bed \
  --h2 ${h2}
  
  
  
for chr in seq `1 22` 
do 
  
### Predict
bfilete=${DIR}/test_dat/test_chr
est=${DIR}/test_dat/out/summary_gemma_chr
InterPred=${DIR}/test_dat/out/internal_pred_chr
## plink 1.9
plink=/usr/cluster/bin/plink-1.9
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 sum --out ${InterPred}${chr}
## plink 2
plink=plink2
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 cols=+scoresums --out ${InterPred}${chr}

done
