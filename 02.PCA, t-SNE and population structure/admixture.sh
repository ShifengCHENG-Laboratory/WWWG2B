#!/bin/sh
# generate random number: shuf -i 1-100 -n 1
export K=$1
export seed=$2
echo "#!/bin/sh
admixture -s $seed --cv /vol3/agis/chengshifeng_group/chenghong/wheat/01.SNPdatasets/02.wgs_exon/Bgenome.maf001.missing02.bed $K -j4 | tee log${K}.out
" > ${K}.admixture.sh

