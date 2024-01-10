#!/bin/sh
cd /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis
awk '$1 !~ /CHROM/' /vol3/agis/chengshifeng_group/chenghong/wheat/04.diversity/01.pi/03.group_1047/AG1.chr1A.1M.windowed.pi |cat |while read line
do
chr=`echo $line | awk '{print $1}'`
start=`echo $line | awk '{print $2}'`
end=`echo $line | awk '{print $3}'`
/public/home/jiangmei/miniconda3/envs/R4.0.3/bin/Rscript /public/home/jiangmei/test/cor.R ./$chr/${chr}-${start}-${end}.cor.input ./$chr/${chr}-${start}-${end}.cor.output
#tail -n 1 ./$chr/${chr}-${start}-${end}.cor.output |awk '{print $1}' >> ${chr}.cor.merge
done
