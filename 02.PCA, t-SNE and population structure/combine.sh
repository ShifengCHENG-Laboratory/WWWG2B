#!/bin/sh
dir=/vol3/agis/chengshifeng_group/chenghong/wheat/01.SNPdatasets/00.input/watkins/1MB_window
cd /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis
awk '$1 !~ /CHROM/' /vol3/agis/chengshifeng_group/chenghong/wheat/04.diversity/01.pi/03.group_1047/AG1.chr1A.1M.windowed.pi |cat |while read line
do
chr=`echo $line | awk '{print $1}'`
start=`echo $line | awk '{print $2}'`
end=`echo $line | awk '{print $3}'`
awk '{print $3}' ./$chr/${chr}-${start}-${end}.p_dis.mat.trans > ./$chr/a
paste ./$chr/a b > ./$chr/${chr}-${start}-${end}.cor.input
done
