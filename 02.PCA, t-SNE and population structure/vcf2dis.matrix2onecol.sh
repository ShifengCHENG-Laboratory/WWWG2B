#!/bin/sh
cd /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis
cat /vol3/agis/chengshifeng_group/chenghong/wheat/04.diversity/01.pi/03.group_1047/AG1.chr1A.1M.windowed.pi |while read line
do
chr=`echo $line | awk '{print $1}'`
start=`echo $line | awk '{print $2}'`
end=`echo $line | awk '{print $3}'`
perl /vol3/agis/chengshifeng_group/chenghong/script/vcf2dis.matrix2onecol.pl vcf2dis.samplelist ./$chr/${chr}-${start}-${end}.p_dis.mat ./$chr/${chr}-${start}-${end}.p_dis.mat.trans
done
