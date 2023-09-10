wrfp=/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
datapath=/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/06.ID_modify_and_GATK_hardfilter/
for chrn in `cat chr.list`
do
    mkdir -p ../$chrn
    grep $chrn split_region.list | while read line
    do 
        region=`echo $line | awk '{print $1}'`
        id=`echo $line | awk '{print $2}'`

        echo "
        ##nos##
        while :
        do
            if [ -e $datapath/$chrn/HaplotypeCaller.${chrn}.${id}.SNP.HARD.Missing-unphasing.ID.vcf.gz.csi ];then 
        ##nos##
                java -Xmx20g -jar /public/agis/chengshifeng_group/fengcong/WGRS/software/snpEff/snpEff.jar \\
                -v  iwgscv1_HC -reg tRNA -reg TE -reg lncRNA -reg miRNA -reg EPI \\
                $datapath/$chrn/HaplotypeCaller.${chrn}.${id}.SNP.HARD.Missing-unphasing.ID.vcf.gz \\
                | bgzip -c > ../$chrn/HaplotypeCaller.${chrn}.${id}.SNP.HARD.Missing-unphasing.ID.ann.vcf.gz
                
                /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../$chrn/HaplotypeCaller.${chrn}.${id}.SNP.HARD.Missing-unphasing.ID.ann.vcf.gz
                break;
        ##nos##
            else
                sleep 1m;
            fi
        done
        ##nos##
        " > $chrn.$id.SNP_ann.sh 
    done 
done
