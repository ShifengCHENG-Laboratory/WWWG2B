wrfp=/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
for chrn in `cat chr.list`
do
    mkdir -p ../$chrn
    gvcflist=""
    for sample in `cat sample.list`
    do
         gvcflist="${gvcflist} --variant ../../02.all_gvcf/$sample/${sample}.HaplotypeCaller.${chrn}.g.vcf.gz "

    done
    grep $chrn split_region.list | while read line
    do 
        region=`echo $line | awk '{print $1}'`
        id=`echo $line | awk '{print $2}'`

        echo "
        /public/agis/chengshifeng_group/fengcong/WGRS/software/gatk-4.1.2.0/gatk --java-options "-XX:ParallelGCThreads=5" CombineGVCFs \
        -R $wrfp \
        $gvcflist \
        -L $region \
        -O ../$chrn/HaplotypeCaller.${chrn}.${id}.g.vcf.gz

        /public/agis/chengshifeng_group/fengcong/WGRS/software/gatk-4.1.2.0/gatk --java-options "-XX:ParallelGCThreads=5" GenotypeGVCFs \
        -R $wrfp \
        -V ../$chrn/HaplotypeCaller.${chrn}.${id}.g.vcf.gz \
        -O ../$chrn/HaplotypeCaller.${chrn}.${id}.vcf.gz
        " > $chrn.$id.jointcall.sh 
    done 
done
