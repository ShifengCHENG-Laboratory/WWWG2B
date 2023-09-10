wrfp="/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta"
for sample in `cat zz.sample_list.txt`
do
    #sample=`echo $i|cut -f 1 -d ","`
    mkdir -p ../$sample
    > $sample.go.sh
    for library in `grep $sample sample_lib.list | cut -f 2 -d ","`
    do
        #library=`echo $i|cut -f 2 -d ","`
        RGID=$library
        datapath=/vol3/agis/chengshifeng_group/chenglab/01.wheat/12.sup_bgi/zz.all_clean_data_89_rename/
        echo "
        wrfp=\"/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta\"
        sample=$sample
        TMP_DIR=`pwd`/tmp
        bwa mem -t 10  -R \"@RG\\tID:$RGID\\tPL:ILLUMINA\\tLB:$library\\tSM:$sample\" $wrfp $datapath/$sample/${sample}_${library}.R1.fastp.fastq.gz $datapath/$sample/${sample}_${library}.R2.fastp.fastq.gz | samtools sort -T ${sample}.${library} -m 1G -@ 4 - > ../${sample}/${sample}.${library}.sorted.bam

        java -Xmx10g -XX:ParallelGCThreads=10  -Djava.io.tmpdir=`pwd`/tmp -jar /public/agis/chengshifeng_group/fengcong/WGRS/software/picard/picard.jar MarkDuplicates \
        I=../${sample}/${sample}.${library}.sorted.bam \
        M=../${sample}/${sample}.${library}.markdup_metrics.txt \
        O=../${sample}/${sample}.${library}.sorted.markdup.bam && rm ../${sample}/${sample}.${library}.sorted.bam

        samtools index ../${sample}/${sample}.${library}.sorted.markdup.bam

        " >> $sample.go.sh

        

    done

    echo "
    ########merge#######
    samtools merge  ../${sample}/merge.rmdup.${sample}.bam  ../${sample}/*sorted.markdup.bam && rm ../${sample}/*sorted.markdup.bam*
    samtools view -h -T $wrfp -C -@ 4 ../${sample}/merge.rmdup.${sample}.bam -o ../${sample}/merge.rmdup.${sample}.cram && rm ../${sample}/merge.rmdup.${sample}.bam
    samtools index ../${sample}/merge.rmdup.${sample}.cram 
    " >> $sample.go.sh

    echo '
    ##nos##
    cat > '$sample.qsubcall.exp' << '\'EOF\''
    #!/usr/bin/expect  -f
    set scriptname [lindex $argv 0] 
    set timeout 10
    spawn ssh login4
    expect "$"
    send "cd '`pwd`' ;QS -m 40 -p 10 -i 1022 -s $scriptname;exit;\r"
    #interact
    expect eof
    EOF
    
    chmod +x '$sample.qsubcall.exp'
    
    for j in `cat chr.list`;do echo "       /public/agis/chengshifeng_group/fengcong/WGRS/software/gatk-4.1.2.0/gatk HaplotypeCaller   -R $wrfp   -I ../${sample}/merge.rmdup.${sample}.cram         --native-pair-hmm-threads 10         -L $j         -O ../${sample}/${sample}.HaplotypeCaller.$j.g.vcf.gz    -ERC GVCF  " > ${sample}.$j.call.sh && /usr/bin/expect ./'$sample.qsubcall.exp' ${sample}.$j.call.sh;    done
    ##nos##
    ' >> $sample.go.sh && QS -m 40 -p10 -i 1023 -s $sample.go.sh
    
done
