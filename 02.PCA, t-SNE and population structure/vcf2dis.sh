
            #!/usr/bin/sh
            #-*# It has been modified by QS.py : vcf2dis.sh.TopHeZBE #-*#
            cat /public/home/fengcong/.chenglab
            #change dir
            cd /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis ;
            #time
            echo Running on host `hostname`
            echo PID $$
            time1=`date +"%Y-%m-%d %H:%M:%S"`
            echo Start Time is `date`
            echo Directory is `pwd`
            #stats,how to del: sed -i "/##--\/--##$/d" your_work.sh
            stats="s"
            stats_num=1
            echo "-------------------------------------------------------------------"
            #origin script content
            #########################################################################start
            
            #!/bin/sh
            dir=/vol3/agis/chengshifeng_group/chenghong/wheat/01.SNPdatasets/00.input/watkins/1MB_window
            cd /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis
            cat /vol3/agis/chengshifeng_group/chenghong/wheat/04.diversity/01.pi/03.group_1047/AG1.chr1A.1M.windowed.pi |while read line
            do
            chr=`echo $line | awk '{print $1}'`
            start=`echo $line | awk '{print $2}'`
            end=`echo $line | awk '{print $3}'`
            /vol3/agis/chengshifeng_group/chenghong/software/VCF2Dis/bin/VCF2Dis -InPut $dir/${chr}/${chr}-${start}-${end}.vcf.gz -OutPut ${chr}/${chr}-${start}-${end}.p_dis.mat
            done
            
            


            #########################################################################end
            echo "-------------------------------------------------------------------"
            echo End Time is `date`
            time2=`date +"%Y-%m-%d %H:%M:%S"`
            timerun1=$(($(date +%s -d "$time2") - $(date +%s -d "$time1")))
            echo $timerun1 | awk '{print "Running time is " $1*1/3600 " hours"}'
            echo $stats
            #qsub:qsub -V -p 0 -q csf_queue -l mem=10G,nodes=1:ppn=2 -M LE0AWrJb /vol3/agis/chengshifeng_group/chenghong/wheat/02.population/08.p_dis/vcf2dis.sh
	