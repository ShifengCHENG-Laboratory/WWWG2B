#coding:utf-8
import sys
import os
import SNPANN_upgrade_3
import filter_codingregion_SNP
import argparse

## script usage:
## python3 block_haplotype_tag_snp.py -c haplotype_cluster.txt -r chr1A:1-100 -v test.vcf.gz -o tag_snp.list

## haplotype_cluster.txt:
## S1   1
## S2   1
## S3   2

## the vcf.gz file should be indexed by tabix

## You must have installed bcftools

## the output file is a list of snp, each line is a snp, the format is:
## chr1A_1
## chr1A_2
bcftools = "/public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools"

def percent_sample_have_correct_hap(hap_combination, selected_snp_index_list ,
    hap_sample_list_d_filter  , sample_genotype_d):

    correct_count = 0 
    sample_list = hap_sample_list_d_filter[hap_combination[0]] + hap_sample_list_d_filter[hap_combination[1]]
    total_count = len(sample_list)

    ## 看看selected_snp_index_list 中的SNP，是否能够区分出两个hap
    ## 现在分为几类？

    new_sample_genotype_d = {} 
    new_hap_info_list = []
    for sample in sample_list:
        new_sample_genotype_d[sample] = []
        for index in selected_snp_index_list:
            new_sample_genotype_d[sample].append( sample_genotype_d[sample][index] )
        
        
        new_sample_genotype_d[sample] = tuple(new_sample_genotype_d[sample])
        new_hap_info_list.append( new_sample_genotype_d[sample] )
    
    new_hap_info_list = list(set(new_hap_info_list)) ## 去重
    new_hap_d = {} ## {hap1:[sample1,sample2],hap2:[sample3,sample4]}
    for i, hap in enumerate(new_hap_info_list):
        new_hap_d[hap] = []
        for sample in sample_list:
            if new_sample_genotype_d[sample] == hap:
                new_hap_d[hap].append(sample)
    
    ## 找到老的分类中 与 新的分类中，交集最多的那个分类
    max_intersection_count = 0
    max_intersection_hap1 = None
    for hap in new_hap_d:
        intersection_count = len( set(new_hap_d[hap]) & set(hap_sample_list_d_filter[hap_combination[0]]) )
        if intersection_count > max_intersection_count:
            max_intersection_count = intersection_count
            max_intersection_hap1 = hap
    
    max_intersection_count = 0
    max_intersection_hap2 = None
    for hap in new_hap_d:
        intersection_count = len( set(new_hap_d[hap]) & set(hap_sample_list_d_filter[hap_combination[1]]) )
        if intersection_count > max_intersection_count:
            max_intersection_count = intersection_count
            max_intersection_hap2 = hap
    
    ## 如果这两个hap 是一样的话，说明有点问题
    if max_intersection_hap1 == max_intersection_hap2:
        print("error hap")
        return 0
    
    ## 看看这两个hap，有多少个样本是正确的
    for i , hap in enumerate( (max_intersection_hap1, max_intersection_hap2) ) :
        correct_count += len( set(new_hap_d[hap]) & set(hap_sample_list_d_filter[hap_combination[i]]) )

    return correct_count / total_count


if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="find block tag SNP /fengcong@caas.cn" )
    cmdparser.add_argument("-r","--region", dest="region",type=str, required=True,
                            help="plot region")
    cmdparser.add_argument("-c","--cluster", dest="cluster",type=str, required=True,
                            help="haplotype cluster file")
    cmdparser.add_argument("-v","--vcf", dest="vcf",type=str, required=True,
                            help="vcf file, should be compressed and indexed by tabix")
    cmdparser.add_argument("-o","--output", dest="output",type=str, required=True,
                            help="output file")
    cmdparser.add_argument("-m","--mac", dest="mac",type=int, required=False, default=5,
                            help="minor hap count, default is 5")

    args = cmdparser.parse_args()

    ## read the cluster file
    sample_list_file = args.output + ".sample.list"
    sample_list_file_f = open(sample_list_file, "w")
    sample_hap_d = {} ## {sample:hap}
    hap_sample_list_d = {} #{hap:[sample1,sample2]}
    with open(args.cluster,"r") as f:
        f.readline()
        for line in f:
            line = line.strip()
            sample,hap = line.split()
            sample_hap_d[sample] = hap
            if hap not in hap_sample_list_d:
                hap_sample_list_d[hap] = [sample]
            else:
                hap_sample_list_d[hap].append(sample)
            
            sample_list_file_f.write(sample + "\n")


    sample_list_file_f.close()


    ## get data and convert format
    inf = os.popen('%s view -r %s -S %s %s'%(bcftools , args.region , sample_list_file ,args.vcf))

    sample_genotype_d = {} ## {sample:[0,-1,]}
    snp_list = [] ## [snp1,snp2,]
    ann_list = [] ## [0,1,] 0:non-coding 1:coding

    header_list = []
    for line in inf:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_list = line.strip().split()
            for sample in header_list[9:]:
                sample_genotype_d[sample] = []
            continue

        vcf_line_list = line.strip().split()
        snp_list.append(vcf_line_list[2])

        ## get ann,
        featureid,annotation,c_xx,p_xx = SNPANN_upgrade_3.SNP_ANN_of_Gene_Structure(vcf_line_list[7])
        if annotation  not in filter_codingregion_SNP.functional_or_not_d:
            sys.stderr.write("function annotaion error\n")
            sys.exit(1)
        else:
            ann_list.append(filter_codingregion_SNP.functional_or_not_d[annotation])

        # output heat matrix
        for index, gt in enumerate( vcf_line_list[9:] ) :
            sample_name = header_list[index+9]
            i = gt.split(":")[0]
            if i == "0/1" or i=="0|1" or i == "1|0": #het
                sample_genotype_d[sample_name].append( 0.9 ) ## same with the HAPPE pipeline
            elif i == "0/0" or i == "0|0" :
                sample_genotype_d[sample_name].append( -1 )
            elif i == "1/1" or i == "1|1" :
                sample_genotype_d[sample_name].append( 1 )
            else:
                sample_genotype_d[sample_name].append( 0 )


    inf.close()



    ## filter the hap , only keep the hap with more than 2 samples
    hap_sample_list_d_filter = {}
    for hap in hap_sample_list_d:
        if len(hap_sample_list_d[hap]) >= args.mac:
            hap_sample_list_d_filter[hap] = hap_sample_list_d[hap]

    ## 生成两两hap 之间的组合
    hap_combination_list = [] ## [(hap1,hap2),(hap1,hap3),]
    for hap1 in hap_sample_list_d_filter:
        for hap2 in hap_sample_list_d_filter:
            if hap1 == hap2:
                continue
            hap_combination_list.append((hap1,hap2))

    hap_combination_distance_list_d = {} ## {(hap1,hap2):[distance1,distance2,]}
    hap_combination_distance_list_index_order_d = {} ## {(hap1,hap2):[index1,index2,]} ## the index is the index of the snp in the snp_list, reverse order of distance

    ## 计算两两hap 之间的 每一个SNP的距离
    for hap1 , hap2 in hap_combination_list:
        hap_combination_distance_list_d[(hap1,hap2)] = []
        for index, snp in enumerate(snp_list):
            ## weight : ann_list[index]
            sample_list1 = hap_sample_list_d_filter[hap1]
            sample_list2 = hap_sample_list_d_filter[hap2]
            
            ## sample1_gt_list:  0  0  0  1  -1
            ## sample2_gt_list:  0  0  0  0   1

            stat_of_sample1 = 0 
            stat_of_sample2 = 0

            count_sample1 = 0
            count_sample2 = 0

            for sample in sample_list1:
                stat_of_sample1 += sample_genotype_d[sample][index]
                count_sample1 += 1
            
            for sample in sample_list2:
                stat_of_sample2 += sample_genotype_d[sample][index]
                count_sample2 += 1
            
            stat_of_sample1 = stat_of_sample1 / count_sample1  ## 先组内求平均
            stat_of_sample2 = stat_of_sample2 / count_sample2  ## 先组内求平均

            ## 计算两个组内的平均值的欧式距离
            distance =  (stat_of_sample1 - stat_of_sample2) ** 2 

            if ann_list[index] == 1:
                distance = distance * 4 ## coding region SNP weight 4
            
            hap_combination_distance_list_d[(hap1,hap2)].append(distance)
        
        ## get the order of the distance
        hap_combination_distance_list_index_order_d[(hap1,hap2)] = sorted( range(len(hap_combination_distance_list_d[(hap1,hap2)])), 
            key=lambda k: hap_combination_distance_list_d[(hap1,hap2)][k], reverse=True)

    
    
    ## 看看如果只挑选距离最大的SNP，正确率的分布
    total_selected_snp_index_list = []

    for hap_combination in hap_combination_distance_list_index_order_d:
        selected_snp_index_list = [ hap_combination_distance_list_index_order_d[hap_combination][0] ]
        # correct_percent = percent_sample_have_correct_hap(hap_combination, selected_snp_index_list ,
        #     hap_sample_list_d_filter  , sample_genotype_d)
        # print(hap_combination, correct_percent)
        total_selected_snp_index_list.extend(selected_snp_index_list)

    total_selected_snp_index_list = list(set(total_selected_snp_index_list))
    total_selected_snp_index_list.sort()

    with open(args.output + ".selectedSNP.list" , "w") as outf:
        for index in total_selected_snp_index_list:
            outf.write( snp_list[index] + "\n" )
    
    os.remove(args.output + ".sample.list")