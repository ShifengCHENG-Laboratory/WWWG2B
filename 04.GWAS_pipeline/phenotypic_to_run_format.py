#fengcong@caas.cn
#2020-6-9 11:03:03

import sys
import os
import argparse

#function:
## 1. ID trans , any ID trans to WatSeq ID
## 2. ID correct ,some elite sample's ID is wrong ,but the phenotype data is correct.
## 3. fill  in all sample , and output the same order of vcf.

#input phenotype data format:
#FID IID trait1 trait2
#1190001 1190002 1  2
#WATDE0058 WATDE0058 3 4
#Iena Iena  1   2

#output phenotype data format:
#FID IID trait1 trait2
#WATDE0011 WATDE0011 1 2
#WATDE0058 WATDE0058 3 4
#ability ability  1   2

vcf_order_file="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/x.vol2_wheat_pop_analysis/33.flowering_gene_final/0x.some_required_file/943_vcf_samplename.txt"
correct_file="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/x.vol2_wheat_pop_analysis/33.flowering_gene_final/0x.some_required_file/corrt.txt"
idtrans_file="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/x.vol2_wheat_pop_analysis/33.flowering_gene_final/0x.some_required_file/ID_trans.txt"
if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="deal phenotypic data to GWAS format" )

    cmdparser.add_argument("-p","--phenotype", dest="phenotype",type=str, required=True,
                            help="input phenotype")
    cmdparser.add_argument("-o","--output", dest="output",type=str, required=True,
                            help="output phenotype")

    args = cmdparser.parse_args()

    #vcf order
    vcf_order=open(vcf_order_file,"r")
    vcf=[]
    for line in vcf_order.readlines():
        vcf.append(line.strip())
    vcf_order.close()
    # print(len(vcf),vcf)
    #correct
    correct_file_reader=open(correct_file,"r")
    correct={}
    for line in correct_file_reader.readlines():
        ls=line.strip().split()
        correct[ls[0]]=ls[1]  # key : error id , value : correct id
    correct_file_reader.close()
    # print(len(correct),correct)

    #new id mapping
    idmapping={} #{correct:[errorid1,errorid2]}
    for i in vcf:
        if i in correct:
            tmp=correct[i]
            idmapping.setdefault(tmp,[]).append(i)
        else:
            idmapping.setdefault(i,[]).append(i)
    ##########output the idmapping ##########
    idm=open("idmapping.txt","w")
    for i in idmapping:
        idm.write(i+"\t")
        idm.write("\t".join(idmapping[i])+"\n")
    
    idm.close()

    #########################################

    #id trans
    idtrans_file_reader=open(idtrans_file,"r")
    idtrans={}
    for line in idtrans_file_reader.readlines():
        ls=line.strip().split()
        idtrans.setdefault(ls[1],[]).append(ls[0])
    idtrans_file_reader.close()
    # print(len(idtrans),idtrans)

    #open file
    inf=open(args.phenotype,"r")
    ouf=open(args.output,"w")
    noseq=open(args.phenotype+".noseq","w")

    #
    header=inf.readline()
    ouf.write(header)
    noseq.write(header)

    phenotype=[]
    line=inf.readline()
    while line:
        ls = line.strip().split()
        if ls[0].startswith("119"): #199 ID
            if ls[0] in idtrans:  # if sequencing 
                ttt=ls[0]
                for item in idtrans[ttt]:
                    ls1=ls[:]
                    ls1[0] = item
                    ls1[1] = item
                    phenotype.append(ls1)
            else: # if not sequencing
                noseq.write(line)
        elif ls[0].startswith("WATDE") or ls[0].startswith("B119"): #WatSeq id
            if ls[0] in vcf:  # if sequencing
                phenotype.append(ls)
            else: # if not sequencing
                noseq.write(line)
        else:  # elite id
            if ls[0] in idmapping: #if seq
                ttt=ls[0]
                for k in idmapping[ttt]:  #add each seqid
                    ls1=ls[:]
                    ls1[0],ls1[1]=k,k
                    phenotype.append(ls1)
            else:#not seq
                noseq.write(line)


        line=inf.readline()
    

    has_phenotype_and_seq=[]
    for i in phenotype:
        has_phenotype_and_seq.append(i[0])
    
    # print(has_phenotype_and_seq)
    
    print("has_phenotype_and_seq: %d"%(len(has_phenotype_and_seq)))
    print("uniq:%d"%(len(list(set(has_phenotype_and_seq)))))
    tmp=["NA" for n in header.split()]
    for t in vcf:
        if t not in has_phenotype_and_seq:
            tmp1=tmp[:]
            tmp1[0]=t
            tmp1[1]=t
            phenotype.append(tmp1)
    
    print("after fill: %d"%(len(phenotype)))

    now_order=[]
    for i in phenotype:
        now_order.append(i[0])

    # print(now_order)
    
    #output 
    for i in vcf:
        idx=now_order.index(i)
        ouf.write("\t".join(phenotype[idx])+"\n")

    ouf.close()
            



    
