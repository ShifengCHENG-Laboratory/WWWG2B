#!/public/home/fengcong/anaconda2/envs/py3/bin/python
#2021-3-11 10:37:21
#usage:python __file__ -v in.vcf/in.vcf.gz/in.vcf.list.file/in.vcf.gz.list.file -o out_dir -t recombinant_to_watseq.csv  *_dfw_genotypes_with_RefSeq1_pos.csv

#version 3 will reduce  the consumption of thread locks
#version 3_test  speed up 
#version 4 consider haplo-block

# -v  --vcf , u can specify a vcf/gzvcf file, or a vcf/gz list file that one vcf/gz file per line
# -o --ouputdir 
# -t --trans , The file has two columns, separated by commas
    # eg:
    # ParW007,WATDE0002
    # ParW032,WATDE0004
    # ParW034,WATDE0005

# *_dfw_genotypes_with_RefSeq1_pos.csv  eg:
    # id,WGA.Chr,WGA.bp,chrom,cM,ParW299_0001,ParW299_0002,ParW299_0003,ParW299_0004,ParW299_0005
    # AX-95154988,1A,14253941,1A,0,H,A,A,B,A
    # AX-94599217,1A,21624369,1A,7.3,H,A,A,B,A
    # AX-95630408,1A,24247344,1A,13,H,A,A,B,A
    # AX-94415384,1A,305004141,1A,26,H,A,A,B,A
    # AX-94481256,1A,117883285,1A,26,H,A,A,B,A
    # AX-95630392,1A,101670410,1A,26,H,A,A,B,A
    # AX-94936813,1A,508319045,1A,38.1,H,A,B,B,A
    # AX-94676999,1A,534265526,1A,53.2,H,B,B,A,B
    # BS00022239,1A,535027890,1A,54.5,H,B,B,A,B

# genotype: A : paragon ,B : Watkins , X or H : heterozygote , - : missing
# note: WGA.Chr or WGA.bp can be "NA"
# note: Imputate is not executed when paragon/watkins is missing or het
# note: if dont have axiom marker on the chromosome , I will not impute this chromosome
block_file = "/vol3/agis/chengshifeng_group/huangzejian/06.wheat_pop_analysis/02.LDblock_identify/01.plot_data/03.functional/09.all_chromosome_blocks_excel_haplotype_matrix/sorted.all.blocks.haplotype.matrix.txt"


import argparse
import sys,os
import numpy as np
import gzip  #read gzip or bgzip file
from bgzip import BGZipWriter # write bgzip file
from multiprocessing import Process
from multiprocessing.managers import BaseManager
from multiprocessing import Manager,Lock
from queue import Queue
from bisect import bisect

import json

import time
bgzip = "/public/home/fengcong/anaconda2/bin/bgzip"
bcftools = "/public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools"
class Manager2(BaseManager):
    pass

Manager2.register('get_priorityQueue', Queue)

class TXTBGZipWriter(BGZipWriter):
    def write(self, text):
        super().write(text.encode("utf-8"))

def find_nearest(sorted_array, value):
    '''
    Find the closest value to args_value in the array
    array: sorted numpy arrary
    value: int 
    return: (nearest_value,index)
    '''
    idx = np.searchsorted(sorted_array, value, side="left")
    if idx >= len(sorted_array):
        nearest = sorted_array[-1]
        idx = len(sorted_array)-1
    elif idx == 0:
        nearest = sorted_array[0]
    else:
        if abs(value - sorted_array[idx-1]) < abs(value - sorted_array[idx]):
            nearest = sorted_array[idx-1]
            idx -= 1
        else:
            nearest = sorted_array[idx]
    return (nearest,idx)

def read_axiom_data(RILs_list,id_trans_file):
    '''
    RILs_list: input axiom data list
    id_trans_file: recombinant to Watkins parent
    
    '''
    #axiom_data_dict: I'm going to put the data in here
    axiom_data_dict = {}

    #read the id trans table
    id_trans_dict = {}
    inf = open(id_trans_file,"r")
    for line in inf.readlines():
        ls = line.strip().split(",")
        id_trans_dict[ls[0]] = ls[1]
    inf.close()
    
    #deal axiom data
    file_reader_list = [open(file,"r") for file in RILs_list]

    #read the axiom data
    for file_reader in file_reader_list:
        #deal header
        header = file_reader.readline().strip().split(",")
        recombinant_list = header[5:] 

        line = file_reader.readline()
        while line:
            ls = line.strip().split(",")
            chromosome = "chr"+ls[1]
            position = ls[2]
            for index,geno in enumerate(ls[5:]):
                if recombinant_list[index].split("_")[0] not in id_trans_dict:
                    continue
                watkins_parent = id_trans_dict[recombinant_list[index].split("_")[0]]
                #go go go
                if ls[1] != "NA" and ls[2] != "NA":
                    pass
                else:
                    chromosome = "chrUnknown"
                    position = 0
                axiom_data_dict.setdefault(watkins_parent,{}) \
                .setdefault(recombinant_list[index],{}) \
                .setdefault(chromosome,{}) \
                .setdefault("id",[]).append( ls[0] )

                axiom_data_dict.setdefault(watkins_parent,{}) \
                .setdefault(recombinant_list[index],{}) \
                .setdefault(chromosome,{}) \
                .setdefault("position",[]).append( int(position) )

                axiom_data_dict.setdefault(watkins_parent,{}) \
                .setdefault(recombinant_list[index],{}) \
                .setdefault(chromosome,{}) \
                .setdefault("genotype",[]).append( geno )

            line = file_reader.readline()

    #sort the data
    # idx_sorted = np.argsort(array)
    # sorted_array = np.array(array[idx_sorted])

    for watkins_parent in axiom_data_dict:
        for recombinant in axiom_data_dict[watkins_parent]:
            for chromosome in axiom_data_dict[watkins_parent][recombinant]:
                tmp_d = axiom_data_dict[watkins_parent][recombinant][chromosome]
                idx_sorted = np.argsort(tmp_d["position"])
                tmp_d["position"] = np.array(np.array(tmp_d["position"])[idx_sorted])
                tmp_d["genotype"] = np.array(np.array(tmp_d["genotype"])[idx_sorted])
                tmp_d["id"] = np.array(np.array(tmp_d["id"])[idx_sorted])


    #close all the axiom file
    for file_reader in file_reader_list:
        file_reader.close()

    return axiom_data_dict


def write_vcf_header(outputfilepath,vcf_sample_list):
    ouf = open(outputfilepath,"wb")
    w1 = TXTBGZipWriter(ouf)
    w1.write('''##fileformat=VCFv4.2
##INFO=<ID=axiomMarker,Number=1,Type=String,Description="imputation marker, reference Paragon">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1A,length=594102056>
##contig=<ID=chr1B,length=689851870>
##contig=<ID=chr1D,length=495453186>
##contig=<ID=chr2A,length=780798557>
##contig=<ID=chr2B,length=801256715>
##contig=<ID=chr2D,length=651852609>
##contig=<ID=chr3A,length=750843639>
##contig=<ID=chr3B,length=830829764>
##contig=<ID=chr3D,length=615552423>
##contig=<ID=chr4A,length=744588157>
##contig=<ID=chr4B,length=673617499>
##contig=<ID=chr4D,length=509857067>
##contig=<ID=chr5A,length=709773743>
##contig=<ID=chr5B,length=713149757>
##contig=<ID=chr5D,length=566080677>
##contig=<ID=chr6A,length=618079260>
##contig=<ID=chr6B,length=720988478>
##contig=<ID=chr6D,length=473592718>
##contig=<ID=chr7A,length=736706236>
##contig=<ID=chr7B,length=750620385>
##contig=<ID=chr7D,length=638686055>
##contig=<ID=chrUn,length=480980714>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT''')
    for sample in vcf_sample_list:
        w1.write("\t%s"%(sample))
    w1.write("\n")
    sys.stderr.write("write %s\n"%(outputfilepath))
    sys.stderr.flush()
    # w1.close()
    return w1

def write_vcf_header2(outputfilepath,vcf_sample_list):
    ouf = open(outputfilepath,"w")
    # w1 = TXTBGZipWriter(ouf)
    w1 = ouf
    w1.write('''##fileformat=VCFv4.2
##INFO=<ID=axiomMarker,Number=1,Type=String,Description="imputation marker, reference Paragon">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1A,length=594102056>
##contig=<ID=chr1B,length=689851870>
##contig=<ID=chr1D,length=495453186>
##contig=<ID=chr2A,length=780798557>
##contig=<ID=chr2B,length=801256715>
##contig=<ID=chr2D,length=651852609>
##contig=<ID=chr3A,length=750843639>
##contig=<ID=chr3B,length=830829764>
##contig=<ID=chr3D,length=615552423>
##contig=<ID=chr4A,length=744588157>
##contig=<ID=chr4B,length=673617499>
##contig=<ID=chr4D,length=509857067>
##contig=<ID=chr5A,length=709773743>
##contig=<ID=chr5B,length=713149757>
##contig=<ID=chr5D,length=566080677>
##contig=<ID=chr6A,length=618079260>
##contig=<ID=chr6B,length=720988478>
##contig=<ID=chr6D,length=473592718>
##contig=<ID=chr7A,length=736706236>
##contig=<ID=chr7B,length=750620385>
##contig=<ID=chr7D,length=638686055>
##contig=<ID=chrUn,length=480980714>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT''')
    for sample in vcf_sample_list:
        w1.write("\t%s"%(sample))
    w1.write("\n")
    sys.stderr.write("write %s\n"%(outputfilepath))
    sys.stderr.flush()
    # w1.close()
    return w1

def thread_func(job_queue,axiom_data_dict,block_boundaries_list_dict,output_file_dict,lock,thread_id,threads_id_list,work_status):
    #job queue
    #(priority,flag,writer_object,paragon_gt,watkins_index,watkins_name,ls)
    sys.stderr.write("thread %d start-jobs:%d\n"%(thread_id,len(output_file_dict)))
    sys.stderr.flush()
    output_writer_dict = {}
    for indx in output_file_dict:
        output_writer_dict[indx] =  write_vcf_header2(output_file_dict[indx][0], list(axiom_data_dict[output_file_dict[indx][1]].keys()) )

    while True:
        if work_status[thread_id] != 1:
            continue

        if len(job_queue) == 0 :
            lock.acquire()
            work_status[thread_id] = 0
            work_status["sum"] = dict_sum(work_status,threads_id_list)
            lock.release()
            continue

        
        para_tuple = job_queue.pop(0)

        
        flag = para_tuple[1]

        if flag == 0:
            for i in output_writer_dict:
                output_writer_dict[i].close()
                os.system("%s %s"%(bgzip,output_file_dict[i][0]))
                os.system("%s index %s.gz"%(bcftools,output_file_dict[i][0]))
                sys.stderr.write("thread %d end:%s\n"%(thread_id,output_file_dict[i][0]))
                sys.stderr.flush()
            break
            
        writer_object = output_writer_dict[para_tuple[2]]
        paragon_gt = para_tuple[3]
        watkins_index = para_tuple[4]
        watkins_name = para_tuple[5]
        ls = para_tuple[6]
        

        chromosome = ls[0]
        position = int(ls[1])
        

        #if no axiom marker on this chromosome , than next
        rnd_recombinant = list(axiom_data_dict[watkins_name].keys())[0]
        if chromosome not in axiom_data_dict[watkins_name][rnd_recombinant].keys():
            continue

        w_gt = ls[watkins_index].split(":")[0]
        if w_gt[0] == ".":
            continue
        elif w_gt[0] == w_gt[-1]:
            if paragon_gt[0] == "0":
                pass
            else:
                if w_gt[0] == "0":
                    w_gt = "1/1"
                else:
                    w_gt = "0/0"
        else:#wgt[0] != wgt[-1]
            continue


        gt_list = []
        info = "axiomMarker="
        paragon_gt_impute = "0/0"

        nearnest_pos,idx = find_nearest(axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"],position)

        ##判断这个snp在哪个block
        block_index = bisect(block_boundaries_list_dict[chromosome],position)
        if block_index % 2 ==1:
            block_start = block_boundaries_list_dict[chromosome][block_index-1]
            block_end = block_boundaries_list_dict[chromosome][block_index]

            idx_start = np.searchsorted(axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"], block_start, side="left")
            idx_end = np.searchsorted(axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"], block_end, side="right")
            temp_pos = axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"][idx_start:idx_end]
            ##判断这个block里面是否有marker分型数据
            if len(temp_pos) != 0:
                nearnest_pos2,idx2 = find_nearest(temp_pos,position)
                idx = np.where(axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"]==nearnest_pos2)[0][0]
            

        # #if in a RIL , the marker is same ,u can just search once
        # nearnest_pos,idx = find_nearest(axiom_data_dict[watkins_name][rnd_recombinant][chromosome]["position"],position)
        for recombinant in axiom_data_dict[watkins_name]:
            # print(recombinant,idx,position,nearnest_pos)
            nearnest_id = axiom_data_dict[watkins_name][recombinant][chromosome]["id"][idx]
            nearnest_geno = axiom_data_dict[watkins_name][recombinant][chromosome]["genotype"][idx]
            info =  "axiomMarker=" + nearnest_id
            if nearnest_geno == "A":
                gt_list.append(paragon_gt_impute)
            elif nearnest_geno == "B":
                gt_list.append(w_gt)
            elif nearnest_geno == "X" or nearnest_geno== "H":
                if paragon_gt_impute == w_gt:
                    gt_list.append("0/0")
                else:
                    gt_list.append("0/1")
                
            else:
                gt_list.append("./.")
        
        write_list = []
        if paragon_gt[0] == "1":
            write_list = [ls[0],ls[1],ls[2],ls[4],ls[3],".",".",info,"GT"]+gt_list
        else:#"0"
            write_list = [ls[0],ls[1],ls[2],ls[3],ls[4],".",".",info,"GT"]+gt_list
        
        # sys.stdout.write("%s\n"%(str(write_list)))
        # sys.stdout.flush()
        # print(str(writer_object))
        writer_object.write("\t".join(write_list)+"\n")

def chunks(arr, m):
    l=[]
    for i in range(m):
        l.append([])
    i = 0
    j = 0
    while True:
        if j >= len(arr):
            break
        l[i].append(arr[j])
        i = (i+1)%m
        j+=1
    return l

def dict_sum(d,keys):
    sum = 0
    for key in keys:
        sum+=d[key]
    return sum

def impute(axiom_data_dict, vcf_file_or_list,output_prefix,threads,block_boundaries_list_dict):
    #deal the vcf file list
    if vcf_file_or_list:
        file_list = []
        if vcf_file_or_list.endswith(".vcf") or vcf_file_or_list.endswith(".vcf.gz"):
            file_list.append(vcf_file_or_list)
        else:
            inf = open(vcf_file_or_list,"r")
            for line in inf.readlines():
                file_list.append(line.strip())
            inf.close()
        
        #open the vcf files
        file_reader_list = [gzip.open(file,"rt") if file.endswith(".gz") else open(file,"r") for file in file_list]
    else:
        file_reader_list = [sys.stdin]

    #which watseq sample will be use
    watseq_index = []
    output_writer_list = []
    paragon_index = None

    #creat a proxy
    m = Manager2()
    m.start()
    #
    # axiom_data_dict_m = Manager().dict(axiom_data_dict)

    used_threads_num = 1

    jobs_queue_d = {}

    threads_list=[]
    threads_id_list = []
    queues_list=[]
    output_file_dict={}
    work_status = Manager().dict({})
    lock = Lock()
    
    

    for file_reader in file_reader_list:
        #skip vcf header
        line = file_reader.readline()
        while line.startswith("##"):
            line = file_reader.readline()
        ouf_list_file = open("%s.list"%(output_prefix),"w")
        header_list = line.strip().split()
        #open output file and write header
        if len(watseq_index) == 0:
            for watkins_parents in axiom_data_dict.keys():
                if watkins_parents in header_list:
                    watseq_index.append(header_list.index(watkins_parents))
                    
                    output_file_dict[header_list.index(watkins_parents)] = \
                        (output_prefix+".Paragonx%s.vcf"%(watkins_parents),watkins_parents)
                    
                    ouf_list_file.write("%s\n"%(output_prefix+".Paragonx%s.vcf.gz"%(watkins_parents)))
                    
                else:
                    sys.stderr.write("%s not in WatSeq project!\n"%(watkins_parents))
                    sys.stderr.flush()
            ouf_list_file.close()

            paragon_index = header_list.index("paragon")
            #calc how many threads should be used
            if len(watseq_index) >= threads:
                used_threads_num = threads
            else:
                used_threads_num = len(watseq_index)

            #Assigning Tasks
            watseq_index_chunk_list = chunks(watseq_index,used_threads_num)
            thread_id = 0
            for chunk in watseq_index_chunk_list:
                tmp_output_file_dict={}
                # tmppq = m.get_priorityQueue()
                tmppq = Manager().list([])
                queues_list.append(tmppq)
                for w_index in chunk:
                    jobs_queue_d[w_index] = tmppq
                    tmp_output_file_dict[w_index] = output_file_dict[w_index]
                threads_list.append(Process(target=thread_func,args=(tmppq,axiom_data_dict,block_boundaries_list_dict,tmp_output_file_dict,lock,thread_id,[id for id in range(0,used_threads_num)],work_status)))
                work_status[thread_id] = 0
                threads_id_list.append(thread_id)
                thread_id +=1
            #start threads 
            work_status["sum"] = 0
            for thread in threads_list:
                thread.start()

        else:
            pass
        
        
        # jobs_queue_d = Manager().dict(jobs_queue_d)
        #read the vcf body
        line = file_reader.readline()
        k=0
        while line:
            ls = line.strip().split()
            
            

            #paragon status
            paragon_gt = ls[paragon_index].split(":")[0]
            if paragon_gt[0] == paragon_gt[-1] and paragon_gt[0] != ".":
                pass
            else: #paragon is heterozygote or missing
                line = file_reader.readline()
                continue

            #traverse watkins parents
            for file_writer_index,w_index in enumerate( watseq_index):
                # print((0,1,output_writer_list[file_writer_index],\
                #     paragon_gt,w_index,header_list[w_index],ls))
                jobs_queue_d[w_index].append((0,1,w_index,paragon_gt,w_index,header_list[w_index],ls))
                


            

            k+=1
            if k % 500 == 0 :
                sys.stderr.write("deal %d SNPs...\n"%(k))
                sys.stderr.flush()
                lock.acquire()
                for id in threads_id_list:
                    work_status[id] = 1
                work_status["sum"] = dict_sum(work_status,threads_id_list)
                lock.release()

                while True:
                    time.sleep(0.1)
                    lock.acquire()
                    tmpv =  work_status["sum"] 
                    lock.release()
                    if tmpv != 0:
                        continue
                    else:
                        break
            
            line = file_reader.readline()

        #send job over signal
        while True:
            time.sleep(0.5)
            lock.acquire()
            tmpv =  work_status["sum"] 
            lock.release()
            if tmpv != 0:
                continue
            else:
                break
        
        lock.acquire()
        for w_index in jobs_queue_d:
            jobs_queue_d[w_index].append((0,0))
        for id in threads_id_list:
            work_status[id] = 1
        work_status["sum"] = dict_sum(work_status,threads_id_list)
        lock.release()
   
    #wait jobs done
    for thread in threads_list:
        thread.join()
    
    #close the ouput file
    # for file_writer in output_writer_list:
    #     file_writer.close()

    #close the vcf files
    for file_reader in file_reader_list:
        file_reader.close()

    #index outpufile
    
    for w_index in output_file_dict.keys():
        os.system("bcftools index %s"%(output_file_dict[w_index][0]))

def read_block_data(block_file,axiom_data_dict):
    ###data 
    block_boundaries_list_dict = {}
    # block_status_dict = {}
    ###data
    with open(block_file,"r") as block_f:
        ##header
        header_line = block_f.readline()
        line = block_f.readline()
        while line:
            ls = line.strip().split()
            block_ls = ls[0].split("_")

            block_id = ls[0]
            block_chr = block_ls[0]
            block_start = int(block_ls[1])
            block_end = int(block_ls[2])

            block_boundaries_list_dict.setdefault(block_chr,[]).append(block_start)
            block_boundaries_list_dict.setdefault(block_chr,[]).append(block_end)
            ##construct block status dict
            # block_status_dict[block_id] = {}


            ##construct block status dict
            
            line = block_f.readline()

    return block_boundaries_list_dict

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="NAM imputation/fengcong@caas.cn" )
    cmdparser.add_argument("-v","--vcf", dest="vcf",type=str, required=False,
                            help="a vcf/gz file or a vcf/gz list file that one file per line")
    cmdparser.add_argument("-t","--trans", dest="trans",type=str, required=True,
                            help="recombinant id to WatSeq ID,The file has two columns, separated by commas")
    cmdparser.add_argument("-o","--outputprefix", dest="outputprefix",type=str, required=True,
                            help="output prefix")   
    cmdparser.add_argument("-b","--blockfile", dest="blockfile",type=str, required=False,default=block_file,
                            help="block file")   
    cmdparser.add_argument("-@","--threads", dest="threads",type=int, default=1,
                            help="Maximum number of threads used during impute")                  
    cmdparser.add_argument('RILs', metavar='RILs', type=str, nargs='+',
                    help='each file represents a RIL population and contains the genotype of each recombinant')

    args = cmdparser.parse_args()

    ## deal args
    output_prefix = args.outputprefix
    id_trans_file = args.trans
    vcf_file_or_list = args.vcf
    RILs_list = args.RILs
    threads = args.threads

    #Design how to store data
    ## I will only traverse vcf once
    ## The entrance is the Watkins parent. 
    ## I can quickly access the corresponding recombinant and all the genotypes of this recombinant. 
    ## Do not consider the memory consumption of axiom data, which will be very small.
    # Watkins_parent
    #               |
    #                ----ParW299_0001
    #               |                |
    #               |                 ----chr1A
    #               |                |         |
    #               |                |          ----id: [ID,ID,...]
    #               |                |          ----pos:[position,position,...]   #Redundant storage
    #               |                |          ----geno:[genotype,genotype,...]  #Redundant storage
    #               |                |
    #               |                 ----chr1B
    #               |                | 
    #               |                 ----chr1D
    #               |                |
    #               |                 ----...
    #               |
    #                ----ParW299_0002
    #               |
    #                ----...

    axiom_data_dict = read_axiom_data(RILs_list,id_trans_file)
    
    block_boundaries_list_dict = read_block_data(block_file,axiom_data_dict)

    with open("block_boundaries_list_dict2","w+") as jsf:
        jsf.write(json.dumps(block_boundaries_list_dict,indent=4))

    impute(axiom_data_dict, vcf_file_or_list,output_prefix,threads,block_boundaries_list_dict)
