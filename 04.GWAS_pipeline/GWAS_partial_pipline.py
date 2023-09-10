#--coding:utf-8--
#GWAS partial ld heat map and manhattanplot

import sys
import gzip
import argparse
import os
from datetime import datetime
import subprocess
import math


script_path = os.path.abspath(__file__)
manplot=os.path.dirname(script_path)+"/"+"manhattanplot_function.r"
qqplot=os.path.dirname(script_path)+"/"+"pretty_QQ_plot_func.r"
mergepy=os.path.dirname(script_path)+"/"+"merge_jpg_dirversion.py"
mergepartpy=os.path.dirname(script_path)+"/"+"merge_jpg_part_dirversion.py"
toheat=os.path.dirname(script_path)+"/"+"vcf2fasta_noref.py"
fakeimpute=os.path.dirname(script_path)+"/"+"fake_impute_for_ldheatmap.py"
pypath="/public/home/fengcong/anaconda2/envs/py3/bin/python"
rpath="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"

def item3(item):
    return -math.log(eval(item[2]),10)

if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="GWAS partial pipline(ldheatmap and manhattan plot)")
    cmdparser.add_argument("-r", "--res", dest="res", type=str, required=True,
                           help="format:chr pos pvalue ,abs path")
    cmdparser.add_argument("-c","--chr", dest="chr",type=str, required=True,
                             help="which chr do you wann plot")
    cmdparser.add_argument("-b","--bed", dest="bed",type=str, 
                            help="assign which region to plot .format: start:end ;default:peak left/right 1Mb")
    cmdparser.add_argument("-t","--thres", dest="threshold",type=str, 
                            help="default:0.1/snps")
    cmdparser.add_argument("-o","--outdir", dest="outdir",type=str, required=True,
                            help="output dir.")
    cmdparser.add_argument("-s","--sample", dest="sample",type=str, required=True,
                            help="sample list ,abs path")
    cmdparser.add_argument("-v","--vcf", dest="vcf",type=str, required=True,
                            help="vcffile ,abs path ")
    cmdparser.add_argument("-p","--pn", dest="phenon",type=str, 
                            help="pheno name ,used for plot ")
    
    args = cmdparser.parse_args()

    #args.bed
    print(args)

    resf=open(args.res,"r")
    lcount=-1
    line=resf.readline()
    while line:
        lcount+=1
        line =resf.readline()

    thres=-math.log(0.1/lcount,10)
    if args.threshold != None:
        thres=args.threshold
    
    #find biger than threshold position
    os.system("mkdir -p %s"%(args.outdir))
    ouf=open(args.outdir+"/%s_bigerthan_%s.tsv"%(args.chr,thres),"w")
    ouf.write("chr\tpos\tpvalue\n")

    big=[]

    resf.seek(0) 
    line = resf.readline()

    line=resf.readline()
    while line:
        ls =line.strip().split()
        if ls[0] == args.chr:
            if -math.log(eval(ls[2]),10) > thres:
                ouf.write(line)
                big.append(ls)
        line=resf.readline()
    ouf.close()
    
    #judge plot pos
    if len(big) == 0:
        print("this chr no peak ")
        exit(-1)
    
    big.sort(key=item3,reverse=True)

    pos = int(big[0][1]) 

    left=pos-1000000
    right=pos + 1000000
    if args.bed != None:
        left=int(args.bed.split(":")[0])
        right=int(args.bed.split(":")[1])
    
    if left <0:
        left=0
    
    #get plot manhattan data
    resf.seek(0)
    ouf=open(args.outdir+"/%s_%d-%d_plotman.tsv"%(args.chr,left,right),"w")
    ouf1=open(args.outdir+"/%s_%d-%d_heatpos.tsv"%(args.chr,left,right),"w")
    ouf2=open(args.outdir+"/%s_%d-%d_heatid.tsv"%(args.chr,left,right),"w")
    ouf.write("chr\tpos\tpvalue\n")
    line = resf.readline()

    line=resf.readline()
    while line:
        ls =line.strip().split()
        if ls[0] == args.chr:
            if int(ls[1]) >= left and int(ls[1]) <= right:
                ouf.write(line)
                ouf1.write(ls[1]+"\n")
                ouf2.write("%s_%s\n"%(ls[0],ls[1]))
        line=resf.readline()

    ouf.close()
    ouf1.close()
    ouf2.close()

    #get plot heatmap data
    os.chdir(args.outdir)
    cmd ='''
    ln -s %s .
    ln -s %s .
    '''%(args.vcf,args.sample)

    

    outsh=open("00.ln_data.sh","w")
    outsh.write(cmd)
    outsh.close()
    os.system(cmd)

    cmd='''
    vcftools --gzvcf ./%s --snps %s --keep %s --recode --recode-INFO-all  --stdout >%s
    ''' %(os.path.basename(args.vcf),os.path.basename(args.outdir+"/%s_%d-%d_heatid.tsv"%(args.chr,left,right)),os.path.basename(args.sample),\
        "%s_%d-%d_heat.vcf"%(args.chr,left,right))

    outsh=open("01.get_snp_data.sh","w")
    outsh.write(cmd)
    outsh.close()

    #exec 
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = p.communicate()  #阻塞父进程
    for line in out.splitlines():
        print(line)
    #septoria_snp.vcf.geno.txt
    geno="%s_%d-%d_heat.vcf"%(args.chr,left,right) + ".geno.txt"
    cmd='''
    %s %s %s %s
    %s %s %s %s
    ''' %(pypath,toheat,"%s_%d-%d_heat.vcf"%(args.chr,left,right),os.path.basename(args.sample),\
        pypath ,fakeimpute, geno, geno+".fakeimpute.txt")

    outsh=open("02.toheat_fakeimpute.sh","w")
    outsh.write(cmd)
    outsh.close()

    #exec 
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = p.communicate()  #阻塞父进程
    for line in out.splitlines():
        print(line)
    
    plotr='''
    args=commandArgs(T)
    snpdata=args[1]
    snppos=args[2]
    filen=args[3]
    ## 导入R包
    library("LDheatmap")
    library("genetics")
    library("grid")
    ## 导入SNP标记基因型信息
    SNPdata <- read.table(snpdata,header=T,sep="\t")
    ## 导入SNP标记位置信息
    SNPpos <- read.table(snppos,header=F,sep="\t")
    ## 将SNP基因型信息转换为genotype格式
    num <- ncol(SNPdata)
    for(i in 1:num){
    SNPdata[,i]<-as.genotype(SNPdata[,i])
    }
    ## 将SNP位置信息转换为vector
    pos <- as.vector(unlist(SNPpos))
    ## 设置热图颜色
    color.rgb <- colorRampPalette(rev(c("white","red")),space="rgb")
    ## 绘制连锁不平衡图

    # pdf("ldheatmap.pdf",width = 5,height = 5)
    #"lgtAUDPCpycnidia_chr3A_26Mto28M_ldheatmap.jpg"
    jpeg(filen,width = 15,height=15,units = 'cm', res = 458)
    LDheatmap(SNPdata,pos,color=color.rgb(20),flip=TRUE)
    grid.edit(gPath("ldheatmap","heatMap","heatmap"),gp=gpar(col="white",lwd=0.1))
    dev.off()
    '''
    outsh=open("ldheatmap.r","w")
    outsh.write(plotr)
    outsh.close()

    plotr='''
    args=commandArgs(T)
    filen=args[1]
    titlenn=args[2]
    source("%s")
    #source("./GAPIT.QQ.FUNC_custom.r")
    source("%s")
    manhattanplot(mydata = filen, key = titlenn, chr = "all",title=titlenn, log10 = T, vline=F,x_tick_labs = "Character",thresholds=%f)
    mydata=read.table(filen,head=T)
    mydata=na.omit(mydata)
    #jpeg(paste( titlenn,"_QQ-Plot.jpg" ,sep = ""),width = 15,height=15,units = 'cm', res = 225)
    #qqunif.plot(mydata[,3])
    #dev.off()
    #GAPIT.QQ(P.values = mydata[,3], name.of.trait = titlenn,DPP=50000)
    '''%(manplot,qqplot,thres)

    outsh=open("manhattan.r","w")
    outsh.write(plotr)
    outsh.close()

    tn="triat"
    if args.phenon != None:
        tn=args.phenon

    cmd ='''
    %s ldheatmap.r  %s %s %s 
    %s manhattan.r  %s %s
    '''%(rpath,geno+".fakeimpute.txt","%s_%d-%d_heatpos.tsv"%(args.chr,left,right),tn+"_%s_%dto%d_ldheatmap.jpg"%(args.chr,left,right),\
        rpath,"%s_%d-%d_plotman.tsv"%(args.chr,left,right),tn+"_%s_%dto%d"%(args.chr,left,right))

    outsh=open("03.go_plot.sh","w")
    outsh.write(cmd)
    outsh.close()

    #exec 
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = p.communicate()  #阻塞父进程
    for line in out.splitlines():
        print(line)


    cmd ='''
    %s %s .
    ''' %(pypath,mergepartpy)
    outsh=open("04.merge.sh","w")
    outsh.write(cmd)
    outsh.close()

    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = p.communicate()  #阻塞父进程
    for line in out.splitlines():
        print(line)

    print("done")


