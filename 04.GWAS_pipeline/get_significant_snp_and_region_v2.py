#python __file__ in.assoc.txt  region gff  outdir

#eg. python __file__ Hd.GWAS.assoc.txt  chr1A:1-200  iwgsc.v1.1.gff   /outdir  
# this version is for natural population
# only support gemma output now

import sys,os

# NAM_gzvcfs = "/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/49.NAM_inputation/04.imputation_output/chr*.merge.vcf.gz"
oneK_gzvcfs = "/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/14.final_SNP_filter2/03.IC_filter/*InbreedingCoeff_retain.vcf.gz"
python3 = "/public/home/fengcong/anaconda2/envs/py3/bin/python"
beta = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/GWAS_make_beta_symbol_right_csv_v2.py"
# beta = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/GWAS_make_beta_symbol_right_csv.py"
csvprepare = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/csv_to_Three-line_table_prepare.py"
csvtoexcel = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/csv_to_Three-line_table.py"
gene_ann = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/gene.func.expression.pl"
plotr_path = "/public/agis/chengshifeng_group/fengcong/WGRS/software/nginx/html/cgi-bin/plot.r"
large_hap_pattern="/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/large_scale_haplotype.py"
sample_inf="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/14.final_SNP_filter2/0x.some_required_file/1059_TEST_group_1_geo.txt"
data_path = "/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/14.final_SNP_filter2/03.IC_filter"
data_suffix="SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz"
mark_list="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/14.final_SNP_filter2/0x.some_required_file/NAM.ID.txt"
NAM_data_path = "/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/49.NAM_inputation/04.imputation_output/"
NAM_data_suffix="merge.vcf.gz"
recombinant_to_watseq="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/49.NAM_inputation/03.axiom_data/recombinant_to_watseq.csv"

def check_abs_path(filename):
    if filename.startswith("/"):
        return filename
    else:
        return os.getcwd()+"/" + filename

if __name__ == "__main__":
    #deal args
    outdir = sys.argv[4]
    inputfile = check_abs_path(sys.argv[1])
    region = sys.argv[2]
    annotation = sys.argv[3]
    if len(sys.argv) == 6:
        mark_list = sys.argv[5]
    

    #
    os.system("mkdir -p %s"%(outdir))
    chromosome = region.split(":")[0]
    chr_start = int(region.split(":")[1].split("-")[0])
    chr_end = int(region.split(":")[1].split("-")[1])

    #cut file
    region_file = outdir+"/"+"%s_%d_%d.txt"%(chromosome,chr_start,chr_end)
    ouf = open(outdir+"/"+"%s_%d_%d.txt"%(chromosome,chr_start,chr_end),"w")
    inf = open(inputfile,"r")
    header = inf.readline()
    ouf.write(header)
    line = inf.readline()
    row_count = 0 
    while line:
        ls = line.strip().split()
        tmp_chr = ls[1].split("_")[0]
        if tmp_chr == chromosome and chr_start <= int(ls[2]) <= chr_end:
            ls[0] = tmp_chr
            ouf.write("\t".join(ls)+"\n")
            row_count +=1
        line = inf.readline()


    inf.close()
    ouf.close()


    ## add beta1 
    ret = os.system("""
    %s %s -i %s   %s > %s
    """%(python3,beta,region_file,oneK_gzvcfs,outdir+"/"+"%s_%d_%d_beta.csv"%(chromosome,chr_start,chr_end)))

    if  ret:
        sys.stdout.write("add beta1 -- retrun code:%d\n"%(ret))
        exit(-1)
    
    #prepare csv 
    ret = os.system("""
    %s %s %s 1 1 "Details of the SNPs in %s:%d-%d"  > %s
    """%(python3,csvprepare,outdir+"/"+"%s_%d_%d_beta.csv"%(chromosome,chr_start,chr_end), \
        chromosome,chr_start,chr_end,outdir+"/"+"%s_%d_%d_beta.csv.pre"%(chromosome,chr_start,chr_end)))
    
    if  ret:
        sys.stdout.write("prepare csv -- retrun code:%d\n"%(ret))
        exit(-1)

    #get gene and annotation
    ret = os.system("""
    awk '{if($3=="gene" && $1 == "%s" && $4 >= %d && $4 <= %d)print}'  %s | sed 's:.*\\tID=\\([^;]*\\);.*:\\1:g' > %s
    
    """%(chromosome,chr_start,chr_end,annotation,outdir+"/"+"%s_%d_%d.gene.list"%(chromosome,chr_start,chr_end)) )

    if  ret:
        sys.stdout.write("get gene list -- retrun code:%d\n"%(ret))
        exit(-1)


        
    #ann
    ret = os.system("""
    perl %s %s  | sed 's:,:;:g' | sed 's:\t:,:g' > %s
    """%(gene_ann, outdir+"/"+"%s_%d_%d.gene.list"%(chromosome,chr_start,chr_end) , \
        outdir+"/"+"%s_%d_%d_gene.csv"%(chromosome,chr_start,chr_end) ) )

    if  ret:
        sys.stdout.write("get gene ann -- retrun code:%d\n"%(ret))
        exit(-1)
    
    #to pre
    ret = os.system("""
    %s %s %s 1 1 "Gene annotation and expression in %s:%d-%d"  > %s
    """%(python3,csvprepare,outdir+"/"+"%s_%d_%d_gene.csv"%(chromosome,chr_start,chr_end), \
        chromosome,chr_start,chr_end,outdir+"/"+"%s_%d_%d_gene.csv.pre"%(chromosome,chr_start,chr_end)) )

    if  ret:
        sys.stdout.write("gene prepare csv -- retrun code:%d\n"%(ret))
        exit(-1)

    #plot local manhattan
    infile = "%s_%d_%d.txt"%(chromosome,chr_start,chr_end)
    ret = os.system("""
    cd %s;
    cat %s  | awk    'BEGIN{OFS="\t"}{print $1,$3,$12}' > %s
    /public/home/fengcong/anaconda2/envs/R/bin/Rscript %s %s  %s.%s_%s_%s  %s  %s %s
    """%(outdir,infile, "%s_%d_%d.addchr.tsv"%(chromosome,chr_start,chr_end),
    plotr_path, "%s_%d_%d.addchr.tsv"%(chromosome,chr_start,chr_end),os.path.basename(outdir), chromosome,chr_start,chr_end, 
    chromosome,chr_start,chr_end
    ) 
    )

    if  ret:
        sys.stdout.write("plot local manhattan -- retrun code:%d\n"%(ret))
        exit(-1)

    os.chdir(outdir)
    # plot large local haplotype pattern : natural population
    ret = os.system(
        """
        %s  %s -v %s/%s.%s \\
        -g %s -r %s:%d-%d -i %s -o %s -t 40 -m %s
        """%(python3, large_hap_pattern, data_path,chromosome , data_suffix,\
            annotation,chromosome,chr_start,chr_end,sample_inf , "%s_%d_%d.1051.large.pattern"%(chromosome,chr_start,chr_end),\
            mark_list
            )
    )

    if  ret:
        sys.stdout.write("hap pattern of natural population -- retrun code:%d\n"%(ret))
        exit(-1)

    # os.chdir("../")
    