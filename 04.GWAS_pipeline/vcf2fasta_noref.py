#--coding:utf-8--
#2020-4-20 10:29:16
import sys
#usage:python __file__ in.vcf sample_out_order.txt
#output file:
#       in.vcf.fasta ×
#       snp_id.txt  ×
#       aln.tab     √

if __name__ == "__main__":
    inf=open(sys.argv[1],"r")
    ordf=open(sys.argv[2],"r")
    #skip header
    line = inf.readline()
    while line.startswith("##"):
        line=inf.readline()
    #header line
    ls=line.split()
    sample_ls=[]
    for i in ls[9:]:
        sample_ls.append(i)
    #order
    sample_ord_ls=[]
    for order in ordf.readlines():
        order = order.strip()
        sample_ord_ls.append(sample_ls.index(order))

    #content
    snpid_ls=[]
    fa_ls = [[] for n in range(len(sample_ls))]
    ref=[]
    line=inf.readline()
    while line:
        ls=line.split()
        snpid_ls.append(ls[2]) #ID
        ref.append(ls[3])
        for index,item in enumerate(ls[9:]) :
            gt=item.split(":")[0]
            gt1=gt[0]
            gt2=gt[-1]
            phased=gt[1]
            if gt1=="." and gt2==".":
                fa_ls[index].append(gt1+phased+gt2)
            else:
                fa_ls[index].append(ls[int(gt1)+3]+phased+ls[int(gt2)+3])

        line=inf.readline()
    inf.close()
    outf_tab=open(sys.argv[1]+".geno.txt","w")
    #header
    #outf_tab.write("Sample\t")
    outf_tab.write("\t".join(snpid_ls))
    outf_tab.write("\n")
    #Ref
    #outf_tab.write("Ref\t")
    #outf_tab.write("\t".join(ref))
    #outf_tab.write("\n")
    #each sample
    for index in sample_ord_ls:
        #outf_tab.write(sample_ls[index]+"\t")
        outf_tab.write("\t".join(fa_ls[index]))
        outf_tab.write("\n")

    outf_tab.close()
