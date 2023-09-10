
import sys


#python __file__ in.geno.txt out.geno.txt

if __name__ == "__main__":
    inf=open(sys.argv[1],"r")
    ouf=open(sys.argv[2],"w")

    line = inf.readline()
    ouf.write(line)

    line=inf.readline()
    preline=[]
    while line:
        if len(preline) == 0:
            ls=line.strip().split()
            for i,k in enumerate(ls) :
                if k=="./.":
                    ls[i]="A/A"
            
            preline = ls
            ouf.write('\t'.join(ls) + "\n")
        else:
            ls=line.strip().split()
            for i,k in enumerate(ls) :
                if k=="./.":
                    ls[i]=preline[i]
                
            preline = ls
            ouf.write('\t'.join(ls) + "\n")
        line=inf.readline()
