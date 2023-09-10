library(bigmemory)
library(FarmCPUpp)

phnotype="/gxyy-local/home/liyongyao/fengcong/01.GWAS/23.farmcpu_disease_new/disease_pheno_for_farcpu.txt"
numgd="/gxyy-local/home/liyongyao/fengcong/01.GWAS/22.farmcpu_28_new/GAPIT.Genotype.Numerical.txt"
numif="/gxyy-local/home/liyongyao/fengcong/01.GWAS/22.farmcpu_28_new/GAPIT.Genotype.map.txt"
#pointer="/gxyy-local/home/liyongyao/fengcong/01.GWAS/04.farmcpupp/00.2-5/prunData_numeric_pointer.desc"
myY <- read.table(phnotype,
                  header = TRUE, stringsAsFactors = FALSE)

# myGD <- read.big.matrix(numgd,
#                         type = "double", sep = "\t", header = TRUE,
#                         col.names = myGM$SNP, ignore.row.names = FALSE,
#                         has.row.names = TRUE)


# Load the data into a big.matrix.
# Note the use of the backingfile and descriptorfile arguments.

myGM <- read.table(numif,
                   header = TRUE, stringsAsFactors = FALSE)

myGD <- read.big.matrix(numgd,
                        type = "double", sep = "\t", header = TRUE,
                        col.names = myGM$SNP, ignore.row.names = FALSE,
                        has.row.names = TRUE, backingfile = "prunData_numeric.bin",
                        descriptorfile = "prunData_numeric.desc")

# Save the pointer for access later
dput(describe(myGD), "prunData_numeric_pointer.desc")

# The big.matrix can be reattached in a different R session using
# desc <- dget(pointer)
# myGD <- attach.big.matrix(desc)
for (i in 2:5)
{
  myResults <- farmcpu(Y = myY[,c(1,i)], GD = myGD, GM = myGM,ncores.glm=100 , ncores.reml=100)
  write_results(myResults)
  #manhattan_plot(myResults$COLEOPCOL_JI15$GWAS, cutoff = 0.01)
  rm(myResults)
}

