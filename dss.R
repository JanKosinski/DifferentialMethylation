#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DSS")
#install.packages("argparse", repos='http://cran.us.r-project.org')
library(DSS)
library(argparse)
require(bsseq)

#preparing bismark.cov files
#awk -v OFS="\t" 'BEGIN{print("chr", "pos", "N", "X")}{print($1, $2, $5+$6, $5)}' KJ078af3_Tn.bismark.cov > KJ078af3_Tn.tsv

parser <- ArgumentParser()
parser$add_argument('-c', '--tc', type="character", nargs='+', dest="tc", help='Tp samples for diferential methylation analysis')
parser$add_argument('-n', '--tn', type="character", nargs='+', dest="tn", help='Tn samples for diferential methylation analysis')
parser$add_argument('-d', '--dir', type="character", dest="directory", help='Directory where input files are stored')
args <- parser$parse_args()
# todo dir should have / at the end

if (length(args$tc)!=length(args$tn)){
  print("Error")
}

files <- list()
# labelsTn contains only Tn samples, labelsTc only Tc samples and labels combined Tn and Tc samples so that the order in files vector corresponds to the order in labels
labelsTn <- labelsTc <- labels <- c()

for (i in 0:(length(args$tn)-1)){
    path <- paste(args$directory,args$tn[i+1],".tsv", sep="")
    files[[i*2+1]] <- read.table(path, header=TRUE)
    labelsTn[i] <- paste("Tn", i, sep = "")
    labels[i*2+1] <- paste("Tn", i, sep = "")
    path <- paste(args$directory,args$tc[i+1],".tsv", sep="")
    files[[i*2+2]] <- read.table(path, header=TRUE)
    labelsTc[i] <- paste("Tc", i, sep = "")
    labels[i*2+2] <- paste("Tc", i, sep = "")
}

path <- file.path(system.file(package="DSS"), "extdata")
BSobj <- makeBSseqData( files, labels )

#(1) estimate mean methylation levels for all CpG site; (2) estimate dispersions at each CpG sites; (3) conduct Wald test.
#authors recommend to use smoothing with whole-genome BS-seq
dmlTest <- DMLtest(BSobj, group1=labelsTn, group2=labelsTc, smoothing=TRUE)    

#differential methylation loci
dmls <- callDML(dmlTest, p.threshold=0.001)
write.table(dmls, paste(args$directory, "dmls.txt", sep="\t"))

#differential methylation regions
dmrs <- callDMR(dmlTest, p.threshold=0.001)
write.table(dmrs, paste(args$directory, "dmrs.txt", sep="\t"))

#coverage depth information plot
showOneDMR(dmrs[1,], BSobj)
