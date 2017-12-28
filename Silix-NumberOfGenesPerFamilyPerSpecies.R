#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Silix fnodes file name required", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Matrix output file name required", metavar="character")
  );
usage = "This R script will produce an abundance matrix (Number of sequences per family per sp) from a silix Fnodes file.\nWe assume that
 the species name is at the begining of the sequence name and separated by an underscore \'_\' (sp1_SeqName).\n It will also report the n
umber of SCO (with different occupancy thresholds) as well as the number of Family." 
opt_parser = OptionParser(option_list=option_list, description=usage);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) & is.null(opt$out)){
  print_help(opt_parser)
  stop("Two argument must be supplied (input and output files names)\n", call.=FALSE)
}


inc <- function(x)
{
 eval.parent(substitute(x <- x + 1))
}

NODES = read.table(opt$file, h=F)
Name_SP = levels(as.factor(t(as.data.frame(strsplit(as.character(NODES$V2), "_")))[,1]))
Name_FAM = levels(as.factor(as.character(NODES$V1)))

# Extract Silix parameters from the filename :
#---------------------------------------------
#Assuming nammed as follow : silix-allTsa.i0.3_r0.7_l250_m0.5.fnodes
n=opt$file
Words = strsplit(n, "_")
Sidentity = as.numeric(sub("^.*silix-allTsa.i","",Words[[1]][1]))
Soverlap = as.numeric(sub("r","",Words[[1]][2]))
Sminlength = as.numeric(sub("l","",Words[[1]][3]))
Sminoverlap = as.numeric(sub(".fnodes","",sub("m","",Words[[1]][4])))


#Create Matrix from Fnodes file :
#--------------------------------
MAT = as.data.frame(matrix(rep(0, c(length(Name_FAM)*length(Name_SP))),nrow=length(Name_FAM), ncol=length(Name_SP)), row.names=Name_FAM)
names(MAT) = Name_SP
for(i in seq(1,length(NODES$V1))){
        SPI = sub("_[0-9]*", "",NODES[i,]$V2)
        inc(MAT[rownames(MAT)==NODES[i,]$V1,names(MAT)==SPI])
}
#Dump matrix
write.table(MAT, quote=F, file=opt$out)

# Compute some metrics :
#-----------------------
#Number of SCO (different level of missing values) :
SCO = MAT[apply(MAT[,1:length(MAT[1,])], 1, function(x) sum(x<2)==length(x)),]
SCO100 = dim(SCO[apply(SCO[,1:length(SCO[1,])], 1, function(x) sum(x)== length(x)),])[1]
SCO75 = dim(SCO[apply(SCO[,1:length(SCO[1,])], 1, function(x) sum(x)> length(x)*0.75),])[1]
SCO50 = dim(SCO[apply(SCO[,1:length(SCO[1,])], 1, function(x) sum(x)> length(x)*0.5),])[1]
SCO375 = dim(SCO[apply(SCO[,1:length(SCO[1,])], 1, function(x) sum(x)> length(x)*0.375),])[1]
NumberFam=dim(MAT)[1]

header=c("FileName", "identity", "overlap", "minlength", "minoverlap", "SCO100", "SCO75", "SCO50", "SCO375", "Nfam")
body=c(n, Sidentity, Soverlap, Sminlength, Sminoverlap, SCO100, SCO75, SCO50, SCO375, NumberFam)
InfosFile=paste(opt$out, '.metrics', sep='')

write.table(rbind(header,body), quote=F, row.names=F, col.names=F,file=InfosFile)

