#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Silix fnodes file name required", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Matrix output file name required", metavar="character")
  );
usage = "This R script will produce an abundance matrix (Number of sequences per family per sp) from a silix Fnodes file.\nWe assume that the species name is at the begining of the sequence name and separated by an underscore \'_\' (sp1_SeqName).\t" 
opt_parser = OptionParser(option_list=option_list, description=usage);
opt = parse_args(opt_parser);

if (is.null(opt$file) & is.null(opt$out)){
  print_help(opt_parser)
  stop("Two argument must be supplied (input and output files names)\n", call.=FALSE)
}

#Function :
inc <- function(x)
{
 eval.parent(substitute(x <- x + 1))
}



NODES = read.table(opt$file, h=F)
Name_SP = levels(as.factor(t(as.data.frame(strsplit(as.character(NODES$V2), "_")))[,1]))
Name_FAM = levels(as.factor(as.character(NODES$V1)))
MAT = as.data.frame(matrix(rep(0, c(length(Name_FAM)*length(Name_SP))),nrow=length(Name_FAM), ncol=length(Name_SP)), row.names=Name_FAM)
names(MAT) = Name_SP

for(i in seq(1,length(NODES$V1))){
	SPI = sub("_[0-9]*", "",NODES[i,]$V2)
	inc(MAT[rownames(MAT)==NODES[i,]$V1,names(MAT)==SPI])
}

write.table(MAT, quote=F, file=opt$out)
