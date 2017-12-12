#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name (generally named Orthogroups.csv)", metavar="character"),
	make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output prefix file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) | is.null(opt$out)){
  print_help(opt_parser)
  stop("Two arguments must be supplied (input and output file name).n", call.=FALSE)
}


data = read.csv(opt$file, row.names=1, sep="\t", h=T)
# data2 = data[,2:length(names(data))]
# row.names(data2) = data[,1]
# looks working : (split by "," and thenn count number of occurance)
#length(unlist(strsplit(as.character(data2[1,1]), ",")))

CountSeq = function(x){
  return(length(unlist(strsplit(as.character(x), ","))))
}
occur = apply(data, MARGIN=c(1,2),FUN = CountSeq)
write.table(file=opt$out , occur, row.names=T, col.names=NA, quote=F, sep="\t")
