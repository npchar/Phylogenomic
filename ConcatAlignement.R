#!/usr/bin/env Rscript
# Exploring saturation :
#https://2infectious.wordpress.com/2014/06/17/concatenating-sequence-alignments/

library("optparse")
library("ape")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="SCOlist file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output alignement prefix", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("One argument must be supplied (input file name).n", call.=FALSE)
}



# Reading list of alignment to concatenate :
listALN = read.table(opt$file, header=FALSE)
bindAln = c()
PartFile = c()
for(i in seq(1, length(listALN$V1))){
	aln = listALN$V1[i]
	new = read.dna(as.character(aln), format="fasta")

	#Test if empty :
	if(dim(new)[2]==0){
		write.table(c(as.character(aln), "Is empty"), file=paste(as.character(aln),".WarningEmpty", sep=""), quote=F, row.names=F)
		next
	}
	row.names(new) = gsub("_.*", "", labels(new))
	if(i==1){ 
		bindAln = new
		PartFile = paste("DNA, gene",as.character(i),"=1-",as.character(dim(new)[2]), sep='') 
	}
	else{ 
		lastPos = dim(bindAln)[2]
		bindAln = cbind(bindAln, new, fill.with.gaps=TRUE,check.names=TRUE) 
		PartFile = c(PartFile, paste("DNA, gene",as.character(i),"=",as.character(lastPos + 1), "-",as.character(lastPos + dim(new)[2]), sep='' ))
	}
	
}


# concatenate the alignments, checking the names of the sequences in case samples order varies between alignments, and padding out the final alignment with gaps if any samples are missing from the individual alignments
#c=cbind(a,b,fill.with.gaps=TRUE,check.names=TRUE)

# have a look at your concatenated aln
pdf(file=paste(opt$out,"visualAlignment.pdf", sep='.'))
image.DNAbin(bindAln)
dev.off()

# write the concat aln to a new fasta, without any column spacers or anything else between nucleotides
write.dna(bindAln, file=paste(opt$out,"fasta", sep='.'), format="fasta", colsep="")
write.table(PartFile, file=paste(opt$out, "partition",sep='.'), quote=F, col.names=F, row.names=F)
