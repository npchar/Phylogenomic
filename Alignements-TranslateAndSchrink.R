#!/usr/bin/env Rscript

library("optparse")
library("ape")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="file with individual alignement file names (or locations)", metavar="character"),
  make_option(c("-a", "--amino"), type="character", default=NULL,
              help="[Work in progress] Translate nt sequence with amino acid code default [1:standard]", metavar="character"),
  make_option(c("-r", "--remove"), type="character", default=NULL,
              help="Remove sites [3:remove third position]", metavar="character"),
  make_option(c("-o", "--prefix"), type="character", default=NULL,
              help="output alignement prefix", metavar="character"),
  make_option(c("-s", "--suffix"), type="character", default=".fasta",
              help="output alignement suffix", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) & is.null(opt$out)){
  print_help(opt_parser)
  stop("Two arguments must be supplied (input and output files names).n", call.=FALSE)
}

# Reading list of alignment to proceed:
#======================================
listALN = read.table(opt$file, header=FALSE)

# Initialize data structure:
#===========================
Infos = c("Name", "Nsp", "Nsites")

# Begin main Loop:
#=================
for(i in seq(1, length(listALN$V1))){
	# Open alignement:
	#-----------------
        aln = listALN$V1[i]
        new = read.dna(as.character(aln), format="fasta")

	# Extract basename of aln:
	#-------------------------
	AlnBNprov =  unlist(strsplit(as.character(aln), '/'))
	BNaln = AlnBNprov[length(AlnBNprov)]

	# get informations:
	#------------------
	Nsites = dim(new)[2]
	Nsp = dim(new)[1]
	Infos = rbind(Infos, c(BNaln, Nsp, Nsites))
	
	print(c(BNaln, Nsp, Nsites))

        # Translate DNA into AA if asked:
	#--------------------------------
        if(is.null(opt$amino)==FALSE & Nsites >= 3){
#                print_help(opt_parser)
#                stop("Work in progress... currently not functionnal option !")

                # print(opt$amino)
                ModAln = trans(new, code= opt$amino )
        }

        # remove Position:
	#-----------------
        if(is.null(opt$remove)==FALSE & Nsites >= 3){
                print(opt$remove)
                ModAln = new[,-seq(as.numeric(opt$remove), ncol(new), by=3)]
        }
	
	# write results:
	#---------------
	write.dna(ModAln, file=paste(opt$prefix,BNaln,opt$suffix, sep=''), format="fasta", colsep="")
	write.table(Infos, file=paste(opt$prefix, "TableAlignements.Infos",sep=''), quote=F, col.names=F, row.names=F)

}
