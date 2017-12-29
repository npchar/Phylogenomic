#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;



# Creation d'une matrice avec en ligne les SCO et en colonnes les sp : (matrice Orthogroup.csv d'Orthofinder)

# Besoin de quatres fichiers :
# i) correspondance codesp -> Sp
# ii) correspondance NumSeq -> id
# iii) silix.fnodes
# iv) liste des SCO d'interet 

# Exemple file :
# -s /data/Multispe/SCO/C95Lv2/OrthoFinder/Results_Mar24/WorkingDirectory/SpeciesIDs.txt 
# -t /data/Multispe/SCO/C95Lv2/OrthoFinder/Results_Mar24/WorkingDirectory/SequenceIDs.txt
# -l /data/Multispe/SCO/C95Lv2/SilixSCO/silix-allTsa.i0.7SCO75.list
# -n /data/Multispe/SCO/C95Lv2/SilixSCO/silix-allTsa.i0.7_r0.9_l300_m0.9.fnodes

my $Spdict ;
my $list ;
my $Tdict ;
my $help ;
my $fnodes ;
my $out = '_MAT' ;
my $verbose ;
my $Shortening = "TRUE" ;
my $outputName ;

GetOptions(
           's=s'          => \$Spdict,
	   'l=s'	=>\$list,
           't=s'         => \$Tdict,
           'h'        => \$help,
	   'o=s'	=>\$outputName,
           'n=s'        =>\$fnodes,
	   'v'		=>\$verbose,
          );



my $usage = "\nUsage: $0 -s <species dictionary> -t <transcripts dictionary> -n <silix fnodes> -l <list of SCO> -h\n";
$usage .= "Description : Build an orthogroup matrix (OrthoFinder like)\n";
$usage .= "Option :    -l   SCO of interest\n";
$usage .= "            -s   <species dictionary> from Orthofinder \n";
$usage .= "            -n   <silix fnodes> linking a transcript to a family\n";
$usage .= "            -t   <transcripts dictionary> from Orthofinder \n";
$usage .= "            -h   Displays this help and exit\n";
$usage .= "	       -o   Matrix file name [<silix fnodes>_MAT]\n" ;
$usage .= "	       -v   Verbose mode : usefull to see what's happening there as well as error ! \n\n" ;


if($help){
    print $usage;
    exit 0;
}

unless(defined($Spdict) && -f $Spdict && defined($Tdict)  && -f $Tdict && defined($list) && -f $list  && defined($fnodes) && -f $fnodes){
    print STDERR "File(s) [$Spdict], [$list] or [$Tdict] [$fnodes] not accessible : $!\n";
    print STDERR $usage;
    exit 1;
}

# Open and store list and dictionary
#===================================
#----list----
open(LIST, "<$list") or die "Unable to open $list !" ;
my %SCOlist ;
while (my $line = <LIST>){
	chomp $line ;
	$SCOlist{$line} = "TODO" ;
	
}
close(LIST) ;

#----SP----
open(SP, "<$Spdict" ) or die "Unable to open $Spdict !" ;
my %DictSP ;
while (my $line = <SP>){
	chomp $line ;
	my @words = split ': ', $line ;
	$DictSP{$words[0]} = $words[1] ;
}
close(SP) ;

#----Transcripts----
open(T, "<$Tdict" ) or die "Unable to open $Spdict !" ;
my %DictT ;
while (my $line = <T>){
	chomp $line ;
	my @words = split ': ', $line ;
	$DictT{$words[0]} = $words[1] ;
}
close(T) ;


# Open Fnodes file and proceed if SCO of interest
#================================================
# Fill the %OccurMat data structure:
# It stores transcripts names associated with a SCO, for every species:
my %OccurMat ;
 
open(FNODES, "<$fnodes" ) or die "Unable to open $fnodes !" ;
while (my $line = <FNODES>){
	chomp $line ;
	my @words = split /\s/, $line ;
	my $SCO = $words[0] ;
	my $transcriptCode = $words[1] ;
	my $Transcript ;

	# Skip SCO if not in the list (Speed improvement !!!)
	unless(exists $SCOlist{$SCO}){next ;}

	# Translate Transcripts
	if(exists $DictT{$transcriptCode}){
		$Transcript = $DictT{$transcriptCode} ;
	} 
	else{ 
		if($verbose){ print "Error with $transcriptCode : no entry in Transcript dictionnary\n" ;}
		 $Transcript = '' ;
	}
	
	#Find the species :
	my $CurrentSp ;
	my @Sp_transcript = split '_', $transcriptCode ;
	if(exists $DictSP{"$Sp_transcript[0]"}){
		$CurrentSp = $DictSP{"$Sp_transcript[0]"} ;		
#		print $transcriptCode . "\t" . $CurrentSp . "\n" ;
	}
	else{
		if($verbose){ print "Error with $transcriptCode : no entry in Sp dictionnary\n"; }
		$CurrentSp = '' ;
	}

	
	# Store the information in the data structure: Initialize it if not.
	# Build the hash even for species which doesn't have an SCO (It will prevent future "Undef")
	unless(exists $OccurMat{$SCO}){
		foreach my $s (sort keys %DictSP){
			$OccurMat{$SCO}{$DictSP{$s}} = '.' ;
		}
	}

	#Store the transcript name for the SCO, SP :
	$OccurMat{$SCO}{$CurrentSp} = $Transcript ;
}
close(FNODES) ;

# Print The Matrix file :
my $outputFile ;
if(defined $outputName){ $outputFile = $outputName}
else{$outputFile = $fnodes . $out ;}
open(OUT, ">$outputFile") or die "Unable to open $outputFile !\n" ;

#Define an order of species, store it and print it on a Header:
my $count=0 ;
my @SpOrder ;
foreach my $sco (sort keys %OccurMat){
	if($count eq "1"){last;}
	print OUT "X" ;
	foreach my $sp (keys $OccurMat{$sco}){
		my $Spname = $sp ;
		push @SpOrder, $Spname ;
		$Spname =~ s/\.\w*$// ;
		print OUT "\t" . $Spname;
	}
	print OUT "\n" ;
	$count = "1" ;
}


#Print the Matrix:
foreach my $sco (sort keys %OccurMat){
	print OUT $sco ;
	for my $i (0 .. $#SpOrder){
		my $sp =$SpOrder[$i] ;
		print OUT "\t" . $OccurMat{$sco}{$sp} ; # . $OccurMat{$sco}{$sp} ;
	}
	print OUT "\n" ;
}
close OUT ;
