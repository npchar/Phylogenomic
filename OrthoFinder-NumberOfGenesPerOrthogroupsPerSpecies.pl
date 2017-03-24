#!/usr/bin/perl
use strict;
use warnings;

#Description :
my $usage ;
$usage= "$0 <OG Matrix from Orthofinder>   [DEPRECATED] \n\n" ;
$usage.= "Description : This script is producing a report of the number of genes per orthogroup and per specie\n" ;
$usage.= "Dependency : <OG Matrix from Orthofinder> is a matrix with Orthogroup as row (i) and species as column (j)\n" ;
$usage.= "	       Mij is a list of genes separated by a comma (,) \n" ;


#my $OrthogroupFile = "Orthogroups.csv";
my $OrthogroupFile = $ARGV[0];
if (not defined $OrthogroupFile){
	die $usage ;
}

open IN, "<$OrthogroupFile" or die "Unable to open \"$OrthogroupFile\" !\n";
my %result ;
# Wrapping row by row ($i)
while(my $i = <IN>){
  chomp $i;
  my @sp = split "\t", $i ;
  print $sp[0] ;
  if($sp[0]!~ /^OG/){
      for my $j (1 .. $#sp) {print "\t",$sp[$j]}
  }
  else{
      # Specie by Specie ($j)
      for my $j (1 .. $#sp) {
	    my @seq = split /,/, $sp[$j] ;
	    my $NumSeq = @seq ;
	    print "\t",$NumSeq ;
      }
  }
  print "\n" ;
}
