# Phylogenomic
## A) OrthoFinder tools :
### i) OrthoFinder-NumberOfGenesPerOrthogroupsPerSpecies.pl
This first perl script is working on OrthoFinder ortohogroup matrix "Orthogroups.csv" :
##### Usage :
```
./OrthoFinder-NumberOfGenesPerOrthogroupsPerSpecies.pl <OG Matrix from Orthofinder>

Description : This script is producing a report of the number of genes per orthogroup and per specie
Dependency : <OG Matrix from Orthofinder> is a matrix with Orthogroup as row (i) and species as column (j)
               Mij is a list of genes separated by a comma (,)
```
