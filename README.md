# Phylogenomic
## A) OrthoFinder tools :
### i) OrthoFinder-NumberOfGenesPerOrthogroupsPerSpecies.pl
This first R script is working on OrthoFinder ortohogroup matrix "Orthogroups.csv" :
##### Usage :
```
./OrthoFinder-CountOrthogroupsPerSpecies.R 
Usage: ./OrthoFinder-CountOrthogroupsPerSpecies.R [options]


Options:
        -f CHARACTER, --file=CHARACTER
                dataset file name

        -o CHARACTER, --out=CHARACTER
                output prefix file name

        -h, --help
                Show this help message and exit

```
Description : This script is producing a report of the number of genes per orthogroup and per species
Dependency : <OG Matrix from Orthofinder> is a matrix with Orthogroup as row (i) and species as column (j)
               Mij is a list of genes separated by a comma (,)

### ii) OrthoFinder--SCOmatrixBuilder.R
This script is filtering orthogroups observed one and only one times per species (Single Copy Orthologs).
It needs an ocurrence matrix (as printed by OrthoFinder-NumberOfGenesPerOrthogroupsPerSpecies.pl), as well as three R library ('ggplot2', 'reshape2' and 'optparse')
```
./OrthoFinder-SCOmatrixBuilder.R 
Usage: ./OrthoFinder-SCOmatrixBuilder.R [options]


Options:
        -f CHARACTER, --file=CHARACTER
                dataset file name

        -o CHARACTER, --out=CHARACTER
                output prefix file name

        -h, --help
                Show this help message and exit
```
Two types of files is produced : a list of orthogroups ids in "*.list" and a graphical representation of genes presence/absence by species.
This two files is produced for SCO observed at least in 75% of species, 50% and 37.5%.

