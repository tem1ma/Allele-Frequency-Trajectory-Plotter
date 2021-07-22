# Allele-Frequency-Trajectory-Plotter

This tool is designed to work with ancient genomic data archived in the [Allen Ancient DNA Resource](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data), which is maintained by Swapan Mallick and David Reich.

AlleleFreqPlotter.sh takes as input the genotype and metadata files from the database, and outputs a plot of minor allele frequencies over time for selected alleles from various geographic subdivisions.

It runs from the command line 

"sh alleleFreqPlotter.sh -f v44.3_1240K_public -r region -g gene -t Homo_sapiens.GRCh37.87.gtf"

3 or 4 arguments are required:

-f: a .anno file

-r: a region to search for (currently supports EAST_ASIA, SOUTHEAST_ASIA, SOUTH_ASIA, ASIA, EUROPE, WESTERN_EUROPE, EASTERN_EUROPE)

You can also provide the name of a gene of interest along with a gtf annotation file, if you are interested in plotting the frequencies of all alleles within and neighboring the gene. 

-g OR -s
-s: a SNP
-g: an entire gene (containing multiple SNPs); if this is provided, a -t argument must also be provided, which is a .gtf file containing gene information

Files and programs required:<br/>
>alleleFreqPlotter.sh<br/>
a .fam file<br/>
a .anno file<br/>
alleleFunction.R<br/>
PLINK<br/>
allelePlotter.R<br/>
a.gtf file<br/>


