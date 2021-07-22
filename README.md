# Allele-Frequency-Trajectory-Plotter

This function takes raw data ancient human remain samples and outputs a plot of minor allele frequencies over time. 
It runs from the command line "sh alleleFreqPlotter.sh -f v44.3_1240K_public -r region -g gene -t Homo_sapiens.GRCh37.87.gtf"

3 or 4 arguments are required:

-f: a .anno file

-r: a region to search for (currently supports EAST_ASIA, SOUTHEAST_ASIA, SOUTH_ASIA, ASIA, EUROPE, WESTERN_EUROPE, EASTERN_EUROPE)

-g OR -s
-s: a SNP
-g: an entire gene (containing multiple SNPs); if this is provided, a -t argument must also be provided, which is a .gtf file containing gene information

Files and programs required:
> alleleFreqPlotter.sh 
	> v44.3_1240K_public.fam
	> v44.3_1240K_public. anno
	> alleleFunction.R
	> PLINK
	> allelePlotter.R
	> Homo_sapiens.GRCh37.87.gtf
