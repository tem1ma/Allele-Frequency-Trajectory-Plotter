# Allele-Frequency-Trajectory-Plotter

This tool is designed to work with ancient genomic data archived in the [Allen Ancient DNA Resource](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data) (AADR), which is maintained by Swapan Mallick and David Reich.

AlleleFreqPlotter.sh takes as input the genotype and metadata files from the AADR, and outputs a plot of minor allele frequencies over 1000-year sliding windows over the last 10,000 years for selected alleles from various provided custom geographic subdivisions.

The tool is run from the command line and requires 3 or 4 arguments:

**-f:** a *.anno* file

**-r:** a region to filter for (currently supports *EAST_ASIA, SOUTHEAST_ASIA, SOUTH_ASIA, ASIA, EUROPE, WESTERN_EUROPE, EASTERN_EUROPE*)

**-s**: a SNP ID
```
sh alleleFreqPlotter.sh -f v44.3_1240K_public -r europe -s rs34536443
```

or<br/>
**-g**: a gene name; if this is provided, a **-t** flag must also be provided with a *.gtf* annotation file containing gene information  
- Use if investigating multiple SNPs in a gene

```
sh alleleFreqPlotter.sh -f v44.3_1240K_public -r europe -g TYK2 -t Homo_sapiens.GRCh37.87.gtf
```


Files and programs required:<br/>
>alleleFreqPlotter.sh<br/>
a *.fam* file<br/>
a *.anno* file<br/>
alleleFunction.R<br/>
PLINK<br/>
allelePlotter.R<br/>
a *.gtf* file<br/>


