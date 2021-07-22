#!/bin/bash

#Must provide an input file prefix, a region, and a SNP or a gene+GTF file

#Ex:  sh alleleFreqPlotter.sh -f v44.3_1240K_public -r EUROPE -g TYK2 -t Homo_sapiens.GRCh37.87.gtf
#Options for regions: EAST_ASIA, SOUTHEAST_ASIA, SOUTH_ASIA, ASIA, EUROPE, WESTERN_EUROPE, EASTERN_EUROPE

#$FILE is a prefix (ex v44.3_1240K_public)
#GTF file: Homo_sapiens.GRCh37.87.gtf

while getopts f:r:s:g:t: option
do
	case "${option}" in
		f) FILE=${OPTARG};;
		r) REGION=${OPTARG};;
		s) SNP=${OPTARG};;
		g) GENE=${OPTARG};;
		t) GTF=${OPTARG};;

	esac
done



# create a temp file of the v44.fam file
cp $FILE.fam temp.fam


# accepts a .anno file, and .fam file as input. Also takes a region argument, and outputs a "list to keep" text file containing the aDNA data for the specified region
Rscript alleleFunction.R $FILE.anno $FILE.fam $REGION list_$REGION.txt




##########
###SNPs###
##########

# if a SNP argument was given, run plink using the SNP
if [ -n "$SNP" ]
then

###PLINK####

# run plink using the v44 file, keeping the afformentioned "list to keep file". Uses the snp arguement, outputs a .frq.strat file which is named using the region and snp 
plink --bfile $FILE --keep list_$REGION.txt --snp $SNP --freq --family --missing --out $REGION_$SNP

###PLOTTING####

#plot the allele frequencies and output a pdf graph
Rscript allelePlotter.R ${REGION}_${SNP}.frq.strat




## if it's gene argument, filter gtf file, run bedtools, generate list_from_${GENE}.txt
## if you have the gene argument, you will need the program to require an argument for the gtf file


##########
###GENE###
##########

# otherwise, if a gene argument was given, run plink using the file of alleles (SNPs) for that gene

elif [ -n "$GENE"  ]
then

# STEP 1: generate a list of gene locations using a GTF file
awk '{ if ($3=="gene") print $0}' $GTF | grep "$GENE"  > filtered_${GENE}.gtf
 

# STEP 2: expand the gene range by 5000 in either direction for the gene previously searched for
awk '{OFS= "\t"; $4 = $4 - 5000; $5 = $5 + 5000 ; print}' filtered_${GENE}.gtf > adjusted_${GENE}.gtf


# STEP .5: #generate a .bed file containing the exact position of all SNPS from the v44.bim file(#$1 = chromosome number; $4 = position of SNP; $2 = rs number)
awk '{FS=" "; OFS="\t"; print $1, $4-1, $4, $2 > "extractedSNPs.bed"}' $FILE.bim


# STEP 3: do an intersect on the gene and the bed file created; save $4, which is the SNP number
bedtools intersect -a extractedSNPs.bed -b adjusted_${GENE}.gtf | awk '{print $4}'  > list_from_${GENE}_SNPs.txt

###PLINK####

#STEP 4: run plink using --extract instead of --SNP
plink --bfile $FILE --keep list_$REGION.txt --extract list_from_${GENE}_SNPs.txt --freq --family --missing --out ${REGION}_${GENE}

###PLOTTING###

#STEP 5: plot the allele frequencies and output graphs into a single pdf file
Rscript allelePlotter.R ${REGION}_${GENE}.frq.strat


###########
###error###
###########

#if no SNP and GENE+GTF file were provided, raise an error
else
echo "No SNP or GENE+GTF_file arguments provided"
exit 1 

fi



#restore the v44.fam file to it's original state (without yearbins)
mv temp.fam $FILE.fam
