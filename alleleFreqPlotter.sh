#!/bin/bash

#Must provide an input file prefix, a region, and a SNP or a gene+GTF file
#Final argument, -m is if want to plot moderns in final output or not (so -m y or -m n for yes or no)

#Ex:  sh alleleFreqPlotter.sh -f v44.3_1240K_public -r europe -g TYK2 -t Homo_sapiens.GRCh37.87.gtf -m y
#Options for regions: east_asia, southeast_asia, south_asia, asia, europe, western_europe, eastern_europe

#$FILE is a prefix (ex v44.3_1240K_public)
#ex GTF file: Homo_sapiens.GRCh37.87.gtf

while getopts f:r:s:g:t:m: option
do
	case "${option}" in
		f) FILE=${OPTARG};;
		r) REGION=${OPTARG};;
		s) SNP=${OPTARG};;
		g) GENE=${OPTARG};;
		t) GTF=${OPTARG};;
		m) MODERN=${OPTARG};;

	esac
done

#if no -m argument was provided, default to including moderns
if [ -n "$MODERN" ]
then
MODERN=$MODERN

else
MODERN="y"
fi


#make a copy of the .fam file for modifying
cp $FILE.fam tempAFP.fam


# input a .anno file, and .fam file, a region argument, and an output file name.
# outputs a "list to keep" text file which contains the year bin and the ID number for the SNPs in the region of interest

Rscript alleleFunction.R $FILE.anno tempAFP.fam $REGION list_$REGION.txt



##########
###SNPs###
##########

# if a SNP argument was given, run plink using the SNP
if [ -n "$SNP" ]
then

###PLINK####

# run plink using the v44 file, keeping the afformentioned "list to keep file". Uses the snp arguement, outputs a .frq.strat file which is named using the region and snp 
# use --fam to specify the tempAFP prefix is used instead of the $FILE prefix
plink --bfile $FILE --fam tempAFP.fam --keep list_$REGION.txt --snp $SNP --freq --family --missing --out ${REGION}_${SNP}

###PLOTTING####

#plot the allele frequencies and output a pdf graph
#second argument states whether to include moder samples in the final graph or not
Rscript allelePlotter.R ${REGION}_${SNP}.frq.strat ${MODERN}

mkdir output_${REGION}_${SNP}


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
plink --bfile $FILE --fam tempAFP.fam --keep list_$REGION.txt --extract list_from_${GENE}_SNPs.txt --freq --family --missing --out ${REGION}_${GENE}

###PLOTTING###

#STEP 5: plot the allele frequencies and output graphs into a single pdf file
Rscript allelePlotter.R ${REGION}_${GENE}.frq.strat ${MODERN}


mkdir output_${REGION}_${GENE}

###########
###error###
###########

#if no SNP and GENE+GTF file were provided, raise an error
else
echo "No SNP or GENE+GTF_file arguments provided"
exit 1 

fi


#remove intermediate files
#move output into the new output directory
rm list_$REGION.txt
rm tempAFP.fam

if [ -n "$SNP" ]
then
rm ${REGION}_${SNP}.nosex
rm ${REGION}_${SNP}.lmiss
rm ${REGION}_${SNP}.imiss
rm ${REGION}_${SNP}.log
rm ${REGION}_${SNP}.frq.strat

#move the output to a directory
mv ${REGION}_${SNP}* output_${REGION}_${SNP}

elif [ -n "$GENE"  ]
then

rm filtered_${GENE}.gtf
rm adjusted_${GENE}.gtf
rm extractedSNPs.bed
rm list_from_${GENE}_SNPs.txt

rm ${REGION}_${GENE}.nosex
rm ${REGION}_${GENE}.lmiss
rm ${REGION}_${GENE}.imiss
rm ${REGION}_${GENE}.log
rm ${REGION}_${GENE}.frq.strat

#move the output to a directory
mv ${REGION}_${GENE}* output_${REGION}_${GENE}

fi 


