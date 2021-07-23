args = commandArgs(trailingOnly = TRUE)

library(dplyr)

file_input = args[1]    #"v44.3_1240K_public.anno"
file_input_2 = args[2]  #tempAFP.fam
region_input = tolower(args[3])  #"east_asia"
file_output = args[4] # "list_east_asia.txt"

#read in the file
anc_meta_raw <- read.delim(file_input, sep = "\t", header = T)

#anc_meta$Index means get the index column from the data
#produces a vector containing numbers from 1 to total number of rows; overwriting the original index column
#so the numbers are accurate
anc_meta_raw$Index <- 1:nrow(anc_meta_raw)

#rename the column names into easier ones
colnames(anc_meta_raw)[6] <- "Year.first.published"
colnames(anc_meta_raw)[10] <- "Date.mean"
colnames(anc_meta_raw)[11] <- "Date.sd"
colnames(anc_meta_raw)[12] <- "Full.date"
# not using the Y and mtDNA haplogroups yet, but these have some potential for helping to classify populations perhaps
colnames(anc_meta_raw)[25] <- "Y.haplogroup"
colnames(anc_meta_raw)[27] <- "mtDNA.haplogroup"

#remove white space
anc_meta_raw$Country <- trimws(anc_meta_raw$Country, which = c("right"))

#############################################
# modify the .fam file so it contains yearBins
#############################################
anc_fam <- read.delim(file_input_2, sep = "", header = F)

#rename the column names into ones that will ,match the .anno file
colnames(anc_fam)[1] <- "Index"
colnames(anc_fam)[2] <- "Version.ID"

anc_fam  <- left_join(anc_fam, anc_meta_raw, by = c("Index", "Version.ID")) 

#merge column 1 &2 to create a unique Index-Version.ID name
library(stringr)
stringr::str_c()
anc_fam <- anc_fam %>% dplyr::mutate(Version.ID = str_c(Index,Version.ID, sep="_" ))

#create a column called yearBins and sort each row into a time period
#yearBins replaces the original first column (which has been merged with column 2)
anc_fam$Index <- cut(anc_fam$Date.mean, breaks = seq(0,10000,by=500), right = TRUE)
colnames(anc_fam)[1]="yearBins"

#replace any "NA" in yearBins with "Modern"
anc_fam$yearBins <- as.character(anc_fam$yearBins) #turn the data from factors (as a result of cut) to characters
for (i in 1:nrow(anc_fam)){
  if (is.na( anc_fam$yearBins[i])){
    anc_fam$yearBins[i] <- "Modern"
  }}

anc_fam.6  <- anc_fam[,1:6]
#modify the .fam file
write.table(anc_fam.6, file = file_input_2, quote = F, row.names =F,
            col.names = F, sep = " ")
##########################################


#remove duplicates from the .anno file
anc_meta <- anc_meta_raw %>% group_by(Master.ID) %>% top_n(1,SNPs.hit.on.autosomal.targets) %>% ungroup()
#remove anyone over 10,000 year
anc_meta <- anc_meta %>% dplyr::filter(Date.mean < 10000)

#merge column 1 &2 to create a unique Index-Version.ID name
library(stringr)
stringr::str_c()
anc_meta_mutated <- anc_meta %>% dplyr::mutate(Version.ID = str_c(Index,Version.ID, sep="_" ))

#create a column called yearBins and sort each row into a time period
#yearBins replaces the original first column (which has been merged with column 2)
anc_meta_mutated$Index <- cut(anc_meta_mutated$Date.mean, breaks = seq(0,10000,by=500), right = TRUE)
colnames(anc_meta_mutated)[1]="yearBins"

#replace any "NA" in yearBins with "Modern"
anc_meta_mutated$yearBins <- as.character(anc_meta_mutated$yearBins) #turn the data from factors (as a result of cut) to characters
for (i in 1:nrow(anc_meta_mutated)){
  if (is.na( anc_meta_mutated$yearBins[i])){
    anc_meta_mutated$yearBins[i] <- "Modern"
  }}

#filter by region

#east Asians
if (region_input=="east_asia"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("China", "Japan", "Mongolia", "Taiwan", "Sri Lanka") )

#southeast asians
} else if(region_input=="southeast_asia") {
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("Vietnam", "Laos","Thailand", "Indonesia","Malaysia") )

#south asians  
} else if(region_input=="south_asia"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("India","Pakistan","Nepal", "South Korea", "Bangladesh") )

#all of asia
} else if(region_input=="asia"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("Indonesia","China","Malaysia",
                                                        "Laos","Sri Lanka","Mongolia",
                                                        "Taiwan","South Korea","Thailand",
                                                        "India","Pakistan","Vietnam",
                                                        "Japan","Nepal","Bangladesh") )
#all of europe
} else if(region_input=="europe"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("Russia", "Spain","Czech Republic", "Great Britain", 
                                                        "Hungary", "Estonia", "Portugal", "Turkey","Belgium", 
                                                        "Italy", "France", "Poland", "Iceland", "Finland", 
                                                        "Greece", "Armenia", "Norway", "Georgia", "Bulgaria", 
                                                        "Albania", "Croatia","Luxembourg", "Germany", "Romania", 
                                                        "Sweden", "Serbia", "North Macedonia", "Ukraine", "Denmark", 
                                                        "Montenegro", "Lithuania","Austria", "Switzerland","Moldova",
                                                        "Ireland", "Netherlands", "Latvia", "Slovakia", "Isle of Man"))
}else if(region_input=="western_europe"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("Spain","Great Britain","Portugal", 
                                                          "Italy", "France", "Iceland", "Finland", 
                                                          "Germany","Sweden", "Denmark", 
                                                          "Switzerland", "Ireland", "Netherlands", "Isle of Man", "Luxembourg" ))
}else if(region_input=="eastern_europe"){
  anc_region <- anc_meta_mutated %>% dplyr::filter(Country %in% c("Russia", "Czech Republic", "Hungary", "Estonia", 
                                                                  "Turkey", "Poland", "Armenia", "Georgia", "Bulgaria", 
                                                                  "Albania", "Croatia","Greece", "Romania", "Serbia", 
                                                                  "North Macedonia", "Ukraine", "Montenegro", "Lithuania",
                                                                  "Austria", "Latvia", "Slovakia"))
}

#get just the Index and Version.ID
listIndex <- anc_region %>% dplyr::select(yearBins, Version.ID)


#Write out into text files
write.table(listIndex, file = file_output, quote= F, col.names = F, row.names = F, sep = "\t")