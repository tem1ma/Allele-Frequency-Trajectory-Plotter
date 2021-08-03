args = commandArgs(trailingOnly = TRUE)
library(ggplot2)
library(epitools)
library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)


file_input = args[1] # ex europe_TYK2.frq.strat" or "europe_rs34536443.frq.strat"
incl_modern = tolower(args[2])  #y or n regarding whether modern samples should be output in the final plot

#gene here is used to mean the gene, or the SNP if only 1 SNP was investigated
#file input contains SNPs in $2--if just one SNP was input, this will be the same for the whole file
#if a gene was input, $2 will contain multiple SNPs

#extract the gene name from the input file
region_gene_name = gsub(".frq.strat","",file_input)
#split on the _ character
split_name = unlist(strsplit(region_gene_name, "_"))
#extract the region name (if compound name like east_asia, paste back together)
region_name = c(split_name[1:length(split_name)-1])
region_name = paste(region_name, collapse = "_")
#extract the gene name
gene_name = split_name[length(split_name)]

freqFile <-read.delim(file_input, sep = "", header = T)

#put the axis in order
freqFile$CLST <- factor(freqFile$CLST, levels = c("Modern","(0,500]", "(500,1e+03]",
                                                  "(1e+03,1.5e+03]",
                                                  "(1.5e+03,2e+03]", "(2e+03,2.5e+03]",
                                                  "(2.5e+03,3e+03]", "(3e+03,3.5e+03]",
                                                  "(3.5e+03,4e+03]", "(4e+03,4.5e+03]",
                                                  "(4.5e+03,5e+03]", "(5e+03,5.5e+03]",
                                                  "(5.5e+03,6e+03]", "(6e+03,6.5e+03]",
                                                  "(6.5e+03,7e+03]", "(7e+03,7.5e+03]",
                                                  "(7.5e+03,8e+03]", "(8e+03,8.5e+03]",
                                                  "(8.5e+03,9e+03]", "(9e+03,9.5e+03]",
                                                  "(9.5e+03,1e+04]"))
#put the dataframe in order on the CLST column
frq_mutated <- freqFile %>% dplyr::arrange(SNP,CLST)

#extract all the rows containing modern samples
#future code excludes them, so df_modern will be used to add them back at the end for plotting
df_modern <- data.frame()
for (i in 1:nrow(frq_mutated)){
  if (frq_mutated$CLST[i] == "Modern"){
    df_modern <- rbind(df_modern, frq_mutated[i,])
  }}

#divide the NCHROBS and MAC columns by 2 (except for the moderns)
for (i in 1:nrow(frq_mutated)){
  if (frq_mutated$CLST[i] !="Modern"){
    frq_mutated$NCHROBS[i] <- (frq_mutated$NCHROBS[i]) /2
    
    if (frq_mutated$MAC[i] !=0){
      frq_mutated$MAC[i] <- frq_mutated$MAC[i] /2
    }}}

#remove the (] brackets
frq_mutated$CLST <- gsub("(","",frq_mutated$CLST, fixed = T)
frq_mutated$CLST <- gsub("]","",frq_mutated$CLST, fixed = T)

#split the CLST date column into two columns "from" and "to"
frq_mutated <- separate(frq_mutated, col = "CLST", into = c("from","to"), sep = ",")

#adjust the Modern column to now be from 0 to 0
frq_mutated$from[is.na(frq_mutated$from)] <- 0
frq_mutated$to[is.na(frq_mutated$to)] <- 0


#change the from and to values into numeric values
frq_mutated$from <- as.numeric(frq_mutated$from)
frq_mutated$to <- as.numeric(frq_mutated$to)


#create a dataframe with the "from" and "to" values
ranges  <- data.frame(from = seq(0,9500, by =  500), to = seq(500,10000, by = 500) )

# split here into list
list_by_snp <- split(frq_mutated, frq_mutated$SNP)

# incorporate left_join after splitting using the from and to column
for (i in 1:length(list_by_snp)){
  list_by_snp[[i]] <- left_join(ranges, list_by_snp[[i]], by = c("from","to"))
}


#replace any NAs in the MAC & NCHROBs columns with zeros
for (i in 1:length(list_by_snp)){
  list_by_snp[[i]]$NCHROBS <- replace(list_by_snp[[i]]$NCHROBS, is.na(list_by_snp[[i]]$NCHROBS), 0)
  list_by_snp[[i]]$MAC <- replace(list_by_snp[[i]]$MAC, is.na(list_by_snp[[i]]$MAC), 0)
}
# Adjusting the MAF
# MAF = (MAC1 + MAC2) /( NCHROBS1 + NCHROBS2)
# Create the sliding windows

updated_list_by_snp <- list()


#for every dataframe in the list
for (i in 1:length(list_by_snp) ){

    new_from <- c()
    new_to <- c()
    new_mac <- c()
    new_nchrobs <- c()
    
    #get the dataframe for the i'th position
    df <- list_by_snp[[i]]
  
    for (j in 1:(nrow(df) - 1) ){
        new_from <- c(new_from, df$from[j])  #keep "from" the same
        new_to <- c(new_to, df$to[j+1])      #adjust "to" to be the next row's "to"
        new_mac <- c(new_mac, sum(df$MAC[j:(j+1)]))  #add the MAC's for two consecutive rows together and append to vector
        new_nchrobs <- c(new_nchrobs, sum(df$NCHROBS[j:(j+1)]))  #do the same for NCHROBS
    }
    #create a new dataframe with the sliding windows and new from, to, MAC, and NCHROBS
    df_new <- data.frame(SNP=df$SNP[1:19],new_from,new_to,new_mac,new_nchrobs)
    
    #create a new column called adjusted_mac. In this column, replace any 0 with 1 
    #(this will be used when creating an uncertainty ribbon so uncertainty is better reflected)
    df_new <- df_new %>% dplyr::mutate(adjusted_mac = new_mac )
    df_new <- df_new %>% dplyr::mutate(adjusted_nchrobs = new_nchrobs )
    df_new$adjusted_mac[df_new$adjusted_mac==0] <- 1
    df_new$adjusted_nchrobs[df_new$adjusted_nchrobs==0] <- 1
    
    #get the new MAF values
    df_new <- df_new %>% dplyr::mutate(new_maf = new_mac / new_nchrobs ) 
    #replace any NaNs with 0s
    df_new$new_maf[is.nan(df_new$new_maf)] <- 0
    
    #combine the characters to create a new column, date.range
    df_new <- df_new %>% dplyr::mutate(Date.range = str_c(new_from,new_to,sep="-"))
    
    df_new$Date.range <- factor(df_new$Date.range, levels = c("0-1000", "500-1500", "1000-2000", "1500-2500",
                                                      "2000-3000",  "2500-3500",  "3000-4000",  "3500-4500",
                                                      "4000-5000",  "4500-5500",  "5000-6000",  "5500-6500", "6000-7000",
                                                      "6500-7500", "7000-8000",  "7500-8500",  "8000-9000",
                                                      "8500-9500",  "9000-10000"))
    
    #get a vector of all the upper and lower uncertainty values
    lower_uncert <- c()
    upper_uncert <- c()
    
    #loop through each row and get the uncertainty values for each row
    #binom.approx outputs lower uncertainty in the 4th position and upper in the 5th
    for (k in 1:(nrow(df_new)) ){
      lower_uncert <- c(lower_uncert, binom.approx(df_new$adjusted_mac[k], df_new$adjusted_nchrobs[k], conf.level = 0.95)[,4])
      upper_uncert <- c(upper_uncert, binom.approx(df_new$adjusted_mac[k], df_new$adjusted_nchrobs[k], conf.level = 0.95)[,5])
    }
    
    df_new <- df_new %>% dplyr::mutate(lower.uncert = lower_uncert, upper.uncert = upper_uncert)
    
    #add this new df to the new list, name that position after the SNP
    SNP_name <- df_new$SNP[1]
    updated_list_by_snp[[SNP_name]] <- df_new
}

#replace any NAs in the SNP column with the SNP name 
#every dataframe in updated_list_by_snp contains info about just one SNP, so replace all
#NAs in that column with that frame's SNP name
for (i in 1:length(updated_list_by_snp)) {
  updated_list_by_snp[[i]]$SNP <- names(updated_list_by_snp)[i]
}

###############################
#adding back the modern columns
###############################

#in the dataframe of moderns, add columns for upper and lower uncertainties
m_lower_uncert <- c()
m_upper_uncert <- c()

#loop through each row and get the uncertainty values for each row
#binom.approx outputs lower uncertainty in the 4th position and upper in the 5th
for (i in 1:(nrow(df_modern)) ){
  m_lower_uncert <- c(m_lower_uncert, binom.approx(df_modern$MAC[i], df_modern$NCHROBS[i], conf.level = 0.95)[,4])
  m_upper_uncert <- c(m_upper_uncert, binom.approx(df_modern$MAC[i], df_modern$NCHROBS[i], conf.level = 0.95)[,5])
}

#add the uncertainties to df_modern
df_modern <- df_modern %>% dplyr::mutate(lower.uncert = m_lower_uncert, upper.uncert = m_upper_uncert)
df_modern$CLST <- as.character(df_modern$CLST)

#create a new data frame with col names that match the ones in updated_list_by_snp
#adjusted_mac & nchrobs are just the same as new_mac &nchrobs
fixed_df_modern <- data_frame("SNP" = df_modern$SNP, "new_from" = df_modern$CLST, 
                              "new_to" = df_modern$CLST, "new_mac" = df_modern$MAC,
                              "new_nchrobs" = df_modern$NCHROBS, "adjusted_mac" = df_modern$MAC,
                              "adjusted_nchrobs" = df_modern$NCHROBS, "new_maf" = df_modern$MAF,
                              "Date.range" = df_modern$CLST, "lower.uncert" = df_modern$lower.uncert,
                              "upper.uncert" = df_modern$upper.uncert)


for (i in 1:nrow(fixed_df_modern)){
  pulled_SNP <- fixed_df_modern$SNP[i]
  updated_list_by_snp[[pulled_SNP]] <- rbind(updated_list_by_snp[[pulled_SNP]], fixed_df_modern[i,])
  
  #if moderns shouldn't be included in the final graph, exclude it from the list of factors
  if (incl_modern=="n"){
    updated_list_by_snp[[pulled_SNP]]$Date.range <- factor(updated_list_by_snp[[pulled_SNP]]$Date.range,
                                                           levels = c("0-1000", "500-1500", "1000-2000", "1500-2500",
                                                                      "2000-3000",  "2500-3500",  "3000-4000",  "3500-4500",
                                                                      "4000-5000",  "4500-5500",  "5000-6000",  "5500-6500", "6000-7000",
                                                                      "6500-7500", "7000-8000",  "7500-8500",  "8000-9000",
                                                                      "8500-9500",  "9000-10000"))
  }else{
  updated_list_by_snp[[pulled_SNP]]$Date.range <- factor(updated_list_by_snp[[pulled_SNP]]$Date.range,
                                                         levels = c("Modern", "0-1000", "500-1500", "1000-2000", "1500-2500",
                                                            "2000-3000",  "2500-3500",  "3000-4000",  "3500-4500",
                                                            "4000-5000",  "4500-5500",  "5000-6000",  "5500-6500", "6000-7000",
                                                            "6500-7500", "7000-8000",  "7500-8500",  "8000-9000",
                                                            "8500-9500",  "9000-10000"))
}}

#output all the graphs into a single PDF
pdf(paste0(region_gene_name,".pdf"), onefile = TRUE)

#generate the plot for each data frame stored in the list
#and create a text file of the data used to generate the plot
for (frame in updated_list_by_snp){
  #get the name of the current SNP
  SNP_name <- frame$SNP[1]
  #name the file after the gene and the SNP (unless only a SNP was given)
  #if gene name = SNP name, it means only a SNP was input originally, not a gene
  if (gene_name != SNP_name){
    output_name <- str_c(region_name, gene_name, SNP_name, sep="_")
  }else{
    output_name <- str_c(region_name, SNP_name, sep="_")
  }
  
  #save the df used into a text file
  #make sure it outputs in levels order
  orderedFrame <- frame[order(frame$Date.range),]
  write.table(orderedFrame, file = paste0(output_name, ".txt"), quote = F, row.names =F, col.names = T, sep = "\t")
  
  #if moderns should not be included in the final graph, remove them from the df
  if (incl_modern=="n"){ 
  new_frame <- frame[-20,]
  }else{
  new_frame <- frame
  }
 

  my_plot <- ggplot(new_frame, aes(x = Date.range, y = new_maf , group = SNP)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.1, lty = 5, col = "red") +
    xlab("Years (BP)") +
    ylab("Allele frequency") +
    ggtitle(output_name) + 
    geom_ribbon(aes(ymin = lower.uncert,  ymax = upper.uncert), alpha = 0.1) +
    #label=paste0("n=",new_frame$new_nchrobs)
    geom_text(label=new_frame$new_nchrobs, nudge_y = 0.010, check_overlap = F) +
    scale_x_discrete(limits = rev(levels(new_frame$Date.range))) +
    theme(axis.text.x = element_text(angle = 45),
          panel.grid = element_blank(),
          axis.text.x.bottom = element_text(hjust = 1),
          panel.border = element_rect(colour = "black",fill = NA)) 
  
  print(my_plot)
  
  
}
dev.off()
