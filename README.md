# BacteroidesGA1

##########################################################################################################
##### DESCRIPTION                                                                                    #####
##### Code for preliminary processing of .fastq files from Illumina sequencing. Processing includes: #####
##### 1. standard DADA2 pipeline                                                                     #####
##### 2. assigning taxonomy using a Bacteroides-specific database                                    #####
##### 3. passing data into phyloseq object                                                           #####
##### 4. correcting a few mis-identified important taxa                                              #####
##### 5. adding columns to mapping file based on calculations from processed sequencing data         #####
##### 6. identify and remove contaminants using decontam                                             #####
##### Use this file to prepare initial raw phyloseq object for all future downstream analysis.       #####
###########################################################################################################################################################
##### This code requires as input:                                                                                                                    #####
##### 1. A working repository (default "~/Desktop/Wannigan/")                                                                                   #####
##### 2. Folder within working repository containing all .fastq files (default "~/Desktop/Wannigan/raw/16S/")                                   #####
##### 3. .csv file containing sample metadata listed with unique sampleIDs in first column (default "~/Desktop/Wannigan/raw/table_map.csv")     #####
##### 4. .fa file containing reference taxonomy database (default "~/Desktop/Wannigan/raw/dbBacteroides.fa")                                     #####                                   
###########################################################################################################################################################

dbBacteroides.fa is included in repository
