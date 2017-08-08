# Preparing data for analysis and creating moisture indexes #####
# owner Gabriel Massaine Moulatlet
# created in November of 2016
# this script is to prepara the data to analyse in the Hydro paper of Jurua                                             
# it is divided by groups: gingers, ferns and palms  because original tables
# were too different. Each dataset required a special treatment.
# as outputs, it creates a wide-format table for community analysis, and
# tables with a moisture index
# don't run the script per parts, but the whole thing. Files will be overwrite, but no worries :)


#load the packages

library(reshape2)
library(vegan)

################# zingibs #########################
# organizing the zingib data
# importing the data and set the working directory

                getwd()
                setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados floristicos originais")
                setwd("C:/workspace gabriel/hidrologia do jurua/Analyses/dados floristicos originais")
                zingdata <- read.csv("Zing_jurua_database.csv",stringsAsFactors=FALSE)
                
                    str(zingdata)
                    zingdata$quant <- as.numeric(zingdata$quant)   # transform the species abundances in numeric
                    zingdata$id_1  <- tolower(zingdata$id_1)       # transform species id to lowercase to avoid spelling errors
                    
                    #fix zings subunits numbers: from 1:100 to 1:20
                    newsubunitnamezing <- matrix(seq(1,100), ncol = 20) # create a matrix with numbers to be related
                    
                    #run the replace looping
                    for(i in seq(1,20)){
                      zingdata$sub_plot[which(zingdata$sub_plot %in% newsubunitnamezing[,i])] <- i
                    }
                    
                zingwide <- dcast(zingdata, transecto + sub_plot ~ id_1, value.var = "quant", fun.aggregate = sum)  #species table: long to wide
                  zingwide <- zingwide[,-which(colnames(zingwide)=="na")]         # delete subplots without species

                  #change transects number to homogenized all datasets
                  transectsseq <- seq(738,777) 
                
                    #looping to change replace values  
                
                        for(i in unique(zingwide$transecto)){
                          nums <- strsplit(i, "[^[:digit:]]")   # split string non-digits
                          nums <- as.numeric(unlist(nums))      # convert string to numeric ("" become NA)
                          nums <- unique(nums[!is.na(nums)])    # remove NA and duplicates
                          zingwide[zingwide$transecto==i,"transecto"] <- transectsseq[nums]
                    }
                    
                          colnames(zingwide)[1:2] <- c("TrNumber", "subunit")       
                  
                  # write the table in .csv in the "analyses" folder
                  
                getwd()
                setwd("..") # this command takes you one folder up
                write.csv(zingwide, "zingdata_wide_gmm_v1.csv")

##########################  Moisture index for zing   #################################################

              # import species classification made by zing's specialist
              
              getwd()  # check the working directory
              setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados floristicos originais")  
              
              zingclass <- read.csv("zing_jurua_baixio.csv", stringsAsFactors = F)
                zingspnames <- colnames(zingwide[-c(1:2)])  # get the species names of the wide table
                zingclass   <- cbind(zingspnames, zingclass[-1]) # and substitute the names that came in the imported file. Names order were previously checked
                    zingclass$MI <- as.numeric(zingclass$MI)     # transform MI values in numeric
                
              # associate MI to the zingdata object
                    
                    zingdata$moistureindex <- NA # creat an empty colunm to the zingdata data.frame
                
              # associate the index to the zingdata data.frame
                    
                for(i in seq(nrow(zingclass))) {
                  zingdata[zingdata$id_1==zingclass$zingspnames[i],"moistureindex"] <- zingclass$MI[i]
                }
              
       
              # cast the zingdata 
                    
                    zingwideMI <- dcast(zingdata, transecto + sub_plot ~ id_1, value.var = "moistureindex", fun.aggregate = sum)  #species table: long to wide
                    zingwideMI <- zingwideMI[,-which(colnames(zingwideMI)=="na")]                             # delete subplots without species
                    
              # to creat the moisture index, first species relative abundances will be calculated and 
              # then weighted by the moisture index. Finally, these values will be summed by subtransect
              # in order to have a index value for each subtransect
                    
                    zingwideMIsp <- zingwideMI[,-c(1:2)]
                    zingwideMIsp[zingwideMIsp>0] <- 2     # all the values higher than 2 will be turned to 2. This is because a species can have two occurrences
                                                          # in a sinlge subtransect, making the index 4.    
                    
                    
                    relabuzing <- decostand(zingwide[,-c(1:2)], method = "total",na.rm=T )
                    moistureindextablezing <- relabuzing*zingwideMIsp
                      moistureindextablezing <- cbind(zingwide[,c(1:2)],moistureindextablezing)
                    
                    MIzing <- cbind(moistureindextablezing[,c(1:2)],rowSums(moistureindextablezing[,-c(1:2)]))
                      colnames(MIzing)[3] <- "MI"   #change the name of the moisture index column
                      
                    # save the files in the folder "Analyses"
                        
                    getwd()
                    setwd("..")
                    write.csv(MIzing, "MIzing_persubunit.csv")  
                    
      
############################# ferns  ##########################
              # organizing ferns data
              # import files      
              # the first part of this code was created by H.T.     
                    
              #load the packages
                    
              library(reshape2)
              library(vegan)
              
              # The file with the old and new species codes should have columns named "Sp_old" and "Sp_new"
              # The file with transect data should have column with name "SpCode"
              
              setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados floristicos originais")

              # change names as appropriate
              sp_data_name <- "Jurua data to analysis_ferns_from hanna.csv"  
              sp_lump_name <- "Jurua fern spp list 2013-10_03_v2.csv"
              
              # open the specified files
              sp_data      <- read.csv(sp_data_name,header = T, stringsAsFactors=FALSE)
              sp_lump_list <- read.csv(sp_lump_name,header = T, stringsAsFactors=FALSE)
              
              
              # replace old SpCodes with new ones
              for(i in seq(nrow(sp_lump_list))) {
                sp_data[sp_data$SpCode==sp_lump_list$Sp_old[i],"SpCode"] <- sp_lump_list$Sp_new[i]
              }
              
              
              # Eliminate unidentified species observations
              sp_data <- sp_data[sp_data$SpCode!="",]
              
              # replace subunits names: there were two subunits with negative values and one with 0
              unique(sp_data$subunit)
              sp_data[sp_data$subunit=="-1","subunit"] <- 19
              sp_data[sp_data$subunit=="-2","subunit"] <- 20
              
              which(sp_data$subunit==0)
              sp_data <- sp_data[-which(sp_data$subunit==0),]
              
              # associate a moisture index and the epiphyte/non-epiphyte classification
              # creat two empty columns
              
              sp_data$moistureindex <- NA
              sp_data$habit         <- NA
              
              for(i in seq(nrow(sp_lump_list))) {
                sp_data[sp_data$SpCode==sp_lump_list$Sp_old[i],"moistureindex"] <- sp_lump_list$index[i]
              }
              
              for(i in seq(nrow(sp_lump_list))) {
                sp_data[sp_data$SpCode==sp_lump_list$Sp_old[i],"habit"] <- sp_lump_list$habit[i]
              }
              
              sp_data <- sp_data[sp_data$moistureindex!="",]
              
############################### creating a moisture index for ferns #############################

              widetable <- dcast(sp_data, TrNumber + subunit ~ SpCode,value.var = "ind"           ,fun.aggregate = sum)
                      widetable <- widetable[,-which(colnames(widetable)=="NA")]
              widetableMI  <- dcast(sp_data, TrNumber + subunit ~ SpCode,value.var = "moistureindex" ,fun.aggregate = sum)
                      widetableMI <- widetableMI[,-which(colnames(widetableMI)=="NA")]
                      
                      
              #get only epiphytes
                      
                      dictionary_habits <- data.frame(dcast(sp_data, SpCode + habit ~ moistureindex)[,1:2])
                      nonepiphytes <- dictionary_habits[which(dictionary_habits$habit=="ne"),]
                      
              #select only epiphytes from the widetables
                      
                      widetableNE <-  widetable[,nonepiphytes$SpCode]
                      widetableNE <-  widetable                      #if this line is on, it does not remove the epiphytes
                      widetableNE <- cbind(widetable[1:2],widetableNE)
                              
                      
                      widetableNEMI <- widetableMI[,nonepiphytes$SpCode]
                      widetableNEMI <- cbind(widetableMI[1:2],widetableNEMI)
                            
              # transform all MI value >2 in 2
              widetableNEMIsp <- widetableNEMI[,-c(1:2)]
              widetableNEMIsp[widetableNEMIsp>0] <- 2
              widetableNEMI <- cbind(widetableNEMI[,c(1:2)], widetableNEMIsp)
              write.csv(widetableNEMI,"MIferns_persubunit_raw.csv",row.names = FALSE)
              
              #calculate the moisture index
              relabu <- decostand(widetableNE[,-c(1:2)], method = "total",na.rm=TRUE )
              moistureindextable <- relabu*widetableNEMIsp
              

              MIfernsNE <- cbind(widetableNEMI[,1:2], rowSums(moistureindextable[,-c(1:2)]))
                colnames(MIfernsNE)[3] <- "MI"  

              # save the files in the folder "Analyses"
              
                getwd()
                setwd("..")  
                
              write.csv(MIfernsNE, "MIferns_persubunit.csv")
              write.csv(widetableNE, "ferns25_widetableNE.csv",row.names = FALSE) 
   

########################### palms #############################
# organizing the palm data
# importing data

                    getwd()
                    setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados floristicos originais")
              
                    palmdata  <- read.csv( "palms_jurua_subunit_gmm1.csv")
                    palmclass <- read.csv( "Jurua_Palms_splist_habitatAff_gmm.csv")
                      palmclass$species <- as.character(palmclass$species)
                    
                    palmdatalong <- melt(palmdata, id.vars = c("TrNumber", "subunit"))
                    
                    # creat a empty column for the moiture value
                    
                    palmdatalong$moistureindex <- NA
                    
                    # associate the index to the palmdata data.frame
                    
                    for(i in seq(nrow(palmclass))) {
                      palmdatalong[palmdatalong$variable==palmclass$species[i],"moistureindex"] <- palmclass$moisture[i]
                    }
                    
                    # cast the palmdata long
                    
                    palmdataMI <- dcast(palmdatalong, TrNumber + subunit ~ variable, value.var = "moistureindex", fun.aggregate = sum)  #species table: long to wide

################### moisture index for palms ##########################################
# to creat the moisture index, first species relative abundances will be calculated and 
# then weighted by the moisture index. Finally, these values will be summed by subtransect
# in order to have a index value for each subtransect

                    palmdataMIsp <- palmdataMI[,-c(1:2)]
                    palmdataMIsp[palmdataMIsp>0] <- 2     # all the values higher than 2 will be turned to 2. This is because a species can have two occurrences
                                                          # in a sinlge subtransect, making the index 4.    
                    palmdataMIsp <- palmdataMIsp[,!is.na(colSums(palmdataMIsp))]

# calculate de relative abundance and the index

                    relabupalm <- decostand(palmdata [,-c(1:2)], method = "total",na.rm=TRUE )
                    moistureindextablepalm <- relabupalm[,colnames(palmdataMIsp)]*palmdataMIsp   # exclude de NA columns of the relaabupalm
                    moistureindextablepalm <- cbind(palmdata[,c(1:2)],moistureindextablepalm)

                    MIpalm <- cbind(moistureindextablepalm[,c(1:2)],rowSums(moistureindextablepalm[,-c(1:2)]))
                    colnames(MIpalm)[3] <- "MI"   #change the name of the moisture index column

# save files in the "Analyses" folder     
                    
                    getwd()
                    setwd("..")  
                    
                    write.csv(moistureindextablepalm,"moistureindextablepalm.csv")  # save a csv file
                    write.csv(MIpalm, "MIpalm_persubunit.csv")  # save a csv file

