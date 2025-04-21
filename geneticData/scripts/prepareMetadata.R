#This script is meant to parse all the NCBI sequencing metadata into a usable format

#Load libraries
library(dplyr)

#biosample file processing -----------------------------------------------------
#Read the biosample file as text
biosampleFile <- "/scratch/dbihnam/lsc585/turtleProject/biosample_result.txt"
biosampleLines <- readLines(biosampleFile)

#Function to parse a single sample entry
parseBiosampleEntry <- function(lines) {
  #Initialize variables
  geographicLocation <- NA
  devStage <- NA
  collectionDate <- NA
  accession <- NA
  
  #Loop through the lines to extract relevant fields
  for (line in lines) {
    if (grepl("geographic location=", line)) {
      geographicLocation <- sub(".*geographic location=\"([^\"]+)\".*", "\\1", line)
    }
    else if (grepl("development stage=", line)) {
      devStage <- sub(".*development stage=\"([^\"]+)\".*", "\\1", line)
    }
    else if (grepl("collection date=", line)) {
      collectionDate <- sub(".*collection date=\"([^\"]+)\".*", "\\1", line)
    }
    else if (grepl("Accession: ", line)) {
      accession <- sub("Accession: ([^[:space:]]+).*", "\\1", line)
    }
  }
  
  #Return a data frame row
  return(data.frame(
    GeographicLocation = geographicLocation,
    DevelopmentStage = devStage,
    CollectionDate = collectionDate,
    Accession = accession,
    stringsAsFactors = FALSE
  ))
}

#Split the file into entries (assuming entries are separated by blank lines)
entries <- split(biosampleLines, cumsum(biosampleLines == ""))
entries <- entries[sapply(entries, length) > 1]  #Remove empty entries

#Parse each entry and combine into a data frame
metadataList <- lapply(entries, parseBiosampleEntry)
metadata <- do.call(rbind, metadataList)

#Check the resulting metadata
head(metadata)

#Bin specific locations into broader sites--------------------------------------
assignBin <- function(location) {
  if (is.na(location)) return("Unknown")
  
  #Convert to lowercase for easier matching
  locationLower <- tolower(location)
  
  #Netanya
  if (grepl("netanya", locationLower)) {
    return("Netanya")
  }
  #Ashkelon
  else if (grepl("ashkelon", locationLower)) {
    return("Ashkelon")
  }
  #Zikim
  else if (grepl("zikim", locationLower)) {
    return("Zikim")
  }
  #Nitzanim/South Nitzanim
  else if (grepl("nitzanim|nizanim", locationLower)) {
    return("Nitzanim")
  }
  #Gdor Nature Reserve
  else if (grepl("gdor", locationLower)) {
    return("Gdor Nature Reserve")
  }
  #Habonim Beach Nature Reserve
  else if (grepl("habonim", locationLower)) {
    return("Habonim Beach Nature Reserve")
  }
  #Hof Galim Megadim
  else if (grepl("hof galim|megadim", locationLower)) {
    return("Hof Galim Megadim")
  }
  #Ceasarea
  else if (grepl("ceasarea", locationLower)) {
    return("Ceasarea")
  }
  #Yavne (including Palmahim, Nahal Poleg, Ashdod)
  else if (grepl("yavne|palmahim|nahal poleg|ashdod", locationLower)) {
    return("Yavne")
  }
  #Nahariyya/Betzet
  else if (grepl("nahariyya|betzet", locationLower)) {
    return("Nahariyya-Betzet")
  }
  #Captive Breeding
  else if (grepl("captive breeding", locationLower)) {
    return("Captive Breeding")
  }
  #Stranded
  else if (grepl("stranded", locationLower)) {
    return("Stranded")
  }
  #Central-Northern Coast (Other)
  else if (grepl("giv'at olga|atlit|alexander river|reading port|north sharon|neve yam", locationLower)) {
    return("Central-Northern Coast (Other)")
  }
  else {
    return("Unknown")
  }
}

#Add the binned geographic region to the metadata
metadata$BinnedLocation <- sapply(metadata$GeographicLocation, assignBin)

#Check the distribution of binned locations
table(metadata$BinnedLocation)

#Swap SAMN accession# for SRR sequencing ID-------------------------------------
#Load the srr_biosample.txt file
srrBiosample <- read.table("/scratch/dbihnam/lsc585/turtleProject/srr_biosample.txt", header = FALSE, stringsAsFactors = FALSE)

#Assign column names
colnames(srrBiosample) <- c("SequencingSampleID", "Accession")

#Check the first few rows
head(srrBiosample)

#Merge the metadata with srrBiosample based on Accession
metadataUpdated <- merge(metadata, srrBiosample, by = "Accession", all.x = TRUE)

#Drop NA value (extra sample that was not sequenced)
metadataUpdated <- na.omit(metadataUpdated)

#Clean up
metadataFin <- metadataUpdated %>% 
  select(SequencingSampleID, BinnedLocation, DevelopmentStage, CollectionDate)

#Save finalized metadata table--------------------------------------------------
write.csv(metadataFin, "/scratch/dbihnam/lsc585/turtleProject/sequencingMetadata.csv", row.names = FALSE)
