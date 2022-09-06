library(optparse)

# read user info
option_list = list(
  make_option(c("-i", "--mafffile"), type="character", default=NULL, help="input maffile", metavar="character"),
  make_option(c("-o", "--tmbfile"), type="character", default=NULL, help="output TMB file to be created", metavar="character"),
  make_option(c("-p", "--proteinAltering"), action="store_true", default=FALSE, help="subset to protein altering mutations"), 
  make_option(c("-c", "--callableSpace"), type="numeric", default=1, help="callable space for the target/capure in Mb", metavar="numeric") 
)

# get options
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


# calculate TMB
library(dplyr)

subsetProteinAlteringMutations <- function(data.mutations){
  include <- c("Missense_Mutation", 
               "In_Frame_Ins", 
               "In_Frame_Del", 
               "Frame_Shift_Ins",
               "Frame_Shift_Del",
               "Splice_Site",
               "Translation_Start_Site",
               "Nonsense_Mutation",
               "Nonstop_Mutation",
               "Silent")
  
  data.mutations <- data.mutations[data.mutations$Variant_Classification %in% include,]
  return (data.mutations)
}

tmbCalc <- function(muts.file, out.file, Callable_space, subsetPA = TRUE){
  data.mutations <- read.csv(muts.file,
                             sep = "\t", 
                             as.is = T)
  if (subsetPA){
    data.mutations <- subsetProteinAlteringMutations(data.mutations)
  }
  
  if (nrow(data.mutations) == 0) {
    TMB <- c()
    TMB[["sample"]] <- c(0,0,Callable_space)
  }
  else {
    sample.names <- unique(data.mutations$Tumor_Sample_Barcode)
    TMB <- c()
    
    # Callable_space <- 37.2855
    for (s in sample.names){
       sample.mutations <- data.mutations[data.mutations$Tumor_Sample_Barcode == s,]
       Total_Mutations <- dim(sample.mutations)[1]
       Mutation_burden <- Total_Mutations/Callable_space
       TMB[[s]] <- c(Total_Mutations, Mutation_burden, Callable_space)
    }
  } 
  TMB <- data.frame(t(data.frame(TMB)))
  colnames(TMB) <- c("Total_Mutations", "Mutation_burden", "Callable_space")
  TMB$Sample_ID <- gsub("X" , "", row.names(TMB))
  TMB <- TMB[,c("Sample_ID", "Total_Mutations", "Mutation_burden", "Callable_space")]
  # write to file
  write.table(TMB, file = out.file, sep = "\t", row.names = F, quote = F)
  return (TMB)
}


# # SCOUT test
# muts.file.scout <- "/Volumes/TGL/gsi/pipeline/data/TGL25/cbioportal/cBioWrap_20191107/output/cbioportal_import_data/data_mutations_extended.txt"
# out.file <- "/Volumes/TGL/gsi/pipeline/data/TGL25/cbioportal/cBioWrap_20191107/output/supplementary_data/tmb.txt"
# Callable_space <- 37.2855
# tmb.scout <- tmbCalc(muts.file.scout, out.file, Callable_space)

tmbCalc(muts.file = opt$mafffile, out.file = opt$tmbfile, Callable_space = opt$callableSpace, subsetPA = opt$proteinAltering)

