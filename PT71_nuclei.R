#### Corelation of IRF1 and pSTAT1 ####
## microscopic data ##
source("/Users/piotrt/Documents/IPPT_PT/R/basic_scripts/normalize_data.R")
source("/Users/piotrt/Documents/IPPT_PT/R/basic_scripts/theme_jetka.R")
#### data preparation- loading and normalization ####
experiment <- "PT71"
path <- paste("/Users/piotrt/Documents/IPPT_PT/ICF/Eksperymenty/",
              experiment, "/R/", sep='')
normaliz <- "input/raw/"
reading <- "second"

setwd(path)
load.bio.nuclei <- function(dye, path, normaliz,
                               front.columns = 5, back.columns = 12){
  
  df.nuc <- read.table(paste(path, normaliz, "ShrinkedNucleiMasked", 
                             '.csv', sep=''), 
                       header=TRUE, sep=",")

  df.nuc <- normalize_data(df.nuc)$data
  
  columns <- c("Intensity_IntegratedIntensity_Alexa488",
               "Intensity_MeanIntensity_Alexa488")
  
  new.column.core <- c("Integrated_Alexa488",
                       "Mean_Alexa488")
  
  df.nuc2 <- variable.subset(df.nuc, columns, 
                             new.columns = (new.column.core))
  
  return(cbind(df.nuc[, head(colnames(df.nuc), front.columns)],
               df.nuc2,
               df.nuc[, tail(colnames(df.nuc), back.columns)]))
}

bio.df.nuclei <- load.bio.nuclei(path=path,
                           normaliz = paste(normaliz, 
                                            reading, "/", sep = ""))
head(bio.df.nuclei)

# pomysł: dodać 4 kolumny: najbliższy, 2nd, 3rd i 4th. Leciec po dołkach, 
# w każdej komórce lecieć po wierszach i dla każdego wiersza wykonać funkcję
# która znajduje najbliższych 4 sąsiadów. Prawdopodobnie obliczenia zajmą 
# dużo czasu, lepiej na PC.

euc.dist <- function(x1, x2) {sqrt(sum((x1 - x2) ^ 2))}
append.list <- function(list, appendix){ 
  list[[length(list)+1]] <- appendix
  return(list)
  }

all.characteristic <- data.frame()

for(well in unique(bio.df.nuclei$well.name)){
  
# bio.well <- bio.df.nuclei[bio.df.nuclei$well.name ==
#                             unique(bio.df.nuclei$well.name)[1], ][1:10, ]
  bio.well <- bio.df.nuclei[bio.df.nuclei$well.name == well, ]
  
well.nuclei <- length(bio.well[, 1])

list.of.used <- list(best = list(), closest = list())
closests.df <-  setNames(data.frame(matrix(ncol = 8, nrow = well.nuclei)),
                         c("distance.best", 
                           "no.best", 
                           "Mean_Alexa488.best", 
                           "benchamrk_difference.best",
                           "distance.closest", 
                           "no.closest", 
                           "Mean_Alexa488.closest", 
                           "benchamrk_difference.closest"))

for(nuclei.start in (1:well.nuclei)){
  coordinates.start <- bio.well[nuclei.start, c("AreaShape_Center_X",
                                             "AreaShape_Center_Y")]
  relation.df <-  setNames(data.frame(matrix(ncol = 2, nrow = well.nuclei)),
                                      c("distance", "no"))
  
  fluo.benchmark <- bio.well[nuclei.start, ]$Mean_Alexa488
  for(nuclei.neighbour in (1:well.nuclei)){
    
    coordinates.neighbour <- bio.well[nuclei.neighbour, 
                                      c("AreaShape_Center_X",
                                        "AreaShape_Center_Y")]
    distance <- euc.dist(coordinates.start,
                         coordinates.neighbour)
    
    relation.df[nuclei.neighbour, ] <- c(distance,
                                         nuclei.neighbour)
  }
    closest.cell <- relation.df[order(relation.df$distance), ][2:3, ]
    closest.cell <- closest.cell[closest.cell$distance<=1000, ]
    if(nrow(closest.cell)==0) {
      nucleus.characteristic.best <- rep(0, 4)
      nucleus.characteristic.close <- rep(0, 4)
    } else {
    closest.cell$Mean_Alexa488 <- bio.well[closest.cell$no, ]$Mean_Alexa488
    closest.cell$benchamrk_difference <- 
      abs(closest.cell$Mean_Alexa488 - fluo.benchmark)
    
    best.neighbout <- 
      which(closest.cell$benchamrk_difference==
              min(closest.cell$benchamrk_difference))
    if(closest.cell[best.neighbout, 2] %in% list.of.used[["best"]] |
       nuclei.start %in% list.of.used[["best"]]) { 
      nucleus.characteristic.best <- rep(0, 4)
    } else {
      nucleus.characteristic.best <- closest.cell[best.neighbout, ]
      list.of.used[["best"]] <- append.list(list.of.used[["best"]], 
                                            closest.cell[best.neighbout, 2])
    }
    if(closest.cell[1, 2] %in% list.of.used[["closest"]] |
       nuclei.start %in% list.of.used[["best"]]){
      nucleus.characteristic.close <- rep(0, 4)
    } else {
      nucleus.characteristic.close <- closest.cell[1, ]
      list.of.used[["closest"]] <- append.list(list.of.used[["closest"]], 
                                               closest.cell[1, 2])
    }
    
    nucleus.characteristic <- c(nucleus.characteristic.best,
                                nucleus.characteristic.close)
    
    closests.df[nuclei.start, ] <- nucleus.characteristic
    }
}
all.characteristic <- rbind(all.characteristic, closests.df)
}
cbind(bio.df.nuclei, all.characteristic)

