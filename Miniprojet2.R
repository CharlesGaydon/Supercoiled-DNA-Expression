rm(list = ls())
setwd("~/Documents/INSA Cours/R/ExercicesStats/TDMiniprojet/miniprojet")
library("seqinr")
##Etude d'un lien statistique entre séquence promotrice et réponse à la relaxation.




acti <- function(){

  ## Importation
  anno <-read.table("annotation_dickeya.csv", header = T, sep = ",")
  anno[,1] <- paste("ECH3937_v6b_",anno[,1], sep = "")

  effets <- read.table("dAdB.csv", header = T, sep  = ",")
  effets <- subset(effets, effets[,4] == "Activation", select = c("gene"))

  seqfile <- "~/Documents/INSA Cours/R/ExercicesStats/TDMiniprojet/miniprojet/sequence.fasta"
  seq <- read.fasta(file = seqfile, seqonly = TRUE,forceDNAtolower = FALSE)
  
  ##Traitement
  dat <- data.frame()
  
  genes <- effets[,1]
  anno <- subset(anno, anno[,1]%in%genes)
  for (i in 1:length(anno[,1])){  
      if (anno[i,4] == 1){
        
        anno[i,5] <- substr(seq, anno[i,2]-300, anno[i,2]+500)
        
        anno[i,6] = strsplit(anno[i,5], "(?<=.{20})", perl = TRUE)[[1]]
        anno[i,7] = c()
        
        for (k in 1:40){
          anno[i,7][k] = prop(anno[i,6][k])
        }
        
      }else{
        
        anno[i,5] <- substr(seq, anno[i,3]-500, anno[i,3]+300)
        anno[i,5] <- rev(anno[i,5]) #A ! complémentaire
        
        anno[i,6] = strsplit(anno[i,5], "(?<=.{20})", perl = TRUE)[[1]]
        anno[i,7] = c()
        
        for (k in 1:40){
          anno[i,7][k] = 1 - prop(anno[i,6][k]) #car complémentaire
        }
        dat <- rbind(dat,anno[i,7])
        
      }
      
      #data.frame
      #dat <- cbind(anno[,1],dat)
      
      }
  }
  

#Pour chaque gène importé, trouver séquence correspondante


#la couper en 40 séquence
#pour chaque écrire la proportion de GC sur un data.frame Ceux qui ont -1 : GC = 1 - GC


datas <- acti()





prop <- function(seq){
  N <- nchar(seq, type = "chars")
  seq <- gsub(pattern = "T" , replacement = "" , seq, fixed = TRUE)
  seq <-gsub(pattern = "A" , replacement = "" , seq, fixed = TRUE)
  n <-nchar(seq, type = "chars")
  return(n/N)
}
