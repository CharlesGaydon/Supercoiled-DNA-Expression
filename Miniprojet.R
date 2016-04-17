rm(list = ls())
setwd("~/Documents/INSA Cours/R/ExercicesStats/TDMiniprojet/miniprojet")

#Note : on a traité séparément les deux réplicats, mais après correction (cf: q7) on aurait pu rassembler les données.


#1 et 2
#Chaque sonde étant différente, on ne peut comparer directement que les données issues d'une même sonde.
#L'effet de chaque sonde étant le même dans les expériences témoin et test, en soustrayant à la fluorescence test
#la fluorescence témoin obtenue avec une même sonde et pour un même gène, on obtient des série de valeurs comparables
#pour chacun des gènes, qu'on compare à 0 pour conclure à un effet.
A1 <- read.table("sans_1.pair", header = T, skip = 2)
A2 <- read.table("relax_1.pair", header = T, skip = 2)
B1 <- read.table("sans_2.pair", header = T, skip = 2)
B2 <- read.table("relax_2.pair", header = T, skip = 2)

dA1 <- subset(A1, A1[,2] != "RANDOM")
dB1 <- subset(B1, B1[,2] != "RANDOM")
dA2 <- subset(A2, A2[,2] != "RANDOM")
dB2 <- subset(B2, B2[,2] != "RANDOM")

dA <- data.frame(gene = dA1[,3],diff = log(dA2[,9])-log(dA1[,9]))
dB <- data.frame(gene = dB1[,3],diff = log(dB2[,9])-log(dB1[,9]))

dAm <- aggregate(dA$diff, by=list(dA$gene), FUN=mean) #gene - moyenne de difference des log
sd <- aggregate(dA$diff, by=list(dA$gene), FUN=sd)
names(dAm) <- c("gene","diff")
dAm$sd <- sd[,2] #gene - moy - sd

dBm <- aggregate(dB$diff, by=list(dB$gene), FUN=mean) 
sd <- aggregate(dB$diff, by=list(dB$gene), FUN=sd)
names(dBm) <- c("gene","diff")
dBm$sd <- sd[,2] #gene - moy - sd

#3
#Comparaison à l'aide d'un test de Student de la moyenne des différences pour chaque gène, qui suit une loi 
#normale de moyenne 0 sous H0 : "La novobiocine n'a pas d'effet sur l'expression du gène".

#### Traitement de A ####

dAm$pvalue <- 2*(1 - pt(abs(dAm$diff)/(dAm$sd/sqrt(15)), df = 14))
dAm <- dAm[order(dAm$pvalue,decreasing = F),]
n <- length(dAm$gene)*0.1
dAm <- dAm[1:n,]
dAm$effet <- ifelse(dAm$diff>0, "Activation", "Repression")

fA <- dAm$pvalue[length(dAm$pvalue)] #pvalue la moins significative : 5.367262e-12
#summary(dAm$diff>0) 432 activations, 23 répressions.
FWERA <- 1 - (1-fA)**length(dAm$diff) #2.442e-09

dAm <- subset(dAm, select = c("gene","pvalue","effet"))
write.csv(dAm,"dA.csv", quote=FALSE) 



#### Traitement de B ####

dBm$pvalue <- 2*(1 - pt(abs(dBm$diff)/(dBm$sd/sqrt(15)), df = 14))
dBm <- dBm[order(dBm$pvalue,decreasing = F),]
n <- length(dBm$gene)*0.1
dBm <- dBm[1:n,]
dBm$effet <- ifelse(dBm$diff>0, "Activation", "Repression")

#summary(dBm$diff>0) #411 répression, 44 activations.
fB <- dBm$pvalue[length(dBm$pvalue)] #pvalue la moins significative : 1.846285e-09
FWERB <- 1 - (1-fB)**length(dBm$diff) #8.400e-07

dBm <- subset(dBm, select = c("gene","pvalue","effet")) 
write.csv(dBm,"dB.csv", quote=FALSE)

#Risque de faire une erreur sur au moins un gène très faible même en majorant!

#Influence des puces est-elle négligeable ? Comparons les résultats des tests standard : on a généré un certain nombre de séquences
#aléatoires, identiques d'une puce à l'autre, qui ne s'hybrideront à aucun ARN exprimé depuis l'organisme étudié. Les valeurs
#de fluorescence obtenues traduisent donc un bruit de fond qui peut varier d'une puce à l'autre. On effectue donc le même type de test
#que précédemment, en étudiant si la différence des log des fluorescences est significativement différente de zero ou non.

dA1 <- subset(A1, A1[,2] == "RANDOM")
dB1 <- subset(B1, B1[,2] == "RANDOM")
dA2 <- subset(A2, A2[,2] == "RANDOM")
dB2 <- subset(B2, B2[,2] == "RANDOM")

dA <- data.frame(gene = dA1[,3],diff = log(dA2[,9])-log(dA1[,9]))
dB <- data.frame(gene = dB1[,3],diff = log(dB2[,9])-log(dB1[,9]))
n <- length(dA[,1])


dAm <- aggregate(dA$diff, by=list(dA$gene), FUN=mean) #gene - moyenne de difference des log
sd <- aggregate(dA$diff, by=list(dA$gene), FUN=sd)
names(dAm) <- c("gene","diff")
dAm$sd <- sd[,2] #gene - moy - sd

dBm <- aggregate(dB$diff, by=list(dB$gene), FUN=mean) 
sd <- aggregate(dB$diff, by=list(dB$gene), FUN=sd)
names(dBm) <- c("gene","diff")
dBm$sd <- sd[,2] #gene - moy - sd

#Comparaison à l'aide d'un test de Student de la moyenne des différences. Sous l'hypothèse nulle H0 : "bpuse = 0" la moyenne des écarts divisée par
#sigma/sqrt(n) suit une loi de Student à n-1 degré de liberté

pvalueA <- 2*(1 - pt(abs(dAm$diff)/(dAm$sd/sqrt(n)), df = n-1))
pvalueB <- 2*(1 - pt(abs(dBm$diff)/(dBm$sd/sqrt(n)), df = n-1))

#Nulles toute deux : on rejette l'hypothèse nulle de manière certaine (hehe).
#cette différence des moyennes des différences de log peut avoir plusieurs origines : 
#une plus grande sensibilité d'une puce à la fluorescence (facteur mulitplicatif appliqué à f), 
#une plus grande fluorescence dirrectement (puce "plus accessible"?),
#l'ajout d'un bruit de fond. 
#Pour simplifier on prend cette dernière hypothèse pour juste et on va corriger les valeurs précédemment obtenue
#de cette moyenne du bruit de fond, tout en rassemblant les valeurs ainsi corrigées.

########### A Nouveau #############

dA1 <- subset(A1, A1[,2] != "RANDOM")
dB1 <- subset(B1, B1[,2] != "RANDOM")
dA2 <- subset(A2, A2[,2] != "RANDOM")
dB2 <- subset(B2, B2[,2] != "RANDOM")

dA <- data.frame(gene = dA1[,3],diff = log(dA2[,9])-log(dA1[,9])-dAm[1,2])
dB <- data.frame(gene = dB1[,3],diff = log(dB2[,9])-log(dB1[,9])-dBm[1,2])

dAB <- rbind(dA,dB)

dAmBm <- aggregate(dAB$diff, by=list(dAB$gene), FUN=mean) #gene - moyenne de difference des log
sd <- aggregate(dAB$diff, by=list(dAB$gene), FUN=sd)
names(dAmBm) <- c("gene","diff")
dAmBm$sd <- sd[,2] #gene - moy - sd



#3
#Comparaison à l'aide d'un test de Student de la moyenne des différences pour chaque gène, qui suit une loi 
#normale de moyenne 0 sous H0 : "La novobiocine n'a pas d'effet sur l'expression du gène".

#### Traitement de A et B ####

dAmBm$pvalue <- 2*(1 - pt(abs(dAmBm$diff)/(dAmBm$sd/sqrt(30)), df = 29))
dAmBm <- dAmBm[order(dAmBm$pvalue,decreasing = F),]
n <- length(dAmBm$gene)*0.1
dAmBm <- dAmBm[1:n,]
dAmBm$effet <- ifelse(dAmBm$diff>0, "Activation", "Repression")

fA <- dAmBm$pvalue[length(dAmBm$pvalue)] #pvalue la moins significative : 0
summary(dAmBm$diff>0) #430 activations, 25 répressions
FWERA <- 1 - (1-fA)**length(dAmBm$diff) # 0 : on gagne en précision en aggrégant les deux données et en les corrigeant!

dAmBm <- subset(dAmBm, select = c("gene","pvalue","effet"))
write.csv(dAmBm,"dAdB.csv", quote=FALSE) #Fichier après correction et rassemblement


