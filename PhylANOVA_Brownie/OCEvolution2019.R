########
#Coded by Cristina Robinson
#Last Modified 5-1-2019
########

#SET UP!!!
rm(list = objects())
#setwd()##set your directory here
dir <- getwd()
source("Source_OCEvolution.R")#Libraries and scripts

#load trees and data (this takes a while)
ConsensusTree <- LoadPrettyTree("HacketTrees.nex") #note that the variable Path is also generated here
ConsensusTree$tip.label[which(ConsensusTree$tip.label == "Philesturnus_carunculatus")] <- "Philesturnus_rufusater"
birbs <- Loadbirds("MainDataset.csv", ConsensusTree)
OCVariants <- LoadOtherOC("StabilityDataDiTriCont.csv", birbs)
remove <- which(is.na(OCVariants$tri))
VariantTree <- drop.tip(ConsensusTree, remove)

#Set variables and filestructure to output data
OpenClose <- birbs$O.C
OpenClose <- factor(OpenClose, labels = c("closed","open"))
names(OpenClose) <- birbs$BirdtreeFormat

OpenCloseTri <- OCVariants$tri[-remove]
OpenCloseTri <- factor(OpenCloseTri, labels = c("closed","delayed-closed","open"))
names(OpenCloseTri) <- OCVariants$BirdtreeFormat[-remove]
#Loop variables for main analysis:
call <- names(birbs)[1:(length(birbs))]
call <- call[!call %in% c("O.C", "BirdtreeFormat", "Family", "FormerFamily")]#remove non-independant vars
mod <- rep("log", length(call))
mod[c(1,19)] <- "linear"
FullRate <- RateMatx(ace(x=OpenClose, phy=ConsensusTree, type="discrete", model="ER"))
FullRateTri <- RateMatx(ace(x=OpenCloseTri, phy=VariantTree, type="discrete", model="ER"))
#Loop variables for the Jacknife; based on significnace from the main analysis
mod2 <- mod[c(2, 5, 8, 11, 14, 18)]
call2 <- call[c(2, 5, 8, 11, 14, 18)]

#Create folder Structure
MakeFolderStructure(dir)








#Parsimony Analyis
#Get the number of transitions.  This takes  while!
ParseTrees <- make.simmap(ConsensusTree, OpenClose, "ER", 10000)
OCRoots <- table(unlist(sapply(1:10000, function(x) names(ParseTrees[[x]]$maps[[1]][1]))))

save <- countSimmap(ParseTrees)
colMeans(save$Tr)
min(save$Tr[,1])









####FIGURES and TABLE DATA:

#Figure 1 V2: make a histogram and boxplot
pdf("Figure 1.pdf")
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
QuickScatterBox(vari=cbind(birbs$Syllable.rep.final,birbs$Song.rep.final),
                OC=OpenClose,title=c("Syllable","Song"), DOIT=FALSE)
QuickScatterBox(vari=cbind(birbs$Syllable.rep.final[-remove]),
                OC=OpenCloseTri,
                labels=c("Song-Stable", "longer Learning", "Song-Plastic"),
                title=c("Syllable"), DOIT=FALSE)
plot(OCVariants$cont[-remove],birbs$Syllable.rep.final[-remove], col=rgb(1,0,1),
     xlab="Years Spent Learning", ylab="Syllable Repertoire",
     log='y', font.lab=2)
linmodel <- lm(x~y, list(x=OCVariants$cont[-remove], y=log(birbs$Syllable.rep.final[-remove])))
summary(linmodel)
abline(linmodel,
       lwd=2)
dev.off()


#Figure 2 V2:
pdf("Figure 2.pdf")
par(mfrow=c(1,2), mar=c(.1,.1,.1,.1))
DataExtraction(OpenClose, birbs, ConsensusTree, call[2], mod=mod[2], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE,Flip="rightwards")
DataExtraction(OpenCloseTri, birbs[-remove,], VariantTree, call[2], mod=mod[2], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE,Flip="leftwards")
dev.off()



#Figure 2AB
#creates RainbowPlots, runs phylANOVA, outputs text, runs Brownielite, plots it
pdf("DoublePlot.pdf")
par(mfrow=c(1,2), mar=c(.1,.1,.1,.1))
DataExtraction(OpenClose, birbs, ConsensusTree, call[2], mod=mod[2], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE,Flip="rightwards")
DataExtraction(OpenClose, birbs, ConsensusTree, call[8], mod=mod[8], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE,Flip="leftwards")
dev.off()






#Table 1 and 2 data
setwd(file.path(dir, "DataWarehouse"))
ANOVAData <- as.list(1:(length(call)-2))
for(i in 18:2){
#get data and rainbow plots
  DataExtraction(OpenClose, birbs, ConsensusTree, call[i], mod=mod[i], fullrate=FullRate)
  ANOVAData[[i-1]] <- ANOVARun #created by DataExtraction()
#plotting brownie data
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], ".Brownie.pdf"))
    par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
    sink(file = "Brownie.txt", append = TRUE, split = FALSE)
      BrowniePlotRates(dataset, paste0(call[i]), Group=c("Stable", "Plastic"))
      BrowniePlotRates(datasetFULL, paste0(call[i],"FULL"), Group=c("Stable", "Plastic"))
    sink(file=NULL)
  dev.off()
}
#get Anova results
ANOVAResults(ANOVAData)

#Figure 3, 7 in one brownie :)
pdf("BrownieFullRatesAlltypes.pdf")
par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
for(i in c(2,8,5,11,14,17,18)){#FUll rates
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
                 sep="", collapse=" ")
  title <- gsub("Final", "", title)
  NAind <- which(is.na(datasetFULL$ARDRate1))
  if(length(NAind)>0){
    MAX <- datasetFULL$ARDRate1[-NAind]
  }else{ MAX <- max(datasetFULL$ARDRate1)}
  BrowniePlotRates(datasetFULL, title, Groups=c("Stable", "Plastic"),
                   Xlim=c(0, max(.2, MAX)))
}
dev.off()














#tristates
setwd(file.path(dir, "DataWarehouse/Tri"))
QuickScatterBox(vari=cbind(birbs$Syllable.rep.final[-remove], birbs$Song.rep.final[-remove]),
                OC=OpenCloseTri,
                labels=c("Song-Stable", "longer Learning", "Song-Plastic"),
                title=c("Syllable","Song"))

pdf("DoublePlotTri.pdf")
par(mfrow=c(1,2), mar=c(.1,.1,.1,.1))
DataExtraction(OpenCloseTri, birbs[-remove,], VariantTree, call[2], mod=mod[2], fullrate=FullRateTri,
               BROWNIE = FALSE, ANOVA = TRUE,DP=TRUE,Flip="rightwards")
DataExtraction(OpenCloseTri, birbs[-remove,], VariantTree, call[8], mod=mod[8], fullrate=FullRateTri,
               BROWNIE = FALSE, ANOVA = TRUE,DP=TRUE,Flip="leftwards")
dev.off()


ANOVAData <- as.list(1:(length(call)-2))
for(i in 18:2){
  #get data and rainbow plots
  DataExtraction(OpenCloseTri, birbs[-remove,], VariantTree, call[i], mod=mod[i], fullrate=FullRateTri)
  ANOVAData[[i-1]] <- ANOVARun #created by DataExtraction()
  #plotting brownie data
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], "Tri.Brownie.pdf"))
  par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
  sink(file = "TriBrownie.txt", append = TRUE, split = FALSE)
  BrowniePlotRates(dataset, paste0(call[i], "Tri"),
                   col = c('blue', 'purple', 'red'),
                   Group=c("Stable", "Longer-Learning", "Plastic"))
  BrowniePlotRates(datasetFULL, paste0(call[i],"FULL-Tri"),
                   col = c('blue', 'purple', 'red'),
                   Group=c("Stable", "Longer-Learning", "Plastic"))
  sink(file=NULL)
  dev.off()
}



#get Anova results
ANOVAResults(ANOVAData)

#Figure 3, 7 in one brownie :)
pdf("BrownieFullRatesAlltypes.pdf")
par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
for(i in c(2,8,5,11,14,17,18)){#FUll rates
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
        sep="", collapse=" ")
  title <- gsub("Final", "", title)
  NAind <- which(is.na(datasetFULL$ARDRate1))
  if(length(NAind)>0){
    MAX <- datasetFULL$ARDRate1[-NAind]
  }else{ MAX <- max(c(datasetFULL$ARDRate1, datasetFULL$ARDRate0, datasetFULL$ARDRate2))}
  BrowniePlotRates(datasetFULL, title,
                   Xlim=c(0, max(.2, MAX)),
                   col=c('blue', 'purple', 'red'),
                   Group=c("Stable", "Delayed", "Plastic"))
}
dev.off()

MakeAllNodePlots(birbs, call, OCVariants, ConsensusTree)




#compare 2 to 3 rates
setwd(file.path(dir, "DataWarehouse/Tri/Di"))
for(i in c(2,5,8)){
  DataExtraction(OpenClose[-remove], birbs[-remove,], VariantTree, call[i], mod=mod[i], fullrate=FullRate,
                 RAIN = FALSE, ANOVA = FALSE)  
}

setwd(file.path(dir, "DataWarehouse"))
sink(file = "BrownieTriDi.txt", append = TRUE, split = FALSE)
for(i in c(2,5,8)){
  
  datasetDi <- read.csv(paste0("Tri/Di/",call[i],"FULL.csv"))
  datasetTri <- read.csv(paste0("Tri/",call[i],"FULL.csv"))
  Mean2 <- mean(datasetDi$ARDloglik)
  Mean3 <- mean(datasetTri$ARDloglik)
  pval <- round(pchisq(2*(Mean3 - Mean2),1,lower.tail=FALSE),digits=3)
  ifelse(pval == 0,pval <- "<0.001", pval <- paste0("=",pval))  
  writeLines(call[i])
  writeLines(paste0("TwoRate=", round(Mean2, digits = 4)))
  writeLines(paste0("ThreeRates=", round(Mean3, digits = 4)))
  writeLines(paste0("pVal", pval))
  writeLines(paste("",sep="\n\n"))
  writeLines(paste("",sep="\n\n"))
}
sink(NULL)

setwd(file.path(dir, "DataWarehouse/Tri/newDi"))
OpenCloseTriSwitch <- OpenCloseTri
OpenCloseTriSwitch[which(OpenCloseTriSwitch=="delayed-closed")] <- "open"
OpenCloseTriSwitch <- droplevels(OpenCloseTriSwitch)
DataExtraction(OpenCloseTriSwitch, birbs[-remove,], VariantTree, call[2], mod=mod[2], fullrate=FullRate,
               RAIN = FALSE, ANOVA = FALSE)  

sink(file = "BrownieTrinewDi.txt", append = TRUE, split = FALSE)
for(i in c(2,8)){
  
  datasetDi <- read.csv(paste0("Tri/newDi/",call[i],"FULL.csv"))
  datasetTri <- read.csv(paste0("Tri/",call[i],"FULL.csv"))
  Mean2 <- mean(datasetDi$ARDloglik)
  Mean3 <- mean(datasetTri$ARDloglik)
  pval <- round(pchisq(2*(Mean3 - Mean2),1,lower.tail=FALSE),digits=3)
  ifelse(pval == 0,pval <- "<0.001", pval <- paste0("=",pval))  
  writeLines(call[i])
  writeLines(paste0("TwoRate=", round(Mean2, digits = 4)))
  writeLines(paste0("ThreeRates=", round(Mean3, digits = 4)))
  writeLines(paste0("pVal", pval))
  writeLines(paste("",sep="\n\n"))
  writeLines(paste("",sep="\n\n"))
}
sink(NULL)


setwd(file.path(dir, "DataWarehouse"))
#Figure 3:
pdf("Figure 4.pdf")
par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
for(i in c(2,8,5)){#FUll rates
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
                 sep="", collapse=" ")
  title <- gsub("Final", "", title)
  NAind <- which(is.na(datasetFULL$ARDRate1))
  if(length(NAind)>0){
    MAX <- datasetFULL$ARDRate1[-NAind]
  }else{ MAX <- max(datasetFULL$ARDRate1)}
  BrowniePlotRates(datasetFULL, title, Groups=c("Stable", "Plastic"),
                   Xlim=c(0, max(.2, MAX)))
  
  datasetFULL <- read.csv(paste0("Tri/",call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
                 sep="", collapse=" ")
  title <- gsub("Final", "", title)
  NAind <- which(is.na(datasetFULL$ARDRate1))
  if(length(NAind)>0){
    MAX <- datasetFULL$ARDRate1[-NAind]
  }else{ MAX <- max(c(datasetFULL$ARDRate1, datasetFULL$ARDRate0, datasetFULL$ARDRate2))}
  BrowniePlotRates(datasetFULL, title,
                   Xlim=c(0, max(.2, MAX)),
                   col=c('blue', 'purple', 'red'),
                   Group=c("Stable", "Delayed", "Plastic"))    

  if(i != 5){
    sink(file = "BrownienewDi.txt", append = TRUE, split = FALSE)
    datasetFULL <- read.csv(paste0("Tri/newDi/", call[i],"FULL.csv"))
    title <- unlist(strsplit(call[i],"[.]"))
    title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
                   sep="", collapse=" ")
    title <- gsub("Final", "", title)
    NAind <- which(is.na(datasetFULL$ARDRate1))
    if(length(NAind)>0){
      MAX <- datasetFULL$ARDRate1[-NAind]
    }else{ MAX <- max(datasetFULL$ARDRate1)}
    BrowniePlotRates(datasetFULL, title, Groups=c("Shorter", "Longer"),
                     Xlim=c(0, max(.2, MAX)))
    sink(NULL)
  }
}
dev.off()


pdf("Figure 5.pdf")
par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
for(i in c(11,14,17,18)){#FUll rates
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
                 sep="", collapse=" ")
  title <- gsub("Final", "", title)
  NAind <- which(is.na(datasetFULL$ARDRate1))
  if(length(NAind)>0){
    MAX <- datasetFULL$ARDRate1[-NAind]
  }else{ MAX <- max(datasetFULL$ARDRate1)}
  BrowniePlotRates(datasetFULL, title, Groups=c("Stable", "Plastic"),
                   Xlim=c(0, max(.2, MAX)))
  
}
dev.off()



#3)Jackknife runs, generates trees we did not show and the table data
#Because Acrocephalidae is paraphyletic, we merged the Lucustellidae with Acrocephalidae
birbs$Family[which(birbs$BirdtreeFormat == "Locustella_naevia")] <- "Acrocephalidae"
#create a list of indicies belonging to each familiy, get indicied and figure out which have 4+ species
Families <- replicate(length(levels(birbs$Family)),NULL)
names(Families) <- levels(birbs$Family)
for(i in 1:length(Families)){Families[[i]]<-which(birbs$Family == names(Families)[i])}
Remove <- names(which(sapply(Families,length)>=4))
Type <- c("", "FULL")


#first loop (i) enters folder for song variable and sets up ANOVA data
#second loop (j) cuts out each of the families in turn and runs the dataextraction protocol
#Third loop (k) run loop 4 with full and partial rates
#Fourth loop (l) generate and plot brownie data
#loop 1
for(i in seq_along(call2)){
  setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
  ANOVAData <- as.list(1:length(Remove))
#loop 2: Jackknife using the ACE values from the tree created after species with NAs for
#a song variable were removed and those which were removed by the jacknife procedure itself
  for(j in 1:length(Remove)){
    ConseJack <- drop.tip(ConsensusTree, Families[[Remove[j]]], root.edge = 0)
    Jacks <- birbs[-Families[[Remove[j]]],]
    OC <- OpenClose[-Families[[Remove[j]]]]
    DataExtraction(OC, Jacks, ConseJack, vari=call2[i], RAIN=FALSE, 
                   mod=mod2[i], cotitle=paste0("No", Remove[j]), fullrate=FullRate)
    ANOVAData[[j]] <- ANOVARun
  }
#loop 3: runs full and partial rates
  sink(file = "Brownie.txt", append = TRUE, split = FALSE)
  for(k in 1:2){
    pdf(paste0(call2[i], Type[k],"Jackknife.pdf"))
    par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
#Plot the original brownie run
    setwd(file.path(dir, "DataWarehouse"))
    dataset <- read.csv(paste0(call2[i],Type[k],".csv"))
    BrowniePlotRates(dataset,paste0(call2[i], " All"))
    setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
#loop 4 Brownie plot the jacknife Runs   
    for(l in 1:length(Remove)){
      dataset <- read.csv(paste0(call2[i],"No",Remove[l],Type[k],".csv"))
      JackLoss <- length(which(is.na(birbs[,call2[i]][Families[[Remove[l]]]])==FALSE))
      BrowniePlotRates(dataset,paste0(call2[i]," No ",Remove[l], "(", JackLoss, ")"))  
    }
    dev.off()
  }
  sink(file = NULL)
  ANOVAPrinter(ANOVAData, CritAlpha[which(call == call2[i])-1])
}




#4) Jacknife with individual Mimids:
#based on the data from above, we decided to repeat the Brownie Analysis
#with each Mimid removed in turn for syll.song
setwd(file.path(dir, "DataWarehouse/MimidJackknife"))
speciesIndex <- Families$Mimidae
NAind <- which(is.na(birbs$Syll.song.final)==TRUE)
SylSong <- birbs$Syll.song.final
names(SylSong) <- birbs$BirdtreeFormat
for(i in 1:length(speciesIndex)){
  ConseMime <- drop.tip(ConsensusTree, c(speciesIndex[i],NAind), root.edge = 0)
  Mime <- birbs[-c(speciesIndex[i], NAind),]
  OCmime <- OpenClose[-c(speciesIndex[i],NAind)]
  sySo <- log(SylSong[-c(speciesIndex[i],NAind)])
  BrownieDataGen(ConseMime, OCmime, sySo, nsim=1300,title=paste(birbs$BirdtreeFormat[speciesIndex[i]]), FullRate)
}

pdf("MimidJackkinfe.pdf")
par(mfrow=c(2,2))
sink(file = "Brownie.txt", append = TRUE, split = FALSE)
for(i in 1:length(speciesIndex)){
  dataset <- read.csv(paste0(birbs$BirdtreeFormat[speciesIndex[i]], ".csv"))
  BrowniePlotRates(dataset, paste("Syl.Song","No",birbs$BirdtreeFormat[speciesIndex[i]], sep=" "))
}
sink(file = NULL)
dev.off()

 




#6) Test with Lincolnii Closed
setwd(file.path(dir, "DataWarehouse/ClosedLink"))
OClink <- OpenClose
OClink[which(birbs$BirdtreeFormat == "Melospiza_lincolnii")] <- "closed"
#creates RainbowPlots, runs phylANOVA, outputs text, runs Brownielite, plots it
ANOVAData <- as.list(1:length(call))
for(i in rev(seq_along(call))){
  DataExtraction(OClink, birbs, ConsensusTree, call[i], mod=mod[i], fullrate=FullRate, RAIN=FALSE)
  ANOVAData[[i]] <- ANOVARun
  #plotting data
  sink(file = "Brownie.txt", append = TRUE, split = FALSE)
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], ".Brownie.pdf"))
  par(mfrow=c(2,2))
  BrowniePlotRates(dataset, paste0(call[i]))
  BrowniePlotRates(datasetFULL ,paste0(call[i],"FULL"))
  dev.off()
  sink(file=NULL)
}

ANOVAResults(ANOVAData)


setwd(file.path(dir, "DataWarehouse/"))
#Transition Plot Bones
pdf("TransitionBones.pdf", width=8.5, height=11)
plot(ConsensusTree, edge.width=2.5, cex=.7,
     label.offset = 1)
tiplabels(pch=ifelse(OpenClose == "closed", 19, 1),
          offset=.5)
dev.off()
