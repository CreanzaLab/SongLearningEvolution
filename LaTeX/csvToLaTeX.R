#Setwds as necessary

if(!require("xtable")){
  install.packages("xtable")
}
library(xtable)
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
NoMate <- function(Tabler){
  remove <- which(Tabler[,1] %in% c("EPP", "Polygyny"))
  if(length(remove > 1)){
    Tabler <- Tabler[-remove,]  
  }
  return(Tabler)
}
SigStar <- function(Tabler, alpha=.05){
  starmie <- which(Tabler[,ncol(Tabler)]<alpha)
  if(length(starmie)>0){
    for(i in starmie){
      Tabler[i,ncol(Tabler)] <- paste0(Tabler[i,ncol(Tabler)],"*")
    }
  }
  return(Tabler)
}
SigStarCorrected <- function(Tabler){
  starmie <- numeric(0)
  for(i in 1:nrow(Tabler)){
    if(Tabler[i,ncol(Tabler)]<Tabler[i,ncol(Tabler)-1]){
      starmie <- c(starmie, i)
    }
  }
  if(length(starmie)>0){
    for(i in starmie){
      Tabler[i,ncol(Tabler)] <- paste0(Tabler[i,ncol(Tabler)],"*")
    }
  }
  return(Tabler)
}


setwd()#"DataWarehouse"
Data <- read.table("ANOVA.txt", stringsAsFactors = FALSE)
Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
Form <- Form[,-6]
Form <- gsub("Closed=", "", Form)
Form <- gsub("Open=", "", Form)
Form <- gsub("Fval=", "", Form)
Form <- gsub("pVal=", "", Form)
Form <- gsub("Alpha=", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
Form <- Form[,c(1:4,6,5)]
colnames(Form) <- c("Song Trait", "Song-Stable", "Song-Plastic", "F-Value", "Corrected alpha","p-Value")
Final <- c(1,4,7,10,13,16,17)
Finals <- Form[Final,]
Finals[,1] <- gsub(" Final", "", Finals[,1])
Finals <- Finals[order(as.numeric(Finals[,6])),]
Finals <- NoMate(Finals)
Finals <- SigStar(Finals)
print(xtable(Finals), include.rownames=FALSE)

MinMax <- Form[-Final,]
MinMax <- MinMax[order(as.numeric(MinMax[,6])),]
MinMax <- NoMate(MinMax)
MinMax <- SigStarCorrected(MinMax)
print(xtable(MinMax), include.rownames=FALSE)




Data <- read.table(file.path("Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- Form[seq(2,nrow(Form),by=2),]
Form <- gsub("FULL", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Two Rates","p-Value")
RowNames <- strsplit(Form[,1], " ")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
Final <- which(!(FinalNames %in% c("Min","Max")))
FormMain <- Form[Final,]
FormMain[,1] <- gsub("Final", " ", FormMain[,1])
FormMain <- FormMain[order(as.numeric(FormMain[,4])),]
FormMain[,1] <- gsub("  ", "", FormMain[,1])
FormMain <- NoMate(FormMain)
FormMain <- SigStar(FormMain)
print(xtable(FormMain), include.rownames=FALSE)

FormMinMax <- Form[-Final,] 
FormMinMax <- FormMinMax[order(as.numeric(FormMinMax[,4])),]
FormMinMax <- SigStar(FormMinMax)
print(xtable(FormMinMax), include.rownames=FALSE)







setwd()#"DataWarehouse/Tri"
Data <- read.table("ANOVA.txt", stringsAsFactors = FALSE)
Form <- matrix(Data[,1], ncol=20,byrow=TRUE)
Form <- Form[,-7]
Form <- gsub("Closed=", "", Form)
Form <- gsub("Open=", "", Form)
Form <- gsub("Fval=", "", Form)
Form <- gsub("pVal=", "", Form)
Form <- gsub("Alpha=", "", Form)
Form <- gsub("Delayed=", "", Form)
Form <- gsub("Tval=", "", Form)
Form <- gsub("Closed", "Stable", Form)
Form <- gsub("Delayed", "Longer", Form)
Form <- gsub("Open", "Plastic", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
post <- Form[,8:19]
Form <- Form[,c(1:5,7,6)]
colnames(Form) <- c("Song Trait", "Song-Stable", "Longer Learning", "Song-Plastic", "F-Value", "Corrected alpha","p-Value")
Final <- c(1,4,7,10,13,16,17)
Finals <- Form[Final,]
Finals[,1] <- gsub(" Final", "", Finals[,1])
reorder <- order(as.numeric(Finals[,7]))
Finals <- Finals[reorder,]
Finals <- NoMate(Finals)
Finals <- SigStar(Finals)
print(xtable(Finals), include.rownames=FALSE)

post <- matrix(t(post), ncol=4, byrow = TRUE)
post <- cbind.data.frame(rep(Form[,1], each=3),
                             post, stringsAsFactors=FALSE)
colnames(post) <- c("Song Trait", "State 1", "State 2",
                        "T-Value","p-Value")
post <- SigStar(post)
post[,1] <- gsub(" Final", "", post[,1])
sig <- which(grepl("[*]",Finals[,7]))
sigindex <- as.vector(sapply(Final[reorder][sig], function(x) x*3-(0:2)))
print(xtable(post[sigindex,]), include.rownames=FALSE)

MinMax <- Form[-Final,]
MMreorder <- order(as.numeric(MinMax[,6]))
MinMax <- MinMax[MMreorder,]
MinMax <- NoMate(MinMax)
MinMax <- SigStarCorrected(MinMax)
print(xtable(MinMax), include.rownames=FALSE)

sig <- which(grepl("[*]",MinMax[,7]))
sigindex <- as.vector(sapply((1:length(Form))[-Final][MMreorder][sig], function(x) x*3-(0:2)))
print(xtable(post[sigindex,]), include.rownames=FALSE)



Data <- read.table(file.path("TriBrownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- Form[seq(2,nrow(Form),by=2),]
Form <- gsub("FULL", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form <- gsub("-Tri", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Three Rates","p-Value")
RowNames <- strsplit(Form[,1], " ")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
Final <- which(!(FinalNames %in% c("Min","Max")))
FormMain <- Form[Final,]
FormMain[,1] <- gsub("Final", " ", FormMain[,1])
FormMain <- FormMain[order(as.numeric(FormMain[,4])),]
FormMain[,1] <- gsub("  ", "", FormMain[,1])
FormMain <- NoMate(FormMain)
FormMain <- SigStar(FormMain)
print(xtable(FormMain), include.rownames=FALSE)

FormMinMax <- Form[-Final,] 
FormMinMax <- FormMinMax[order(as.numeric(FormMinMax[,4])),]
FormMinMax <- SigStar(FormMinMax)
print(xtable(FormMinMax), include.rownames=FALSE)





setwd()#"DataWarehouse/"

Data <- read.table(file.path("BrownieTriDi.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- gsub("FULL", "", Form)
Form <- gsub("ThreeRates=", "", Form)
Form <- gsub("TwoRate=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- gsub(" final", "", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "Two Rates", "Three Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
Form <- SigStar(Form)
print(xtable(Form), include.rownames=FALSE)


Data <- read.table(file.path("BrownieTrinewDi.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- gsub("FULL", "", Form)
Form <- gsub("ThreeRates=", "", Form)
Form <- gsub("TwoRate=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- gsub(" final", "", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "Two Rates", "Three Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
Form <- SigStar(Form)
print(xtable(Form), include.rownames=FALSE)


Data <- read.table(file.path("BrownienewDi.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- gsub("FULL", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- gsub(" final", "", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Two Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
Form <- SigStar(Form)
print(xtable(Form), include.rownames=FALSE)


setwd()#"DataWarehouse/Jackknife"
for(i in list.files()){
  Data <- read.table(file.path(i,"ANOVA.txt"), stringsAsFactors = FALSE)
  Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
  Form <- Form[,-6]
  Form <- gsub("Closed=", "", Form)
  Form <- gsub("Open=", "", Form)
  Form <- gsub("Fval=", "", Form)
  Form <- gsub("pVal=", "", Form)
  Form <- gsub("Alpha=", "", Form)
  
  Form[,1] <- gsub("[.]", " ", Form[,1])
  RowNames <- strsplit(Form[,1], " ")
  FinalNames <- character(length(RowNames))
  for(j in seq_along(FinalNames)){
    FinalNames[j] <- tail(RowNames[[j]], n=1)
  }
  Form[,1] <- FinalNames
  Form <- gsub("finalNo", "", Form)
  Form <- gsub("rateNo", "", Form)
  Form <- Form[,c(1:4,6,5)]
  colnames(Form) <- c("Removed Family", "Song-Stable", "Song-Plastic", "F-Value", "Corrected alpha","p-Value")
  Form <- Form[order(as.numeric(Form[,6])),]
  Form <- NoMate(Form)
  Form <- SigStarCorrected(Form)
  print(i)
  print(xtable(Form), include.rownames=FALSE)
}

for(i in list.files()){
  Data <- read.table(file.path(i,"Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
  Form <- matrix(Data[(nrow(Data)/2+1):nrow(Data),1], ncol=4,byrow=TRUE)
  Form <- gsub("OneRate=", "", Form)
  Form <- gsub("TwoRates=", "", Form)
  Form <- gsub("pVal[=,<]?", "", Form)
  Form <- gsub("\\s*\\([^\\)]+\\)", "", Form)
  Form[,1] <- gsub("[.]final", " ", Form[,1])
  RowNames <- strsplit(Form[,1], " ")
  FinalNames <- character(length(RowNames))
  for(j in seq_along(FinalNames)){
    FinalNames[j] <- tail(RowNames[[j]], n=1)
  }
  Form[,1] <- FinalNames
  colnames(Form) <- c("Removed Family", "One Rate", "Two Rates","p-Value")
  Form <- Form[-1,]
  Form <- Form[order(as.numeric(Form[,4])),]
  Form <- SigStar(Form)
  print(xtable(Form, caption = i), include.rownames=FALSE)
}


ssetwd()#"DataWarehouse/MimidJackknife"
Data <- read.table(file.path("Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
RowNames <- strsplit(Form[,1], " ")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
Form[,1] <- FinalNames
Form <- gsub("_", " ",Form)
colnames(Form) <- c("Removed Mimid", "One Rate", "Two Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
Form <- SigStar(Form)
print(xtable(Form), include.rownames=FALSE)


setwd()#"DataWarehouse/ClosedLink"
Data <- read.table("ANOVA.txt", stringsAsFactors = FALSE)
Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
Form <- Form[,-6]
Form <- gsub("Closed=", "", Form)
Form <- gsub("Open=", "", Form)
Form <- gsub("Fval=", "", Form)
Form <- gsub("pVal=", "", Form)
Form <- gsub("Alpha=", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
Form <- Form[,c(1:4,6,5)]
colnames(Form) <- c("Song Trait", "Song-Stable", "Song-Plastic", "F-Value", "Corrected alpha","p-Value")
Final <- c(1,4,7,10,13,16,17)
Finals <- Form[Final,]
Finals[,1] <- gsub(" Final", "", Finals[,1])
Finals <- Finals[order(as.numeric(Finals[,6])),]
Finals <- NoMate(Finals)
Finals <- SigStarCorrected(Finals)

print(xtable(Finals), include.rownames=FALSE)



Data <- read.table(file.path("Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
RowNames <- strsplit(Form[,1], "[.]")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
dump <- which(FinalNames %in% c("final", "finalFULL"))
Form <- Form[dump,]
Form <- Form[seq(2,nrow(Form),by=2),]
Form <- gsub("FULL", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Two Rates","p-Value")
Form <- gsub(" Final", "", Form)
Form <- Form[order(as.numeric(Form[,4])),]
Form <- NoMate(Form)
Form <- SigStar(Form)
print(xtable(Form), include.rownames=FALSE)


setwd()#PGLS location
PGLS <- read.csv("OCpglsAll_noIntercept.csv", stringsAsFactors = FALSE)
PGLS[,1] <- gsub("[.]final", " ", PGLS[,1])
PGLS[,1] <- gsub("[.]", " ", PGLS[,1])
PGLS[,1] <- sapply(PGLS[,1], simpleCap)
PGLS[,7] <- 0
names(PGLS) <- c("Song Trait", "Slope", "Std Error", "T-Value", "p-Value", "lambda", "Corrected alpha")
PGLS <- PGLS[,c(1:3,6,4,7,5)]
PGLS <- PGLS[order(PGLS$`p-Value`),]
PGLS[,6] <- .05/(7:1)
for(i in 2:7){
  PGLS[,i] <- as.character(round(PGLS[,i], digits=4))
}
PGLS <- SigStar(PGLS)
print(xtable(PGLS), include.rownames=FALSE)