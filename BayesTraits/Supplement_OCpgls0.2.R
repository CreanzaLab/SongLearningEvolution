##Running BayesTraits Analysis##
########
#Coded by Kate T. Snyder
#Last Modified 5-2-2019
#Built using RStudio Version 1.1.456
#R Version 3.3.1
#
#ape_4.1  phytools_0.5-38   maps_3.1.0   nlme_3.1-137
#
########
########
########


pglsfunc <- function() { 
  require(ape)
  require(geiger)
  require(nlme)
  require(phytools)
df <- read.csv("OCPaperDataAllbtw.csv", stringsAsFactors = FALSE)
OCtree <- read.nexus("OCtreeHack1000.nex")

df$BirdtreeFormat <- as.character(df$BirdtreeFormat)
df$BirdtreeFormat[which(df$BirdtreeFormat == "Philesturnus_rufusater")] <- "Philesturnus_carunculatus"



songfeats <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Continuity", "Song.rate")
pglsresults <- set.seed(10)
for (i in c(1:7)) {
  trait1 <- "cont"
  trait2 <- songfeats[i]
  print(trait2)
subsetdf <- df[which(df[,trait1] != "-" & df[,trait2] != "-"),]
speciesToDrop <- df$BirdtreeFormat[which(df[,trait1] == "-" | df[,trait2] == "-")]
tipsToDrop <- which(OCtree$tip.label %in% speciesToDrop)
subsettree <- drop.tip(OCtree,tipsToDrop)

 trait1vec <- log(as.numeric(subsetdf[,trait1]))
# trait1vec <- as.numeric(subsetdf[,trait1])
 trait2vec <- log(as.numeric(subsetdf[,trait2]))
 subsetdf[,trait1] <- trait1vec
 subsetdf[,trait2] <- trait2vec
 row.names(subsetdf) <- subsetdf$BirdtreeFormat
if (i == 1) {
#pglsModel <- gls(cont ~ Syllable.rep.final, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
  pglsModel <- gls(Syllable.rep.final ~ cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
} else if (i == 2) {
#pglsModel <- gls(cont ~ Syll.song.final, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
  pglsModel <- gls(Syll.song.final ~ cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
} else if (i == 3) {
  cor <- corPagel(1,phy = subsettree)
#pglsModel <- gls(cont ~ Song.rep.final, data = subsetdf, correlation = cor, method = "ML")  #Setting lambda equal to 1 but not fixing lambda caused an error in the search for the most likely lambda for this song feature. We manually iterated through lambdas to generate this starting point and thus avoid the error.
pglsModel <- gls(Song.rep.final ~ cont, data = subsetdf, correlation = cor, method = "ML") 
} else if (i == 4) {
pglsModel <- gls(Duration.final ~cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
#pglsModel <- gls(cont~Duration.final, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
} else if (i == 5) {
pglsModel <- gls(Interval.final ~ cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
#  pglsModel <- gls(cont~Interval.final, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
} else if (i == 6) {
pglsModel <- gls(Continuity ~ cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
#  pglsModel <- gls(cont~Continuity, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
} else if (i == 7) {
pglsModel <- gls(Song.rate ~ cont, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
#pglsModel <- gls(cont~Song.rate, data = subsetdf, correlation = corPagel(1,phy = subsettree), method = "ML")
}
 

pglssumm <- summary(pglsModel)

lambda <- pglssumm$modelStruct$corStruct[1]
tableout <- pglssumm$tTable
temppglsresults <- cbind(rep(trait1,length(rownames(tableout))),rep(trait2,length(rownames(tableout))), rownames(tableout), tableout, lambda)
pglsresults <- rbind(pglsresults,temppglsresults)

}
write.csv(pglsresults, file = paste(Sys.Date(),"OCpglsLinearLearning.csv"))


}


testlambda <- function() {  ###Testing lambdas for song rep
lambda <- seq(0.01,1.00,0.01)
loglik <- c()
lambdavec <-c()
for (j in 1:length(lambda)) {
  cor <- corPagel(lambda[j], phy = subsettree, fixed = TRUE)
  fit <- gls(cont ~ Song.rep.final, data = subsetdf, correlation = cor, method = "ML")
  loglik <- c(loglik,fit$logLik)
  lambdavec <- c(lambdavec,fit$modelStruct$corStruct[1])
}

cor = corPagel(phy = subsettree)
pglsModel <- gls(cont ~ Song.rep.final, data = subsetdf, correlation = cor, method = "ML")

}