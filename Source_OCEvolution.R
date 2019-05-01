#Coded by Cristina Robinson
#Last Modified 5-1-2019

#Packages

#load packages
Packs <- c("maps", "ape", "phytools", "R.utils")
for(i in Packs){
  if(!eval(parse(text=paste0("require(",i,")")))){
    install.packages(i)
  }
  eval(parse(text=paste0("library(",i,")")))
}
rm(Packs)




#Functions that Load Data/Create the File Structure
LoadPrettyTree <- function(filename, Pretty=TRUE){
  #what it says: runs the fuctions necessary to create a nice consensus tree
  RawTree <- read.nexus(filename)
  ConsensusTree <- consensus.edges(RawTree)
  ConsensusTree <- PrepareConsensus(ConsensusTree)
  path <- PathFinder2(ConsensusTree)
  ConsensusTree <- PrettyConsensus(ConsensusTree, path)
  assign("Path", path, pos=1)#this will become a global variable
  return(ConsensusTree)
}
PrepareConsensus <-function(conse){
  #phytools becomes irate about edge lengths of 0
  conse$edge.length[conse$edge.length==0] = 0.0000001
  #the tree needs to be unrooted for phytools.
  #We used sayornis phoebe (suboscine) as the root and remove it
  #to create our unrooted tree
  conse=drop.tip(conse, 1, root.edge = 0)
  return(conse)
}
PrettyConsensus <- function(conse, path){
  if(missing(path)==TRUE){
    path <- PathFinder2(conse)
  }
#get the full path length based on indicies above, get tip index, load into database
  speciespath <- matrix(NA, ncol=3, nrow=length(path))
  for(i in 1:length(path)){
    endedge <- max(path[[i]])
    pathleng <- sum(conse$edge.length[path[[i]]])
    finaledge <- conse$edge.length[max(path[[i]])]
    speciespath[i,]<- cbind(endedge, pathleng, finaledge)
  }
  tipstop <- max(speciespath[,2])
  increments <- (speciespath[,2]-tipstop)*-1
  conse$edge.length[speciespath[,1]] <- conse$edge.length[speciespath[,1]] + increments
  return(conse)
}
PathFinder2 <- function(conse){
  #there is a minor issue with the consensus tree code in that it does not
  #ensure that the full length of each path is the same (so that all of tips end at the same place).
  
  #This code first recreates the node paths from root to tip for each species as a vector of nodes (the loops).
  #Once this exist, it calculates the length of each path and finds the max length.
  #It then subtracts each path length from the max length to get a value of how much
  #length needs to be added to have the tips line up.
  #That value is added to the last edge in the path (the one connecting to the tip).
  
  #Based on my testing, this does not affect future calculations and makes the trees look better.
  #If you are unsure, comment out the "PrettyConsensus" line of LoadPrettyTree;
  #the rest of the code will run just fine. :)
  
  #Pull the edges from the consensus and pull the first path
  edges <- conse$edge
  specindex <- which(edges[,2] < edges[1,1])
  path <- rep(list(NA),length(specindex))
  path[[1]] <- 1:specindex[1]
  names(path)[[1]] <- conse$tip.label[conse$edge[max(path[[1]]),2]]
  for(i in 2:length(specindex)){
#find new nodes to next tip and merge with last path
    new <- (specindex[i-1]+1):specindex[i]
    twopath <- c(edges[path[[i-1]]],edges[new])
#find the last common node and use it to make full node path
    common <- which(table(twopath) == 2)
#case 1 deals with the fork at the root as no path is shared
    if(names(common) == edges[1,1]){
      path[[i]] <- c(new)
      names(path)[[i]] <- conse$tip.label[conse$edge[max(path[[i]]),2]]
    }else{
      path[[i]] <- c(path[[i-1]][1:(common-1)],new)
      names(path)[[i]] <- conse$tip.label[conse$edge[max(path[[i]]),2]]
  }}
return(path)
}
Loadbirds <- function(filename, conse){
  #reading in the info  
  birdorder <- c(conse$tip.label)
  allinfo <- read.table(filename, sep=",", header=TRUE)
  #reorder to match the concensus tree
  allinfo[,1] <- gsub(" ", "_", allinfo[,"BirdtreeFormat"])
  reorder <- 1:length(birdorder)
  for(i in 1:length(birdorder)){
    reorder[i] <- which(allinfo[,"BirdtreeFormat"] == birdorder[i])
  }
  allinforeorder <- allinfo[reorder,]
  return(allinforeorder[,c(1:3,5:22,4)])
}
LoadOtherOC <- function(fileName, birbs){
  OCdata <- read.csv(fileName)
  reorder <- numeric(length=length(birbs$BirdtreeFormat))
  for(i in 1:length(birbs$BirdtreeFormat)){
    reorder[i] <- which( OCdata$BirdtreeFormat == birbs$BirdtreeFormat[i])
  }
  return(OCdata[reorder,1:4])
}
MakeFolderStructure <- function(dir){
  #create output directories and ANOVA printouts
  dir.create(file.path(dir, "DataWarehouse"))
  dir.create(file.path(dir, "DataWarehouse/Tri"))
  dir.create(file.path(dir, "DataWarehouse/Jackknife"))
  dir.create(file.path(dir, "DataWarehouse/MimidJackknife"))
  dir.create(file.path(dir, "DataWarehouse/SpecRates"))
  dir.create(file.path(dir, "DataWarehouse/ClosedLink"))
  for(i in seq_along(call2)){
    setwd(file.path(dir, "DataWarehouse/Jackknife"))
    dir.create(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
    setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
    if(file.exists("ANOVA.txt") == FALSE){
      file.create("ANOVA.txt")
    }
    if(file.exists("Brownie.txt") == FALSE){
      file.create("Brownie.txt")
    }
  }
  setwd(file.path(dir, "DataWarehouse"))
  if(file.exists("ANOVA.txt") == FALSE){
    file.create("ANOVA.txt")
  }
  if(file.exists("Brownie.txt") == FALSE){
    file.create("Brownie.txt")
  }
  setwd(file.path(dir, "DataWarehouse/Tri"))
  if(file.exists("ANOVA.txt") == FALSE){
    file.create("ANOVA.txt")
  }
  if(file.exists("TriBrownie.txt") == FALSE){
    file.create("TriBrownie.txt")
  }
  
    setwd(file.path(dir, "DataWarehouse/MimidJackknife"))
  if(file.exists("Brownie.txt") == FALSE){
    file.create("Brownie.txt")
  }
  setwd(file.path(dir, "DataWarehouse"))
}





#Functions that make the main Figures
QuickScatterBox <- function(vari, OC, title, labels=c("Song-Stable","Song-Plastic"), DOIT=TRUE){
  if(DOIT){
    pdf("QuickScatter.pdf")
    par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0) )
    
  }

  for(i in 1:ncol(vari)){
    current <- vari[,i]
    data <- sapply(levels(OC), function(type) current[which(OC == type)])
    for(j in seq_along(data)){
      remove <- which(is.na(data[[j]])==TRUE)
      if(length(remove)>0){
        data[[j]] <- data[[j]][-remove]
      }
    }
    ScatterBox(data,ylab=paste(title[i], "Repertoire Size"),
             xlab="",labels=labels,log='y',
             col1=rgb(.7,.1,.7,.5), col2 = rgb(.7,.1,.7,.5))
  }
  if(DOIT){
    dev.off()
    
  }
}
ScatterBox <- function(data, main="", sub="", xlab="", ylab="",
                       labels=FALSE,col1=NA,col2=rgb(0,0,0,.7),log=""){
  #make data a list, creat empty plot, get quantiles and get offset
  if(is.list(data)==FALSE){
    data <- list(data)
  }
  offset <- 1:length(data)-1
  plot(1,type='n',xlim=c(0,length(data)),ylim=c(min(unlist(data)),max(unlist(data))),xaxt='n',
       xlab=xlab, ylab=ylab, main=main, sub=sub,log=log, font.lab=2)
  axis(1,at=.5+offset,labels = labels, font=2)
  quan <- sapply(data, quantile, c(.25,.5,.75),type=8)
  #box
  rect(.25+offset,quan[1,],.5+offset,quan[3,],col = col1)
  #mean
  segments(.25+offset,quan[2,],.5+offset,quan[2,], lwd=2)
  
  
  #Scatter colors
  if(length(col2)!=length(data)){
    col2 <- rep(col2,length(data))
  }
  
  #get data for whiskers, scatter, and outliers
  up <- vector("numeric",length(data))
  down <- vector("numeric",length(data))
  index <- vector("list",length(data))
  jity <- vector("list",length(data))
  for(i in seq_along(data)){
    up[i] <- min(max(data[[i]]), quan[3,i] + 1.5*(quan[3,i]-quan[1,i]))
    down[i] <- max(min(data[[i]]), quan[1,i] - 1.5*(quan[3,i]-quan[1,i]))
    index <- which(data[[i]] > up[i] | data[[i]] < down[i])
    jity <- jitter(rep(.625,length(data[[i]])),8)+offset[i]
    type <- rep(16,length(data[[i]]))
    if(length(index) > 0){
      type[index] <- 1
    }
    #plot points
    points(jity, data[[i]], pch=type, col=col2[i])
  }
  
  #whisckerlength and block
  segments(.5+offset,c(quan[3,],quan[1,]),.5+offset,c(up,down),lty = 2)
  segments(.5+offset,c(up,down),.375+offset,c(up,down),lwd = 1)
}
DataExtraction <- function(OC, data, conse, vari, mod = "linear", cotitle="",
                           fullrate, RAIN=TRUE, ANOVA=TRUE, BROWNIE=TRUE,
                           DP=FALSE, Flip="rightwards"){

#cut NAs
  NAind <- which(is.na(data[,vari])==TRUE)
  if(length(NAind) != 0){
    data <- data[-NAind,]
    OC <- OC[-NAind]
    conse <- drop.tip(conse, NAind, root.edge = 0)
  }

#get Vari
  testvar <- data[,vari]
  names(testvar) <- data$BirdtreeFormat
  if(mod == "log"){
    testvar <- log(testvar)
    testvarpretty <- Squish(testvar)
  }else(testvarpretty <- testvar)
#plot it
  if(RAIN==TRUE){
    rates <- Rainbowplot(conse, testvarpretty, OC, title=paste0(vari,cotitle),DP,Flip)
    print("RainbowPlot Finished!")
  }else(rates <- RateMatx(ace(x=OC, phy=conse, type="discrete", model="ER")))
#test for diff
  if(ANOVA==TRUE){
    if(length(unique(testvar))>2){
      set.seed(49) #lucky number :)
      output<-phylANOVA(conse, OC, testvar, nsim = 1000)
      assign("ANOVARun", list(paste0(vari,cotitle),testvar,OC,output), pos=1)#this will become a global variable
      print("ANOVA Complete!")
    }else{
      warning(paste0(vari, " is binary; phylANOVA will not be generated"))
      assign("ANOVARun", NULL, pos=1)}
  }
#Brownie 
  if(BROWNIE == TRUE){
    BrownieDataGen(conse, OC, testvar, title=paste0(vari,cotitle), nsim=1300, rater=rates)
    if(missing(fullrate) == FALSE){
      print("Full Rate Sims started!")
      BrownieDataGen(conse, OC, testvar, title=paste0(vari,cotitle,"FULL"), nsim=1300, rater=fullrate)
    }
  }
}
Rainbowplot <- function(conse, testvar, OCdata, title="RainbowPlot",
                        doubleplot=FALSE, Flip="rightwards", ANC=FALSE,
                        CONT=FALSE){

  #cheating way to offset species labels
  if(Flip == "rightwards"){
    conse$tip.label <- paste("     ", conse$tip.label, sep = "")
    names(testvar) <- paste("     ", names(testvar), sep = "")
    Scootch <- 1.5
    Shove <- 0
    ShoveLess <- 0
  }else{
    conse$tip.label <- paste(conse$tip.label, "     ",sep = "")
    names(testvar) <- paste(names(testvar), "     ", sep = "")
    Scootch <- -.5
    Shove <- 23
    ShoveLess <- 14
  }
  names(OCdata) <- conse$tip.label
  #make pretty map
  if(length(levels(OCdata))==2){
    DI <- TRUE
    co <- c("black","white")
    levelLabel <- c("Song Stable","Song Plastic")
  }else{
    DI <- FALSE
    co <- c("black","grey50", "white")
    levelLabel <- c("Song Stable", "Longer Learning",  "Song Plastic")
  }
  if(length(unique(testvar))>2){
    Bin <- FALSE
    #continuous data
    AB<-contMap(conse,x=testvar,type="fan",fsize=0.05,lwd=2,plot=FALSE)
    names(AB$cols)=rev(names(AB$cols))
    #AB<-setMap(AB,colors=rev(c("#000000", "#6a3a2d", "#b84a29", "#e66731", "#ed9445", "#f6c15b", "#fdeb73", "#ffffd1", "#FFFFFF")))
    AB<-setMap(AB,colors=c("white","orchid1", "#530547"))
                 #c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))
  }else{
#Binary Data
    Bin <- TRUE
#test which model is better and pull the rates for the simmaps
    er <- ace(x=testvar, phy=conse, type="discrete", model="ER")
    ard <- ace(x=testvar, phy=conse, type="discrete", model="ARD")
    results <- anova(er,ard)
    ifelse(results$`Pr(>|Chi|)`[2] >= .05, rates <- RateMatx(er), rates <- RateMatx(ard))
    AB <- make.simmap(conse, testvar, Q=rates, nsim=1000, message=FALSE)
  }
#get the OC data
  testopenclosed = ace(x=OCdata, phy=conse, type="discrete", model="ER")
  returnrates <- RateMatx(testopenclosed)
  title <- gsub("[.]", " ", title)
  title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", title, perl=TRUE)
  if(!doubleplot && !ANC){
    pdf(paste(title, ".pdf",sep="")) 
  }else{
    if(doubleplot){
      if(title == "Syllable Rep Final"){
        title <- "Syllables"
      }else{
        title <- "Songs"}      
    }
  }
  if(Bin){
    AB <- densityMap(AB,fsize=0.6,legend=FALSE, direction=Flip)
    add.color.bar(leg=15,cols=rev(AB$cols),lims=c(0,1), font=2, fsize=.8,
                  digits=3,prompt=FALSE,x=3, title = "",
                  states = c("poly","mono"), subtitle=title,
                  y=(length(conse$tip.label)-.8*length(conse$tip.label)),
                  lwd=4,fsize=.7)
  }else{
    plot(AB,lwd=2,fsize=0.6, label.offset=.5, legend=FALSE, direction=Flip)
    add.color.bar(leg=15,rev(AB$cols), subtitle=title, title="",
                  lims=AB$lims,digits=3,prompt=FALSE,x=3.25+Shove, font=2, fsize=.8,
                  y=(length(conse$tip.label)-.8*length(conse$tip.label)),
                  lwd=4,fsize=.7,labels=exp(AB$lims))
  }
  if(!Bin){
    #place unlogged bar labels
    rect(-2.25+Shove,length(conse$tip.label)-.79*length(conse$tip.label),
         24+Shove, length(conse$tip.label)-.75*length(conse$tip.label),
         col="white", border=NA)
    text(c(2,20)+Shove,length(conse$tip.label)-.78*length(conse$tip.label),
         round(exp(AB$lims), digits=1), cex=.7) 
  }
  
  if(!ANC){
    tiplabels(pch = 21, bg = co[as.numeric(OCdata)], cex = 0.7, adj = Scootch)
    nodelabels(pie=testopenclosed$lik.anc, piecol=co, cex=.6)
    #legendthermo
    if(DI){
      rect(.4+ShoveLess,(length(conse$tip.label)-.935*length(conse$tip.label)),
           1.5+ShoveLess,(length(conse$tip.label)-.97*length(conse$tip.label)))
      rect(.4+ShoveLess,(length(conse$tip.label)-.90*length(conse$tip.label)),
           1.5+ShoveLess,(length(conse$tip.label)-.935*length(conse$tip.label)), col="black")
    }else{
      rect(.4+ShoveLess,(length(conse$tip.label)-.9625*length(conse$tip.label)),
           1.5+ShoveLess,(length(conse$tip.label)-.9925*length(conse$tip.label)))
      
      rect(.4+ShoveLess,(length(conse$tip.label)-.93*length(conse$tip.label)),
           1.5+ShoveLess,(length(conse$tip.label)-.9625*length(conse$tip.label)), col="grey50")
      
      rect(.4+ShoveLess,(length(conse$tip.label)-.90*length(conse$tip.label)),
           1.5+ShoveLess,(length(conse$tip.label)-.93*length(conse$tip.label)), col="black")
    }
    legend(2+ShoveLess,(length(conse$tip.label)-.9*length(conse$tip.label)), levelLabel,
           pch=21, cex=.75, col='black', pt.lwd = .5, pt.bg = co)
  }else{
    if(CONT){
      tiplabels(round(as.numeric(OCdata), digits=1), bg='white', frame='none', cex = 0.7, adj = .5, offset=2) 
    }else{
      tiplabels(pch = 21, bg = co[as.numeric(OCdata)], cex = 0.7, adj = Scootch)
      #legendthermo
      if(DI){
        rect(.4+ShoveLess,(length(conse$tip.label)-.935*length(conse$tip.label)),
             1.5+ShoveLess,(length(conse$tip.label)-.97*length(conse$tip.label)))
        rect(.4+ShoveLess,(length(conse$tip.label)-.90*length(conse$tip.label)),
             1.5+ShoveLess,(length(conse$tip.label)-.935*length(conse$tip.label)), col="black")
      }else{
        rect(.4+ShoveLess,(length(conse$tip.label)-.9625*length(conse$tip.label)),
             1.5+ShoveLess,(length(conse$tip.label)-.9925*length(conse$tip.label)))
        
        rect(.4+ShoveLess,(length(conse$tip.label)-.93*length(conse$tip.label)),
             1.5+ShoveLess,(length(conse$tip.label)-.9625*length(conse$tip.label)), col="grey50")

        rect(.4+ShoveLess,(length(conse$tip.label)-.90*length(conse$tip.label)),
             1.5+ShoveLess,(length(conse$tip.label)-.93*length(conse$tip.label)), col="black")
      }
      legend(2+ShoveLess,(length(conse$tip.label)-.9*length(conse$tip.label)), levelLabel,
             pch=21, cex=.75, col='black', pt.lwd = .5, pt.bg = co)
    }
    nodelabels(1:(length(testvar)-1), bg = "white", cex=.7)
  }
  if(!doubleplot  && !ANC){
    dev.off()
  }
  return(returnrates)
}
BrownieDataGen <- function(conse, OC, testvar, nsim=10, title, rater){
  print("Making simmaps for brownie.") 
  multisimmap <- make.simmap(tree=conse,x=OC, nsim=nsim, message=FALSE, Q=rater) 
  print(paste0("Simmaps generated", "  Starting a for loop with ", nsim, " loops."))
  BrownieWithWarnings(multisimmap, testvar, title, ncol(rater))
  print("End brownie loop")
}






BrownieWithWarnings <- function(multisimmap, testvar, title="title", levels=2){
  #generate empty dataframe for the data
  nsim <- length(multisimmap)
  if(nsim < 50){
    SaveIt <- nsim
  }else{
    SaveIt <- c(seq(50,nsim,by=50), nsim)
  }
  browniedata <- data.frame(Pval=numeric(nsim), ERRate=numeric(nsim), ERloglik=numeric(nsim),
                            ERace=numeric(nsim), ARDloglik=numeric(nsim),ARDace=numeric(nsim),
                            stringsAsFactors = FALSE)
  for(i in 1:levels){
    browniedata[paste0("ARDRate", i-1)] <- numeric(nsim)
  }
  browniedata['convergence'] <- character(nsim)
  #for each simmap, run the brownie code and extract results.  Save every 50 runs and at the end
  for (i in 1:nsim) {
    Results <- NULL
    withTimeout(expr={
      Results <- brownie.lite(multisimmap[[i]],testvar, maxit=750000)
      browniedata[i,1:(ncol(browniedata)-1)] <- c(Results[[11]], Results[[1]], Results[[4]],  Results[[2]],
                              Results[[9]], Results[[7]], Results[[6]])
      browniedata[i,ncol(browniedata)] <- Results[[12]]},
      timeout = 5, substitute = FALSE, onTimeout = "warning")
    #if a run takes too long, quit the brownie and throw NAs
    if(is.null(Results)==TRUE){print(paste("Timeout warning; data not generated for simmap", i, sep=" "));browniedata[i,1:8] <- NA}
    if(i %in% SaveIt){
      write.csv(file = paste(title,".csv",sep=""), browniedata);print("Saved data")
      print(paste("Brownie loop iteration",i,Sys.time()))
    }
  }
}#Coded with help from Kate Snyder

BrowniePlotRates <- function(dataset, title="BrowniePlot", col=c("blue","red"), Groups=c("Rate0", "Rate1"),
                             Xlim =c(min(XMins),
                                     max(XMaxes)) ){
  #data cleaning
  dataset<- dataset[!is.na(dataset$Pval),]
  dataset<- dataset[dataset$convergence == "Optimization has converged.",]
  
  #plots
  ARDs <- length(grep("ARDRate", colnames(dataset)))
  D <- list()
  YMaxes <- numeric(length=ARDs)
  YMins <- numeric(length=ARDs)
  XMaxes <- numeric(length=ARDs)
  XMins <- numeric(length=ARDs)
  for(i in 1:ARDs){
    D[[i]] <- density(dataset[,paste0('ARDRate', i-1)])
    YMaxes[i] <- max(D[[i]]$y)
    YMins[i] <- min(D[[i]]$y)
    XMaxes[i] <- max(D[[i]]$x)
    XMins[i] <- min(D[[i]]$x)
  }
  plot(D[[1]],col=col[1],
       xlim=Xlim,
       ylim=c(min(YMins),
              max(YMaxes)),
       main=title, xlab="Rates", font.lab=2,
       ylab="Number of Observations",
       cex.main=1, lwd=2)
  for(i in 2:ARDs){
    lines(D[[i]], col=col[i], lwd=2)
  }
  abline(v=dataset$ERRate[1])
  legend("topright", legend=Groups, col=col, lty=1, lwd=2)
  #stats
  MeanArd <- mean(dataset$ARDloglik)
  MeanER <- mean(dataset$ERloglik)
  pval <- round(pchisq(2*(MeanArd - MeanER),1,lower.tail=FALSE),digits=3)
  ifelse(pval == 0,pval <- "<0.001", pval <- paste0("=",pval))
  writeLines(title)
  writeLines(paste0("OneRate=", round(MeanER, digits = 4)))
  writeLines(paste0("TwoRates=", round(MeanArd, digits = 4)))
  writeLines(paste0("pVal", pval))
  writeLines(paste("",sep="\n\n"))
  writeLines(paste("",sep="\n\n"))
}

ANOVAResults <- function(ANOVAData){
#Prints out the ANOVA data with an alpha adjusted for
#multiple testing via the Holmes-Bonferroni Method
#for related independant variables
  AllInd <- length(ANOVAData)
  pval <- vector(mode="numeric", AllInd)
#Pull the pvals
  for(i in seq_along(ANOVAData)){
    pval[i] <- ANOVAData[[i]][[4]]$Pf
    names(pval)[i]<- ANOVAData[[i]][[1]]
  }
  Extra <- c(grep("[.]min", names(pval)), grep("[.]max", names(pval)))
#Get the unique tests (omit the max and min double checks)
  Unique <- c(1:AllInd)[-Extra]
#Calculate the corrected alphas
  newAlpha <- HolmesBonfer(Unique)
  PlaceMax <- cbind(order(pval[Unique]), 1:length(Unique))
  PlaceMax <- PlaceMax[order(PlaceMax[,1]),]
  newAlpha <- newAlpha[PlaceMax[,2]]
  AlphaCol <- rep(NA, length(pval))
  AlphaCol[Unique] <- newAlpha
  for(i in seq_along(AlphaCol)){
    if(is.na(AlphaCol[i])==TRUE){
      AlphaCol[i] <- AlphaCol[i-1]
    }
  }
  assign("CritAlpha", AlphaCol, pos=1)#this will become a global variable
  ANOVAPrinter(ANOVAData, AlphaCol)
  warning("This code assumes cols are in order of final, min, max as a unit for any vars.")
}





#Functions that are accessories for the above
meanminmaxNA <- function(variable, Birdtemp){
  NAind <- which(is.na(Birdtemp[,paste0(variable, ".min")])==TRUE)
  Birdtemp[NAind,paste0(variable, ".min")] <- Birdtemp[NAind,paste0(variable, ".final")]
  NAind <- which(is.na(Birdtemp[,paste0(variable, ".max")])==TRUE)
  Birdtemp[NAind,paste0(variable, ".max")] <- Birdtemp[NAind,paste0(variable, ".final")]
  return(Birdtemp)
}
Squish <- function(variable){
  quant <- quantile(variable)
  low<-quant[2] - (quant[3]-quant[2])
  high<-quant[4] + (quant[4]-quant[3])
  variable[which(variable > high)] <- high
  variable[which(variable < low)] <- low
  return(variable)
}
RateMatx <- function(rawrate){
  matsize <- sqrt(length(rawrate$index.matrix))
  
  #create matrix of correct size; make diag negative; multiple by rates
  qrates <- matrix(rep(1,matsize^2),matsize,matsize)
  diag(qrates) <- -1
  qrates <- qrates*rawrate$rates  
  #pulls the state names and sets as col and row names
  rownames(qrates) <- dimnames(rawrate$lik.anc)[[2]]
  colnames(qrates) <- dimnames(rawrate$lik.anc)[[2]]
  return(qrates)
}

SpeciesSimmap <- function(Conse, Fam, Path){
#make a simmap object
  SimPaths <- Conse
  numrow <- length(SimPaths$edge.length)
  SimPaths[[5]] <- as.list(rep(NA, numrow))
  SimPaths[6] <- list(cbind(Conse$edge.length, rep(0,numrow)))
  edgenames <-  paste(Conse$edge[,1],Conse$edge[,2], sep=",")
  rownames(SimPaths[[6]]) <- edgenames 
  colnames(SimPaths[[6]]) <- c(0,1)
  SimPaths[c(7:8)] <- NA
  names(SimPaths)[5:8] <- c("maps", "mapped.edge", "Q", "LogL")
  class(SimPaths) <- c("simmap", "phylo")
#make a holder for the map data
  edgelen <- SimPaths$edge.length
  names(edgelen) <- rep(0, numrow)

#turn all of nodes for Fam into the 1 state and load into simmap object
  OneIndex <- CladeHighlight(Conse, Fam, Path)
  names(edgelen)[OneIndex] <- "1"
  SimPaths[[6]][OneIndex,] <- cbind(rep(0,length(OneIndex)),SimPaths[[6]][OneIndex,1])
  for(i in 1:numrow){SimPaths[[5]][[i]] <- edgelen[i]}
return(SimPaths)
}
CladeHighlight <- function(Conse, Fam, Path){
#set the nodes assocaited with a family to 1
#---get the index for species of interest. pull the path values and count
#---the largest value present num of species times is the first common node
  Speci <- Conse$tip.label[Fam]
  numSpec <- length(Fam)
  Counts <- table(unlist(Path[Speci]))
  CommonPath <- which(Counts==numSpec)
  ComnNode <- max(as.numeric(names(CommonPath)))
#pull all nodes after the common node
  PostCNodes <- rep(list(NA),numSpec)
  for(i in 1:numSpec){
    WorkingData <- unlist(Path[Speci[i]])
    NewNodes <- WorkingData[which(WorkingData > ComnNode)]
    PostCNodes[[i]] <- NewNodes
  }
  AllOneNodes <- c(unique(unlist(PostCNodes)), ComnNode)
  return(AllOneNodes)
}

ANOVAPrinter <- function(ANOVAData, CritAlpha=.05){
  if(length(CritAlpha) != length(ANOVAData)){
    CritAlpha <- rep(CritAlpha, length(ANOVAData)) 
  }
  sink(file = "ANOVA.txt", append = TRUE, split = FALSE)
  
  if(length(levels(ANOVAData[[i]][[3]])) == 2){
    for(i in seq_along(ANOVAData)){
      writeLines(ANOVAData[[i]][[1]])
      writeLines(paste0("Closed=", round(mean(ANOVAData[[i]][[2]][which(ANOVAData[[i]][[3]]=="closed")]), digits = 4)))
      writeLines(paste0("Open=", round(mean(ANOVAData[[i]][[2]][which(ANOVAData[[i]][[3]]=="open")]), digits = 4)))
      writeLines(paste0("Fval=", round(ANOVAData[[i]][[4]]$F, digits = 4)))
      writeLines(paste0("pVal=", ANOVAData[[i]][[4]]$Pf))
      writeLines(paste0("Corrected Alpha=", round(CritAlpha[i], digits = 4)))
      writeLines(paste("",sep="\n\n"))
      writeLines(paste("",sep="\n\n"))
    }    
  }else{
    for(i in seq_along(ANOVAData)){
      writeLines(ANOVAData[[i]][[1]])
      writeLines(paste0("Closed=", round(mean(ANOVAData[[i]][[2]][which(ANOVAData[[i]][[3]]=="closed")]), digits = 4)))
      writeLines(paste0("Delayed=", round(mean(ANOVAData[[i]][[2]][which(ANOVAData[[i]][[3]]=="delayed-closed")]), digits = 4)))
      writeLines(paste0("Open=", round(mean(ANOVAData[[i]][[2]][which(ANOVAData[[i]][[3]]=="open")]), digits = 4)))
      writeLines(paste0("Fval=", round(ANOVAData[[i]][[4]]$F, digits = 4)))
      writeLines(paste0("pVal=", ANOVAData[[i]][[4]]$Pf))
      writeLines(paste0("Corrected Alpha=", round(CritAlpha[i], digits = 4)))
      
      writeLines("")
      writeLines("Closed")
      writeLines("Delayed")
      writeLines(paste0("Tval=", round(ANOVAData[[i]][[4]]$T[2], digits = 4)))
      writeLines(paste0("pVal=", ANOVAData[[i]][[4]]$Pt[2]))
      writeLines("")
      
      writeLines("")
      writeLines("Closed")
      writeLines("Open")
      writeLines(paste0("Tval=", round(ANOVAData[[i]][[4]]$T[3], digits = 4)))
      writeLines(paste0("pVal=", ANOVAData[[i]][[4]]$Pt[3]))
      writeLines("")
      
      writeLines("")
      writeLines("Open")
      writeLines("Delayed")
      writeLines(paste0("Tval=", round(ANOVAData[[i]][[4]]$T[6], digits = 4)))
      writeLines(paste0("pVal=", ANOVAData[[i]][[4]]$Pt[6]))
      ANOVAData[[i]][[4]]$T
      writeLines("")
      
      writeLines(paste("",sep="\n\n"))
      writeLines(paste("",sep="\n\n"))
    }
  }

  sink(file = NULL)
}
HolmesBonfer <- function(pVals){
  Alpha <- vector(mode="numeric", length(pVals))
  for(i in seq_along(pVals)){
    Alpha[i] <- (.05/(length(pVals) - i + 1))
  }
  return(Alpha)
}





####more stuff I had to add at the end
NodePlot <- function(tree, testVar, openClose, title,
                     OCAce, testVarAce, cont=FALSE, Log=TRUE){
  
  if(cont){
    Values <- round(OCAce$ace, digits=2)
    AceLength <- length(OCAce$ace)
    States <- "Length of Plasticity"
    OC <- Values
  }else{
    Values <- round(OCAce$lik.anc, digits=2)
    AceLength <- nrow(OCAce$lik.anc)
    States <- "State Likelihood"
    if(ncol(OCAce$lik.anc) == 2){
      OC <- paste0('S: ', Values[,1], ", P: ", Values[,2])  
    }else{
      OC <- paste0('E: ', Values[,1], ", D: ", Values[,2], ", P: ",Values[,3])
    }
  }
  Ylength <- seq(1,-.02,length=AceLength+1)[2:(AceLength+1)]
  
  
  Rainbowplot(tree, testvar=Squish(testVar), OCdata=openClose, ANC=TRUE, title=title,
              CONT=cont)
  #plot
  plot(0, type="n", ylim=c(0,1), xlim=c(0,1), xaxt='n', yaxt='n')
  #header
  text(x=c(0, .2, .8), 1.02,
       c("Node", States, "Trait"), adj=0, font=2, cex=.8)
  text(x=rep(c(.05, .25, .8), each=AceLength), 
       y=Ylength,
       c(paste0(1:AceLength, ': '),#node number
         OC, #OC
         round(exp(testVarAce$ace), digits=2)),#trait labels
       adj=0, cex=.7)  
}
MakeNodePlots <- function(testVar, names, trait, OCVariants,
                          conse, vConse, Log=TRUE){
  
  if(Log){
    testVar <- log(testVar)
  }
  names(testVar) <- names
  
  #make the OC vectors
  OpenClose <- factor(OCVariants$di)
  OpenCloseTri <- factor(OCVariants$tri)
  OpenCloseCont <- OCVariants$cont
  remove <- which(is.na(OCVariants$tri))
  OpenCloseTri <- OpenCloseTri[-remove]
  OpenCloseCont <- OpenCloseCont[-remove]
  
  #get the rep aces
  DITestVar <- ace(x=testVar, phy=conse, type="continuous", model="ER")
  CONTTRITestVar <- ace(x=testVar[-remove], phy=vConse, type="continuous", model="ER")
  
  #do the main aces
  DI <- ace(x=OpenClose, phy=conse, type="discrete", model="ER")
  TRI <- ace(x=OpenCloseTri, phy=vConse, type="discrete", model="ER")
  CONT <- ace(x=OpenCloseCont, phy=vConse, type="continuous", model="ER")
  
  
  
  
  
  NodePlot(conse, testVar, OpenClose, title=trait, OCAce=DI,
           testVarAce=DITestVar, cont=FALSE, Log=Log)
  NodePlot(vConse, testVar[-remove], OpenCloseTri, Log=Log,
           title=trait, OCAce=TRI, testVarAce=CONTTRITestVar, cont=FALSE)
  NodePlot(vConse, testVar[-remove], OpenCloseCont, title=trait,
           OCAce=CONT, testVarAce=CONTTRITestVar, cont=TRUE, Log=Log)

}
MakeAllNodePlots <- function(birbs, call, OCVariants, ConsensusTree){
  pdf('Nodes.pdf', width=7, height=7)
  layout(matrix(c(1,2), ncol=2, nrow=1), widths = c(.6,.45))
  for(i in c(2,5,8,11,14,17,18)){
    OCVariantTemp <- OCVariants
    CTree <- ConsensusTree
    TestVar <- birbs[,call[i]]
    names <- birbs$BirdtreeFormat
    
    remove <- which(is.na(TestVar))
    if(length(remove) > 0){
      TestVar <- TestVar[-remove]
      names <- names[-remove]
      OCVariantTemp <- OCVariantTemp[-remove,] 
      CTree <- drop.tip(CTree, remove)
    }
    remove2 <- which(is.na(OCVariantTemp$tri))
    if(length(remove2) > 0){
      VTree <- drop.tip(CTree, remove2)
    }else{
      VTree <- CTree
    }
    MakeNodePlots(testVar=TestVar,
                  names=names,
                  trait=call[i],
                  OCVariants=OCVariantTemp,
                  conse=CTree,
                  vConse=VTree)
  }
  dev.off()
  
}
