##Running BayesTraits Analysis##
########
#Coded by Kate T. Snyder
#Last Modified 5-2-2019
#Built using RStudio Version 1.0.136
#R Version 3.3.1
#
#ape_4.1  phytools_0.5-38   maps_3.1.0  btw_V1.0
#BayesTraitsV2 
########
########


#Must have BayesTraitsV2 installed. Must be performed on a Mac, as btw_V1.0 is only compatible with Mac.
#Must have .BayesTraitsPath set to the location of the program BayesTraitsV2 on computer, e.g.: .BayesTraitsPath <- "~/Documents/BayesTraits/BayesTraitsV2"
#Must have data files "OCPaperDataAll.csv" and "OCPaperDataAllbtw.csv" and tree file "OCtreeHack1000.nex"

source("SupplementOC_btwfunction0.1.R")
source("SupplementOC_subsetbirddata0.1.R")
source("SupplementOC_BayesPlots_choosebin0.1.R")
source("SupplementOC_BayesPlots0.2.R")
source("Supplement_OCpgls0.2.R")

#Run ALL standard BayesTraits analyses:
#Will return two PDFs (single 3-bin plot and combined file with plots for 2-5 bins) per test and one CSV per test
#7 tests total (OC + each of: Song Repertoire,Syllable repertoire,Syllables per song, Song duration, Intersong interval,Song rate,Song continuity)
cyclesongs(jackknife = FALSE,csvsout = TRUE,plotbt = TRUE,simnum = 50)

#If you only want to run a specific test:
#jackknife = TRUE to cycle through all families, removing each one in turn
#csvsout = FALSE to suppress .csv file outputs
#nsim = number of simulations. Warning that jackknife tests can take a long time - ~30min for 5 simulations. 
btwfunction("OC",SongParam = "Syllrep",plot = TRUE,jackknife = FALSE,csvsout = TRUE,nsim = 5)


#Perform PGLS test on all song features
#Will return one CSV file
pglsfunc()