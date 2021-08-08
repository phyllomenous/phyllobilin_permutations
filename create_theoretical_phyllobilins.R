### Creates list of theoretical phyllobilins by adding combinations of modifications to a pyro RCC


##### LOAD LIBRARIES AND SET WORKING DIRECTORY #####
library(openxlsx)
library(stringr)
library(tidyverse)
library(tibble)
library(dplyr)

## set your working directory

yourwd <- ""
setwd(yourwd)

##### DEFINE INPUT FILES AND IMPORT DATA #####

# List of known phyllobilin modifications
modifications_file <- "20200609_PhyllobilinModifications.xlsx"
modifications.df <- read.xlsx(modifications_file)



##### DEFINE FUNCTIONS #####

allcombinations <- function(dataframe, rowstocombine, variabletocombine) { # create combinations of modifications
  modgroups.list <- c()
  for(i in rowstocombine) { 
    tempsubset.df <- subset(dataframe, Position == i)
    formulalist <- list(tempsubset.df[,variabletocombine])
    modgroups.list <- c(modgroups.list, formulalist)
  }
  mydataframe.df <- expand.grid(modgroups.list) # create data frame with catabolites, all with different combinations of modifications
  names(mydataframe.df) <- rowstocombine
  return(mydataframe.df)
}


create.catabolites <- function(dataframe, combinationvariable) {
  rowstocombine <- unique(dataframe$Position)[-1] # define columns with elements to combine (w/o C3-2/C12-3)
  createdcatabolites.df <- allcombinations(dataframe, rowstocombine, combinationvariable) # create first data frame of combinations (w/o C3-2/C12-3)
  createdcatabolites.df <- add_column(createdcatabolites.df, "C3-2/C12-3" = NA, .after = 2) # add empty NA column for C3-2/C12-3
  
  rowstocombine <- unique(dataframe$Position)[c(-2,-3)] # define columns with elements to combine (with C3-2/C12-3, but w/o C3-2 & C12-3)
  secondcombinations <- allcombinations(dataframe, rowstocombine, combinationvariable) # create second data frame of combinations (with C3-2/C12-3, but w/o C3-2 & C12-3)
  secondcombinations <- add_column(secondcombinations, "C3-2" = NA, .after = 0) # add empty NA column for C3-2
  secondcombinations <- add_column(secondcombinations, "C12-3" = NA, .after = 0) # add empty NA column for C12-3
  
  createdcatabolites.df <- rbind(createdcatabolites.df, secondcombinations) # combine the two data frames
  return(createdcatabolites.df)
}

## add a function to create phyllobilins from different skeleton masses, RCC, pFCC,...
#add.skeletonmass <- function(skeletonmass)




###    #########

## count how many atoms of an element are in a chemical formula
## used in function read.chemformula()
count.atoms <- function(chemformula, i) { 
  chemformula <- unlist(strsplit(chemformula, split = ""))
  numberofatoms <- c()
  while(TRUE) {
    if(isTRUE(str_detect(chemformula[i+1], pattern = "[:upper:]") & # if atom is the only one
              str_detect(chemformula[i], pattern = "[:upper:]"))) {
      numberofatoms <- 1
      break
    } else if(isTRUE(str_detect(chemformula[i+1], pattern = "\\d"))) { # if character is followed by number, record number
      digit <- chemformula[i+1]
      numberofatoms <- append(numberofatoms, digit)
      i <- i+1
    }
    else {
      break
    }
  }
  return(numberofatoms)
}

## read a chemical formula and return a dataframe stating which elements and how many atoms of them are in the formula
## used in calculate.mass()
read.chemformula <- function(chemformula) {
  if(length(chemformula) == 0) { # if chemformula is NA, add Atoms = NA and Number = NA
    masscalculation.df <- data.frame(Atoms = NA,
                                     Number = NA)
    return(masscalculation.df)
  }
  if(is.na(chemformula)) { # if chemformula is NA, add Atoms = NA and Number = NA
    masscalculation.df <- data.frame(Atoms = NA,
                                     Number = NA)
    return(masscalculation.df)
  }
  else{
    if(str_detect(chemformula,"-")) {
      negative <- -1
      chemformula <- str_replace(chemformula, "-", "")
    } else {
      negative <- 1
    }
    chemformula <- unlist(strsplit(as.character(chemformula), split = ""))
    masscalculation.df <- data.frame(Atoms = as.character(),
                                     Number = as.numeric())
    nameofcolumns <- names(masscalculation.df)
    for(i in 1:length(chemformula)) {
      if(isTRUE(str_detect(chemformula[i], pattern = "\\d|[:lower:]"))) { # skip iteration if it's a number or a lower case letter.
        next
      }
      if(isTRUE(str_detect(chemformula[i+1], pattern = "\\d|[:upper:]"))) { # if letter is followed by a number, it's an atom.
        atom <- chemformula[i]
        numberofatoms <- count.atoms(chemformula, i)
      }
      if(isTRUE(str_detect(chemformula[i+1], pattern = "[:lower:]"))) { # if letter is followed by lower case letter, it's an atom with two letter abbreviation.
        firstletter <- chemformula[i]
        secondletter <- chemformula[i+1]
        atom <- paste(firstletter, secondletter, sep = "")
        i <- i+1
        numberofatoms <- count.atoms(chemformula, i)
      }
      if(isTRUE(str_detect(chemformula[i], pattern = "[:upper:]") & i == length(chemformula))) { # if letter is followed by an upper case letter, it's an atom (only one).
        atom <- chemformula[i]
        numberofatoms <- 1
      }
      numberofatoms <- paste(numberofatoms, collapse = "") # collapse vector to combine seperate strings
      numberofatoms <- as.numeric(numberofatoms) # turn string of digits into number
      numberofatoms <- numberofatoms * negative
      temp.df <- data.frame(atom,
                            numberofatoms)
      names(temp.df) <- names(masscalculation.df)
      masscalculation.df <- rbind(masscalculation.df, temp.df)
    }
    if(nrow(masscalculation.df) == 0) {
      masscalculation.df <- data.frame(Atoms = NA,
                                       Number = NA)
      return(masscalculation.df)
    }
    masscalculation.df <- aggregate(masscalculation.df$Number, by = list(masscalculation.df$Atoms), FUN = sum)
    names(masscalculation.df) <- nameofcolumns
    return(masscalculation.df)
  }
}



## write chemical formula as string from data frame
write.chemformula <- function(dataframe) {
  if(length(which(is.na(dataframe))) > 0) {
    print("No valid input")
    return(NA)
  }
  else {
    formula <- NULL
    for(i in 1:nrow(dataframe)) {
      element <- toString(dataframe$Atoms[i])
      number <- toString(dataframe$Number[i])
      if(number == 1) {
        atoms <- element
      } else if(number == 0) {
        atoms <- ""
      } else {
        atoms <- paste(element, number, sep = "")
      }
      formula <- paste(formula, atoms, sep = "")
    }
    return(formula)
  }
}



## THIS FUNCTION WORKS ALSO FOR NA AND INVALID INPUT
## summarize several chemical formulas 
add.chemformula <- function(formulavector) {
  if(is_empty(formulavector)) {
    myresult <- NA
  } else if(isTRUE(unique(is.na(formulavector)))) {
    myresult <- NA
  } else {
    myresult <- data.frame()
    mynames <- c("Atoms", "Number")
    for(i in formulavector) {
      tempformula <- read.chemformula(i)
      myresult <- rbind(myresult, tempformula)
    }
    names(myresult) <- mynames
    myresult <- aggregate(myresult$Number, by = list(myresult$Atoms), FUN = sum)
    names(myresult) <- mynames
    myresult <- write.chemformula(myresult)
  }
  return(myresult)
}


##             #####



##   ###













##### CALL FUNCTIONS #####

masscombo.df <- create.catabolites(modifications.df, "Mass.difference") # add up masses of all modifications within a combination of modifications 
modcombo.df <- create.catabolites(modifications.df, "Modification") # create a "name combination" of all names of modifications within a combination of modifications
formulacombo.df <- create.catabolites(modifications.df, "Formula.difference") # create a "formula combination" of all formulas of modifications within a combination of modifications

effectivemasses.df <- add_column(masscombo.df,"C32H36N4O5" = 556.268570, .after = 0) # add a column at the first position with the basic skeletion mass
effectivemasses.df$TotalMass <- rowSums(effectivemasses.df[,1:ncol(effectivemasses.df)], na.rm = TRUE)

effectiveformula.df <- add_column(formulacombo.df, "Skeleton" = "C32H36N4O5", .after = 0) # add a column at the first position with the basic skeleton formula
effectiveformula.df$TotalFormula <- NA
effectiveformula.df[] <- lapply(effectiveformula.df, as.character)
for(i in 1:nrow(effectiveformula.df)) {
  tempvec <- as.character(as.vector(effectiveformula.df[i,]))
  effectiveformula.df$TotalFormula[i] <- add.chemformula(tempvec)
}

HypotheticalPhyllobilins.df <- data.frame(effectivemasses.df$TotalMass) # make new dataframe with first column being the total mass of a hypothetical phyllobilin
HypotheticalPhyllobilins.df <- cbind(HypotheticalPhyllobilins.df, effectiveformula.df$TotalFormula)
names(HypotheticalPhyllobilins.df) <- c("CataboliteMass", "CataboliteFormula") # change name of column
HypotheticalPhyllobilins.df <- cbind(HypotheticalPhyllobilins.df, modcombo.df)
HypotheticalPhyllobilins.df <- distinct(HypotheticalPhyllobilins.df) # remove duplicates

# Rename columns
myorder <- c("CataboliteMass", "CataboliteFormula", "C1", "C2-1", "C2/C4", "C3-2", "C6/C7", "C8-2", "C12-3", "C15", "C18", "C3-2/C12-3") # define order of the columns
HypotheticalPhyllobilins.df <- HypotheticalPhyllobilins.df %>% # reorder columns, so it has the same order as "HypotheticalFragments" data frame
  select(myorder)

write.xlsx(HypotheticalPhyllobilins.df, "20200609_HypotheticalPhyllobilins.xlsx")

view(HypotheticalPhyllobilins.df)

##### END OF SCRIPT #####


