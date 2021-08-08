### Creates data frame of hypothetical MSMS fragment ions by combinations of
### modifications to fragments of an unmodifiec NCC


library(openxlsx)
library(stringr)
library(tidyverse)
library(plyr)


## set your working directory

yourwd <- ""
setwd(yourwd)


##### DEFINE INPUT FILES AND IMPORT DATA #####

# List of known phyllobilin modifications
modifications_file <- "20200609_PhyllobilinModifications.xlsx"
modifications.df <- read.xlsx(modifications_file)

reportedfragments_file <- "Chlorophyll catabolite MSMS fragments.xlsx"
reportedfragments.df <- read.xlsx((reportedfragments_file))
reportedfragments.df_columns <- c("Mass.loss", "Fragment", "Chemical.composition")

##### DEFINE PARAMETERS #####

rings <- c("A", 
           "A",
           "D",
           "A&B",
           "C&D")

ringformula <- c("C7H11NO", # non-hydroxylated DNCC, fragment d
                 "C8H11NO", # non-hydroxylated DNCC, fragment h
                 "C7H9NO",
                 "C16H20N2O2", # Rings A&B of non-hydroxylated pyro DNCC
                 "C16H18N2O3") # Rings C&D of RCC

ringmasses <- c(125.0840639804, 
                137.0840639804,
                123.0684139162,
                272.1524778966,
                286.1317424545)

fragmentationpositions <- c("d",
                            "h",
                            "c",
                            "g'",
                            "g")

## create dataframe with most basic losses / forms of catabolite fragment ions
fragmentskeletons.df <- data.frame(rings, fragmentationpositions, ringformula, ringmasses)



## create periodic table for mass calculation of chemical formulas
elementletters <- c("C",
                    "H",
                    "N",
                    "O")
elementmasses <- c(12.0000000000,
                   1.0078250321,
                   14.0030740052,
                   15.9949146221)

periodictable.df <- data.frame(elementletters, elementmasses) # combine the two vectors into data frame
names(periodictable.df) <- c("Element", "Mass")
str(periodictable.df)





##### DEFINE FUNCTIONS #####


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



## make all possible combinations between objects
makecombos <- function(dataframe, mygroups, combovariable) {
  modgroups.list <- c()
  for(j in mygroups) { # mygroups is the list of all positions at which modifications can occur
    positionsubset <- subset(dataframe, Position == j) # subset the data frame by modification positions
    templist <- list(positionsubset[,combovariable]) # make a list with all the remaining values within the chosen column (combovariable) 
    modgroups.list <- c(modgroups.list, templist) # append the main list which contains smaller lists (templists)
  }
  combo.df <- expand.grid(modgroups.list) # create data frame with fragments, all with different combinations of modifications
  colnames(combo.df) <- mygroups
  return(combo.df)
}









## create combinations of modifications for each ring
modificationcombos <- function(modifications.df, fragmentskeletons.df) {
  mydataframe.df <- data.frame(Ring = factor(),
                               Mass = numeric(),
                               Modifications = factor(),
                               # "C3-2/C12-3" = factor(),
                               # "C12-3" = factor(),
                               # "C3-2" = factor(),
                               # "C8-2" = factor(),
                               # "C6/C7" = factor(),
                               # "C2/C4" = factor(),
                               # "C1" = factor(),
                               # "C18" = factor(),
                               # "C15" = factor(),
                               Formula = character())
  fragments <- unique(fragmentskeletons.df$rings)
  for(i in fragments) {
    if(str_detect(i, pattern = "&")) {
      firstindex <- str_locate(i, pattern = "[:alpha:](?=&)")[1]
      secondindex <- str_locate(i, pattern = "(?<=&)[:alpha:]")[1]
      mystring <- unlist(strsplit(i, split = ""))
      firstring <- unlist(strsplit(mystring, split = ""))[firstindex]
      secondring <- unlist(strsplit(mystring, split = ""))[secondindex]
      ringsubset <- subset(modifications.df, modifications.df$Ring == firstring | 
                             modifications.df$Ring == secondring |
                             modifications.df$Ring == i)
    } else {
      ringsubset <- subset(modifications.df, modifications.df$Ring == i)
    }
    modpositions <- unique(ringsubset$Position)
    
    # make combinations of modification masses and sum them up
    combo.df <- makecombos(ringsubset, modpositions, "Mass.difference")
    if(ncol(combo.df) == 1) { # in case of fragment D (Ring D) there is only one position with two states = no combinations = dataframe has only one column
      combo.df$TotalMass <- combo.df[,1]
    } else {
      combo.df$TotalMass <- rowSums(combo.df[,1:ncol(combo.df)], na.rm =TRUE) # sum up all masses within one combination = Total mass
    }
    temp.df <- data.frame(i, # initiate a data frame and add fragment name in the first column, summed up mass of all mod combinations in the second column
                          combo.df$TotalMass)
    
    # make combinations of modification names
    moddataframe <- data.frame("C3-2/C12-3" = factor(),
                                "C12-3" = factor(),
                                "C3-2" = factor(),
                                "C8-2" = factor(),
                                "C6/C7" = factor(),
                                "C2/C4" = factor(),
                                "C1" = factor(),
                                "C18" = factor(),
                                "C15" = factor(),
                                "C2" = factor()  ## newly added 3.12.2020
                                )
    combo.df <- makecombos(ringsubset, modpositions, "Modification")
    modnameslist <- c()
    for(k in 1:nrow(combo.df)) {
      newname <- NULL
      for(l in 1:ncol(combo.df)) {
        modname <- combo.df[k,l]
        if(is.na(modname)) {
          modname <- ""
          newname <- paste(newname, modname, sep = "")
        } else {
          modname <- paste(modname, ", ", sep = "")
          newname <- paste(newname, modname, sep = " ")
        }
      }
      modnameslist <- append(modnameslist, newname)
    }
    combo.df$Modifications <- modnameslist
    temp.df <- cbind(temp.df, combo.df)
    

    # make combination of modification chem formulas
    combo.df <- makecombos(ringsubset, modpositions, "Formula.difference")
    modformulalist <- c()
    for(k in 1:nrow(combo.df)) {
      newformulavec <- c()
      for(l in 1:ncol(combo.df)) {
        modformula <- as.character(combo.df[k,l])
        newformulavec <- append(newformulavec, modformula)
      }
      newformula <- add.chemformula(newformulavec)
      modformulalist <- append(modformulalist, newformula)
    }
    combo.df$Formula <- as.character(modformulalist)
    temp.df <- cbind(temp.df, combo.df$Formula)
    mydataframe.df <- rbind.fill(mydataframe.df, temp.df)
  }
  # names(mydataframe.df) <- c("Ring", "MassOfModification", "NameOfModification", "Formula")
  return(mydataframe.df)
}






## add modification combinations to ring skeletons to create ring fragments
create.fragments <- function(fragmentinfo, modinfo) { 
  mydataframe.df <- data.frame(FragmentMass = numeric(),
                               Fragmentationposition = factor(),
                               Ring = factor(),
                               Modification = character(),
                               Formula = character())
  for(i in 1:nrow(fragmentinfo)) {
    skeletonmass <- fragmentinfo$ringmasses[i]
    skeletonformula <- as.character(fragmentinfo$ringformula[i])
    fragposition <- fragmentinfo$fragmentationpositions[i]
    ring <- fragmentinfo$rings[i]
    ringmodifications.df <- subset(modinfo, Ring == ring)
    for(j in 1:nrow(ringmodifications.df)) {
      modificationmass <- ringmodifications.df$MassOfModification[j]
      modificationformula <- as.character(ringmodifications.df$Formula[j])
      modificationposition.df <- ringmodifications.df[j,4:ncol(ringmodifications.df)]
      # fragmentmodification <-  ringmodifications.df$NameOfModification[j]
      fragmentmass <- skeletonmass + modificationmass
      fragmentformula <- add.chemformula(c(skeletonformula, modificationformula))
      temp.df <- data.frame(fragmentmass,
                            fragposition,
                            ring,
                            # fragmentmodification,
                            fragmentformula)
      temp.df <- cbind(temp.df, modificationposition.df)
      mydataframe.df <- rbind(mydataframe.df, temp.df)
    }
  }
  return(mydataframe.df)
}


##### CALL FUNCTIONS #####

# create data frame with different combinations of modifications
modinfo.df <- modificationcombos(modifications.df = modifications.df, 
                                 fragmentskeletons.df = fragmentskeletons.df)

modinfo.df <- modinfo.df %>% # reorder columns, so it has the same order as "HypotheticalPhyllobilins" data frame
  select(5,6,11,9,10,8,7,14,13,15,16,12)

colnames(modinfo.df)[1:3] <- c("Ring", "MassOfModification", "Formula")

# add modification combos to a skeleton catabolite
HypotheticalFragments <- create.fragments(fragmentskeletons.df, modinfo = modinfo.df)

# add fragments that have already been reported
reportedfragments.df <- reportedfragments.df[reportedfragments.df_columns]
names(reportedfragments.df) <- c("fragmentmass", "fragmentmodification", "fragmentformula")
HypotheticalFragments <- rbind.fill(HypotheticalFragments, reportedfragments.df)

write.xlsx(HypotheticalFragments, "20200615_NEW_HypotheticalFragments.xlsx")
view(HypotheticalFragments)


##### END OF SCRIPT #####





