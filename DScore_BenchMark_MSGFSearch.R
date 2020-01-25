########NQR+0.984016 search; Main (First) search##########

###########general parameter
fas<-file.choose() #choose .fasta file to run against
basename(fas)

library("MSGFplus")
par <- msgfPar(
  database = fas,               #fasta file
  tolerance = '20 ppm',         ##Parent mass tolerance
  isotopeError = c(-1,2),       #isotope error range, e.g. "-t 20ppm -ti -1,2" tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2.
  tda = TRUE,                   #Target/Decoy search mode
  fragmentation = 3,            #1.CID,2.ETD,3.HCD
  enzyme = 1,                   #Trypsin
  instrument = 3,               #1. High-res LTQ, 2.TOF, 3.QExactive
  protocol = 0,                 #no protocol
  ntt = 1,                      ##Number of tryptic termini,The number of peptide termini that must have been cleaved by the enzyme (default 1)
  lengthRange = c(6,50),        #(Minimum peptide length,Maximum peptide length) to consider
  chargeRange = c(2,5),         #(Minimum precursor charge,Maximum precursor charge) to consider
  matches = 1                   ##Number of matches per spectrum to be reported
)                            #setup par
show(par)

#Amino Acid Modification Examples
# Specify static modifications using one or more StaticMod= entries
# Specify dynamic modifications using one or more DynamicMod= entries
# Modification format is:
# Mass or CompositionString, Residues, ModType, Position, Name (all the five fields are required).
# CompositionString can only contain a limited set of elements, primarily C H N O S or P
#
# Examples:
#   C2H3N1O1,  C,  fix, any,         Carbamidomethyl    # Fixed Carbamidomethyl C (alkylation)
#   O1,        M,  opt, any,         Oxidation          # Oxidation M
#   15.994915, M,  opt, any,         Oxidation          # Oxidation M (mass is used instead of CompositionStr)
#   H-1N-1O1,  NQ, opt, any,         Deamidated         # Negative numbers are allowed.
#   CH2,       K,  opt, any,         Methyl             # Methylation K
#   C2H2O1,    K,  opt, any,         Acetyl             # Acetylation K
#   HO3P,      STY,opt, any,         Phospho            # Phosphorylation STY
#   C2H3NO,    *,  opt, N-term,      Carbamidomethyl    # Variable Carbamidomethyl N-term
#   H-2O-1,    E,  opt, N-term,      Glu->pyro-Glu      # Pyro-glu from E
#   H-3N-1,    Q,  opt, N-term,      Gln->pyro-Glu      # Pyro-glu from Q
#   C2H2O,     *,  opt, Prot-N-term, Acetyl             # Acetylation Protein N-term

#######specify modifications
###Static
mods(par)[[1]]<-msgfParModification(name='Carbamidomethyl', #Carbamidomethyl
                                    composition='C2H3N1O1',
                                    residues='C',           #Cysteine C (alkylation)
                                    type='fix',             #fixed modification
                                    position='any')         #modification position


###Variable
#common
mods(par)[[2]] <- msgfParModification(name = 'Oxidation',   #Oxidation
                                      mass = 15.994915,
                                      residues = 'M',       #methionine
                                      type = 'opt',         #variable modification
                                      position = 'any')

mods(par)[[3]] <- msgfParModification(name='Acetyl',
                                      composition='C2H2O1',
                                      residues='*',
                                      type='opt',
                                      position='N-term')    # Acetylation Protein N-term

mods(par)[[4]] <- msgfParModification(name='Gln->pyro-Glu',
                                      composition='H-3N-1',
                                      residues='Q',
                                      type='opt',
                                      position='N-term')    # Pyro-glu from Q

#Deamidation & Citrullination
mods(par)[[5]] <- msgfParModification(name = 'Deamidated',
                                      mass = 0.984016,
                                      residues = 'NQR',
                                      type = 'opt',
                                      position = 'any')    # Deamidation & Citrullination

nMod(par) <- 3 #Set max number of modifications per peptide #nmod():Get and set the modifications in msgfPar objects
mods(par)      #Get and set the modifications in msgfPar objects
show(par)
###Modifications specification done>>>>

#specify peak-list file
mymzf<-file.choose()
mymzf

#MSGF+ database search
idres <- runMSGF(par, mymzf, memory=2000)  #run msgfplus


########Non-NQR+0.984016 search; Reference (Ref) search##########

###########general parameter
fas<-file.choose() #choose .fasta file to run against
basename(fas)

library("MSGFplus")
par <- msgfPar(
  database = fas,               #fasta file
  tolerance = '20 ppm',         ##Parent mass tolerance
  isotopeError = c(-1,2),       #isotope error range, e.g. "-t 20ppm -ti -1,2" tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2.
  tda = TRUE,                   #Target/Decoy search mode
  fragmentation = 3,            #1.CID,2.ETD,3.HCD
  enzyme = 1,                   #Trypsin
  instrument = 3,               #1. High-res LTQ, 2.TOF, 3.QExactive
  protocol = 0,                 #no protocol
  ntt = 1,                      ##Number of tryptic termini,The number of peptide termini that must have been cleaved by the enzyme (default 1)
  lengthRange = c(6,50),        #(Minimum peptide length,Maximum peptide length) to consider
  chargeRange = c(2,5),         #(Minimum precursor charge,Maximum precursor charge) to consider
  matches = 1                   ##Number of matches per spectrum to be reported
)                            #setup par
show(par)

#######specify modifications
###Static
mods(par)[[1]]<-msgfParModification(name='Carbamidomethyl', #Carbamidomethyl
                                    composition='C2H3N1O1',
                                    residues='C',           #Cysteine C (alkylation)
                                    type='fix',             #fixed modification
                                    position='any')         #modification position


###Variable
#common
mods(par)[[2]] <- msgfParModification(name = 'Oxidation',   #Oxidation
                                      mass = 15.994915,
                                      residues = 'M',       #methionine
                                      type = 'opt',         #variable modification
                                      position = 'any')

mods(par)[[3]] <- msgfParModification(name='Acetyl',
                                      composition='C2H2O1',
                                      residues='*',
                                      type='opt',
                                      position='N-term')    # Acetylation Protein N-term

mods(par)[[4]] <- msgfParModification(name='Gln->pyro-Glu',
                                      composition='H-3N-1',
                                      residues='Q',
                                      type='opt',
                                      position='N-term')    # Pyro-glu from Q

#Deamidation & Citrullination
# mods(par)[[5]] <- msgfParModification(name = 'Deamidated',
#                                       mass = 0.984016,
#                                       residues = 'NQR',
#                                       type = 'opt',
#                                       position = 'any')    # Deamidation & Citrullination

nMod(par) <- 3 #Set max number of modifications per peptide #nmod():Get and set the modifications in msgfPar objects
mods(par)      #Get and set the modifications in msgfPar objects
show(par)
###Modifications specification done>>>>

#specify peak-list file
mymzf<-file.choose()
mymzf

#MSGF+ database search
idres <- runMSGF(par, mymzf, memory=2000)  #run msgfplus


########NQR+0.984016&NQR+1.22694 search; deamidation/citrullination FDR analysis#####

###########general parameter
fas<-file.choose() #choose .fasta file to run against
basename(fas)

library("MSGFplus")
par <- msgfPar(
  database = fas,               #fasta file
  tolerance = '20 ppm',         ##Parent mass tolerance
  isotopeError = c(-1,2),       #isotope error range, e.g. "-t 20ppm -ti -1,2" tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2.
  tda = TRUE,                   #Target/Decoy search mode
  fragmentation = 3,            #1.CID,2.ETD,3.HCD
  enzyme = 1,                   #Trypsin
  instrument = 3,               #1. High-res LTQ, 2.TOF, 3.QExactive
  protocol = 0,                 #no protocol
  ntt = 1,                      ##Number of tryptic termini,The number of peptide termini that must have been cleaved by the enzyme (default 1)
  lengthRange = c(6,50),        #(Minimum peptide length,Maximum peptide length) to consider
  chargeRange = c(2,5),         #(Minimum precursor charge,Maximum precursor charge) to consider
  matches = 1                   ##Number of matches per spectrum to be reported
)                            #setup par
show(par)

#Amino Acid Modification Examples
# Specify static modifications using one or more StaticMod= entries
# Specify dynamic modifications using one or more DynamicMod= entries
# Modification format is:
# Mass or CompositionString, Residues, ModType, Position, Name (all the five fields are required).
# CompositionString can only contain a limited set of elements, primarily C H N O S or P
#
# Examples:
#   C2H3N1O1,  C,  fix, any,         Carbamidomethyl    # Fixed Carbamidomethyl C (alkylation)
#   O1,        M,  opt, any,         Oxidation          # Oxidation M
#   15.994915, M,  opt, any,         Oxidation          # Oxidation M (mass is used instead of CompositionStr)
#   H-1N-1O1,  NQ, opt, any,         Deamidated         # Negative numbers are allowed.
#   CH2,       K,  opt, any,         Methyl             # Methylation K
#   C2H2O1,    K,  opt, any,         Acetyl             # Acetylation K
#   HO3P,      STY,opt, any,         Phospho            # Phosphorylation STY
#   C2H3NO,    *,  opt, N-term,      Carbamidomethyl    # Variable Carbamidomethyl N-term
#   H-2O-1,    E,  opt, N-term,      Glu->pyro-Glu      # Pyro-glu from E
#   H-3N-1,    Q,  opt, N-term,      Gln->pyro-Glu      # Pyro-glu from Q
#   C2H2O,     *,  opt, Prot-N-term, Acetyl             # Acetylation Protein N-term

#######specify modifications
###Static
mods(par)[[1]]<-msgfParModification(name='Carbamidomethyl', #Carbamidomethyl
                                    composition='C2H3N1O1',
                                    residues='C',           #Cysteine C (alkylation)
                                    type='fix',             #fixed modification
                                    position='any')         #modification position


###Variable
#common
mods(par)[[2]] <- msgfParModification(name = 'Oxidation',   #Oxidation
                                      mass = 15.994915,
                                      residues = 'M',       #methionine
                                      type = 'opt',         #variable modification
                                      position = 'any')

mods(par)[[3]] <- msgfParModification(name='Acetyl',
                                      composition='C2H2O1',
                                      residues='*',
                                      type='opt',
                                      position='N-term')    # Acetylation Protein N-term

mods(par)[[4]] <- msgfParModification(name='Gln->pyro-Glu',
                                      composition='H-3N-1',
                                      residues='Q',
                                      type='opt',
                                      position='N-term')    # Pyro-glu from Q

#Deamidation & Citrullination
mods(par)[[5]] <- msgfParModification(name = 'Deamidated',
                                      mass = 0.984016,
                                      residues = 'NQR',
                                      type = 'opt',
                                      position = 'any')    # Deamidation & Citrullination
#Mock Deamidation & Citrullination
mods(par)[[2]] <- msgfParModification(name = 'Mock_Deamidated',
                                      mass = c(1.022694),
                                      residues = 'NQR',
                                      type = 'opt',
                                      position = 'any')    #Mock Deamidation & Citrullination

nMod(par) <- 3 #Set max number of modifications per peptide #nmod():Get and set the modifications in msgfPar objects
mods(par)      #Get and set the modifications in msgfPar objects
show(par)
###Modifications specification done>>>>

#specify peak-list file
mymzf<-file.choose()
mymzf

#MSGF+ database search
idres <- runMSGF(par, mymzf, memory=2000)  #run msgfplus








