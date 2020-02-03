library(dplyr)
library(stringr)
###########MSGF+##############
##
###Extract PSMs from mzid file of Mod_search, keep all
##
mzids<-file.choose() #Choose Results from Mod_MSGF+ search (.mzid file; MCP_cit_concatenated_mod.mzid)
#or
mzids<-"/Users/xavierwang/Desktop/Deamidation paper/revision/datasets from MCP citrullination paper/mzid_file/MCP_cit_concatenated_mod.mzid"
library("MSnID")
msnid <- MSnID(".")
msnid <- read_mzIDs(msnid, mzids) #tansfer .mzid file to S4 subject
#The following two function call create the new  numMisCleavages and numIrrCleabages columns in the MSnID object
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
#extact experimental mass error
msnid$ParentMassErrorPPM <- mass_measurement_error(msnid)
#-log10 transformed MS-GF+ Spectrum E-value,
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
PSM_Ori<-msnid@psms          #retrive psms which data type is dataframe from msnid file.
#######Deal with exported dataframe file
library(dplyr)
library(stringr)
PSM_Ori<-as_tibble(PSM_Ori)
PSM_Ori<-PSM_Ori%>%
  arrange(PSM_Ori$`spectrum title`) #sorting by scan number

#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_all_MSGF<-PSM_Ori %>%
  distinct(spectrumID, .keep_all = TRUE)   #!Output_A

##
###Extract PSMs from mzid file of ref_search, keep all PSMs
##
##
mzids<-file.choose() #Choose Results from Ref_MSGF+ search (.mzid file; MCP_cit_concatenated_ref.mzid)
#or
mzids<-"/Users/xavierwang/Desktop/Deamidation paper/revision/datasets from MCP citrullination paper/mzid_file/MCP_cit_concatenated_ref.mzid"
library("MSnID")
msnid <- MSnID(".")
msnid <- read_mzIDs(msnid, mzids) #tansfer .mzid file to S4 subject
#The following two function call create the new  numMisCleavages and numIrrCleabages columns in the MSnID object
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
msnid$ParentMassErrorPPM <- mass_measurement_error(msnid) #add DelM_ppm
#-log10 transformed MS-GF+ Spectrum E-value,
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
without_Ori<-msnid@psms          #retrive psms which data type is dataframe from msnid file.
#####tiding and sorting data
library(dplyr)
library(stringr)
without_Ori<-as_tibble(without_Ori)
without_Ori<-without_Ori%>%
  arrange(without_Ori$`spectrum title`) #sorting by scan number
##
UC_NoMod_Fhit<-without_Ori%>%
  distinct(spectrumID, .keep_all = TRUE)  #!Output_B

##
##############Linking PSMs from ref_search to mod_search by scan number&Delta Score Calculation
##
#calculation of systamic mass shift
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit%>%
  filter(abs(ParentMassErrorPPM)<=20&
           UC_NoMod_Fhit$`MS-GF:SpecEValue`<=1e-10)
SMS<-median(Temp$ParentMassErrorPPM)
hist(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit%>%
  filter(spectrumID%in% UC_ThrMod_Fhit_all_MSGF$spectrumID)

Combo_AB<-left_join(UC_ThrMod_Fhit_all_MSGF,In_A_B,by="spectrumID")

Temp<-Combo_AB %>%
  mutate(Del_Score=msmsScore.x-msmsScore.y,
         MassError_A=ParentMassErrorPPM.x-SMS
  )%>%
  select(spectrumID,
         `spectrum title`=`spectrum title.x`,
         Charge=chargeState.x,
         PrecursorMZ=experimentalMassToCharge.x,
         ParentMassErrorPPM_A=ParentMassErrorPPM.x,
         modification_A=modification.x,
         peptide_A=peptide.x,
         modification_B=modification.y,
         peptide_B=peptide.y,
         protein_A=accession.x,
         protein_B=accession.y,
         description_A=description.x,
         numIrregCleavages_A=numIrregCleavages.x,
         numMissCleavages_A=numMissCleavages.x,
         IsotopeError_A=IsotopeError.x,
         SpecEvalue_A=`MS-GF:SpecEValue.x`,
         SpecEvalue_B=`MS-GF:SpecEValue.y`,
         msmsScore_A=msmsScore.x,
         msmsScore_B=msmsScore.y,
         Del_Score,
         idFile=idFile.x,
         spectrumFile=spectrumFile.x,
         databaseFile=databaseFile.x,
         pepSeq_A=pepSeq.x,
         pepSeq_B=pepSeq.y,
         ParentMassErrorPPM_B=ParentMassErrorPPM.y,
         MassError_A,
         start_A=start.x,
         end_A=end.x
  )

#Extract candidate deamidated/citrullinated PSMs
Delta_Score_MSGF<-Temp
colnames(Delta_Score_MSGF)




####Xtandem##############
#####import results of Xtandem
library(rTANDEM)
MCP_mod_Xtandem<- GetResultsFromXML(mod_Xtandem<-file.choose()) #import mod_Xtandem
MCP_ref_Xtandem<- GetResultsFromXML(ref_Xtandem<-file.choose()) #import ref_Xtandem

##
#####Extract PSMs from Mod_search, Keep all
##
Temp1<-MCP_mod_Xtandem@peptides
Temp2<-MCP_mod_Xtandem@ptm
colnames(Temp1)
colnames(Temp2)
Temp<-left_join(Temp1,Temp2,by="pep.id")
####abs(DelM_ppm)<20
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=(Temp$delta/Temp$mh)*1e6)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_all_Xtandem<-Temp%>%
  distinct(spectrum.id, .keep_all = TRUE)   #!Output_A

##
###Extract PSMs from mzid file of ref_search, keep all PSMs#
##
Temp1<-MCP_ref_Xtandem@peptides
Temp2<-MCP_ref_Xtandem@ptm
colnames(Temp1)
colnames(Temp2)
Temp<-left_join(Temp1,Temp2,by="pep.id")
####keep all
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=(Temp$delta/Temp$mh)*1e6)%>%
  distinct(spectrum.id, .keep_all = TRUE)

#generate first-hit PSMs for without_Ori
UC_NoMod_Fhit_Xtandem<-Temp%>%
  distinct(spectrum.id, .keep_all = TRUE)   #!Output_B

##
####calculation of systamic mass shift
##
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit_Xtandem%>%
  filter(abs(ParentMassErrorPPM)<=20&
           expect.value<=0.05)
SMS<-median(Temp$ParentMassErrorPPM)
hist(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit_Xtandem%>%
  filter(spectrum.id%in% UC_ThrMod_Fhit_all_Xtandem $spectrum.id)

Combo_AB<-left_join(UC_ThrMod_Fhit_all_Xtandem,In_A_B,by="spectrum.id")
colnames(Combo_AB)

Temp<-Combo_AB %>%
  mutate(Del_Score=-log10(expect.value.x)-(-log10(expect.value.y)),
         MassError_A=ParentMassErrorPPM.x-SMS
  )%>%
  select(spectrum.id,
         Charge=spectrum.z.x,
         spectrum.mh=spectrum.mh.x,
         ParentMassErrorPPM_A=ParentMassErrorPPM.x,
         modification_A=modified.x,
         type_A=type.x,
         pepSeq_A=sequence.x,
         modification_B=modified.y,
         type_B=type.y,
         pepSeq_B=sequence.y,
         XTandemEvalue_A=expect.value.x,
         XTandemEvalue_B=expect.value.y,
         Del_Score,
         ParentMassErrorPPM_B=ParentMassErrorPPM.y,
         MassError_A
  )

#Extract candidate deamidated/citrullinated PSMs
Delta_Score_Xtandem<-Temp
colnames(Delta_Score_Xtandem)

##



####comet##############
#####import results of comet
MCP_mod_comet #<-  import mod_comet by readr, MCP_cit_concatenated.comet.txt
MCP_ref_comet #<-  #import ref_comet by readr
colnames(MCP_mod_comet)
colnames(MCP_ref_comet)
##
#####Extract PSMs from Mod_search, mass error<=20ppm
##

####abs(DelM_ppm)<20
Temp<-MCP_mod_comet
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=((Temp$exp_neutral_mass-Temp$calc_neutral_mass)/Temp$calc_neutral_mass)*1e6)
#generate first-hit all PSMs for PSM_Ori
UC_ThrMod_Fhit_all_comet<-Temp%>%
  distinct(scan, .keep_all = TRUE)   #!Output_A

##
###Extract PSMs from mzid file of ref_search, keep all PSMs#
##
####keep all
Temp<-MCP_ref_comet
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=((Temp$exp_neutral_mass-Temp$calc_neutral_mass)/Temp$calc_neutral_mass)*1e6)%>%
  distinct(scan, .keep_all = TRUE)

#generate first-hit PSMs for without_Ori
UC_NoMod_Fhit_comet<-Temp%>%
  distinct(scan, .keep_all = TRUE)   #!Output_B

##
####calculation of systamic mass shift
##
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit_comet%>%
  filter(abs(ParentMassErrorPPM)<=20&
           `e-value`<=0.05)
SMS<-median(Temp$ParentMassErrorPPM)
hist(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit_comet%>%
  filter(scan%in% UC_ThrMod_Fhit_all_comet$scan)

Combo_AB<-left_join(UC_ThrMod_Fhit_all_comet,In_A_B,by="scan")
colnames(Combo_AB)

Temp<-Combo_AB%>%
  mutate(msmsScore_A=-log10(Combo_AB$`e-value.x`),
         msmsScore_B=-log10(Combo_AB$`e-value.y`),
         Del_Score=msmsScore_A-msmsScore_B,
         MassError_A=ParentMassErrorPPM.x-SMS)
colnames(Temp)
Temp<-Temp%>%
  select(
    scan,
    Charge=charge.x,
    exp_neutral_mass_A=exp_neutral_mass.x,
    ParentMassErrorPPM_A=ParentMassErrorPPM.x,
    modification_A=modifications.x,
    peptide=modified_peptide.x,
    pepSeq_A=plain_peptide.x,
    pepSeq_B=plain_peptide.y,
    cometEvalue_A=`e-value.x`,
    cometEvalue_B=`e-value.y`,
    Del_Score,
    ParentMassErrorPPM_B=ParentMassErrorPPM.y,
    MassError_A
  )

#Extract candidate deamidated/citrullinated PSMs
Delta_Score_comet<-Temp
colnames(Delta_Score_comet)




####maxquant##############
#####import results of maxquant
MCP_mod_maxquant #<-  import mod_maxquant, evidence.txt
MCP_ref_maxquant #<-  #import ref_maxquant
colnames(MCP_mod_maxquant)
colnames(MCP_ref_maxquant)
##
#####Extract PSMs from Mod_search, mass error<=20ppm
##

####get all PSMs
Temp<-MCP_mod_maxquant
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=1e6*(Temp$`MS/MS m/z`-Temp$`m/z`)/Temp$`m/z`)
colnames(Temp)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_all_maxquant<-Temp%>%
  distinct(`MS/MS scan number`, .keep_all = TRUE)   #!Output_A

##
###Extract PSMs from mzid file of ref_search, keep all PSMs#
##
####keep all
Temp<-MCP_ref_maxquant
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=1e6*(Temp$`MS/MS m/z`-Temp$`m/z`)/Temp$`m/z`)

#generate first-hit PSMs for without_Ori
UC_NoMod_Fhit_maxquant<-Temp%>%
  distinct(`MS/MS scan number`, .keep_all = TRUE)  #!Output_B

##
####calculation of systamic mass shift
##
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit_maxquant%>%
  filter(abs(ParentMassErrorPPM)<=20&
           PEP<=0.05)
SMS<-median(Temp$ParentMassErrorPPM)
hist(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit_maxquant%>%
  filter(`MS/MS scan number`%in% UC_ThrMod_Fhit_all_maxquant $`MS/MS scan number`)

Combo_AB<-left_join(UC_ThrMod_Fhit_all_maxquant,In_A_B,by="MS/MS scan number")
colnames(Combo_AB)

Temp<-Combo_AB%>%
  mutate(msmsScore_A=-log10(Combo_AB$PEP.x),
         msmsScore_B=-log10(Combo_AB$PEP.y),
         Del_Score=msmsScore_A-msmsScore_B,
         MassError_A=ParentMassErrorPPM.x-SMS)
colnames(Temp)
Temp<-Temp%>%
  select(
    scan=`MS/MS scan number`,
    Charge=Charge.x,
    Mass_A=Mass.y,
    ParentMassErrorPPM_A=ParentMassErrorPPM.x,
    modification_A=Modifications.x,
    peptide=`Modified sequence.x`,
    pepSeq_A=Sequence.x,
    pepSeq_B=Sequence.y,
    maxquentPEPValue_A=PEP.x,
    maxquentPEPValue_B=PEP.y,
    Del_Score,
    ParentMassErrorPPM_B=ParentMassErrorPPM.y,
    MassError_A
  )

#Extract candidate deamidated/citrullinated PSMs
Delta_Score_maxquant<-Temp
colnames(Delta_Score_maxquant)

##


###!Step: import meta_data_ID>>>>>>>
#!Step: import Supplementary Table S1 (candidate_CitPep_MCP) by readr from paper: https://doi.org/10.1074/mcp.RA118.000696,>>
#>>title: Mining the Human Tissue Proteome for Protein Citrullination>>>
colnames(candidate_CitPep_MCP)

###append Delta Score from different probability score to each spectrum
###MSGF+ SpecEvalue
Temp1<-Delta_Score_MSGF%>%
  select(spectrumID,`spectrum title`,
         peptide_A,modification_A,
         pepSeq_A,pepSeq_B,
         SpecEvalue_A,
         MassError_A,
         Del_Score_From_MSGFSpecEvalue=Del_Score)
Temp2<-meta_data_ID%>%
  select(mgf_spectrum_index,spectrumID)
Temp3<-left_join(Temp1,Temp2,by="spectrumID")
colnames(Temp3)
###X!tandem Evalue
colnames(Delta_Score_Xtandem)
Temp4<-Delta_Score_Xtandem%>%
  select(spectrum.id,
         Del_Score_From_XTandemEvalue=Del_Score,
         Xtandem_pepSeq_A=pepSeq_A)
Temp5<-left_join(Temp3,Temp4,by=c("mgf_spectrum_index" ="spectrum.id"))
###comet Evalue
colnames(Delta_Score_comet)
Temp6<-Delta_Score_comet%>%
  select(scan,
         Del_Score_From_cometEvalue=Del_Score,
         comet_pepSeq_A=pepSeq_A)
Temp7<-left_join(Temp5,Temp6,by=c("mgf_spectrum_index" ="scan"))
###maxquant PEPvalue
colnames(Delta_Score_maxquant)
Temp8<-Delta_Score_maxquant%>%
  select(scan,
         Del_Score_From_maxquantPEPvalue=Del_Score,
         maxquant_pepSeq_A=pepSeq_A)
Temp9<-left_join(Temp7,Temp8,by=c("mgf_spectrum_index" ="scan"))
##
Temp<-Temp9
Delta_Score_From_DiffProbScore<-Temp
#output
write.csv(Temp,file.path(file.path("Output_DScore_Benchmark","Delta_Score_From_DiffProbScore.csv")))

##correlation plot
Temp<-Delta_Score_From_DiffProbScore%>%
  filter(str_detect(modification_A,"0.984016"),
         str_count(modification_A,"0.984016")==1,
         pepSeq_A==pepSeq_B,
         SpecEvalue_A<=1e-10,
         MassError_A<=5#,
         #pepSeq_A==Xtandem_pepSeq_A,
         #pepSeq_A==comet_pepSeq_A
  )
plot(Temp$Del_Score_From_MSGFSpecEvalue,Temp$Del_Score_From_XTandemEvalue)
plot(Temp$Del_Score_From_MSGFSpecEvalue,Temp$Del_Score_From_cometEvalue)
plot(Temp$Del_Score_From_XTandemEvalue,Temp$Del_Score_From_cometEvalue)
#output
write.csv(Temp,file.path(file.path("Output_DScore_Benchmark","Delta_Score_From_DiffProbScore_DeaSameSeqOnesiteMassError.csv")))
Delta_Score_From_DiffProbScore_filter<-Temp
########stat
###MSGF+ vs X!Tandem
colnames(Delta_Score_From_DiffProbScore_filter)
Temp<-Delta_Score_From_DiffProbScore_filter%>%
  filter(!is.na(Del_Score_From_MSGFSpecEvalue)&
           !is.na(Del_Score_From_XTandemEvalue))
Total<-length(Temp$spectrumID)
#
Temp1<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue>0&
           Del_Score_From_XTandemEvalue>0)
I<-length(Temp1$spectrumID)
#
Temp2<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue<0&
           Del_Score_From_XTandemEvalue>0)
II<-length(Temp2$spectrumID)
#
Temp3<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue<0&
           Del_Score_From_XTandemEvalue<0)
III<-length(Temp3$spectrumID)
#
Temp4<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue>0&
           Del_Score_From_XTandemEvalue<0)
IV<-length(Temp4$spectrumID)
##stat
stat<-cbind(Total,I,II,III,IV)


###MSGF+ vs Comet
colnames(Temp)
Temp<-Delta_Score_From_DiffProbScore_filter%>%
  filter(!is.na(Del_Score_From_MSGFSpecEvalue)&
           !is.na(Del_Score_From_cometEvalue))
Total<-length(Temp$spectrumID)
#
Temp1<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue>0&
           Del_Score_From_cometEvalue>0)
I<-length(Temp1$spectrumID)
#
Temp2<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue<0&
           Del_Score_From_cometEvalue>0)
II<-length(Temp2$spectrumID)
#
Temp3<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue<0&
           Del_Score_From_cometEvalue<0)
III<-length(Temp3$spectrumID)
#
Temp4<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue>0&
           Del_Score_From_cometEvalue<0)
IV<-length(Temp4$spectrumID)
##stat
stat<-cbind(Total,I,II,III,IV)

Temp2<-Temp%>%
  filter(Del_Score_From_MSGFSpecEvalue<0&
           Del_Score_From_cometEvalue>2)

Temp3<-Temp%>%
  filter(
    Del_Score_From_cometEvalue>2)


