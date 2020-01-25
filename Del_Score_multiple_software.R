library(dplyr)
library(stringr)
####Xtandem##############
#####import results of Xtandem
library(rTANDEM)
MCP_mod_Xtandem<- GetResultsFromXML(file.choose()) #import mod_Xtandem
MCP_ref_Xtandem<- GetResultsFromXML(file.choose()) #import ref_Xtandem

##
#####Extract PSMs from Mod_search, mass error<=20ppm
##
Temp1<-MCP_mod_Xtandem@peptides
Temp2<-MCP_mod_Xtandem@ptm
colnames(Temp1)
colnames(Temp2)
Temp<-left_join(Temp1,Temp2,by="pep.id")
####abs(DelM_ppm)<20
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=(Temp$delta/Temp$mh)*1e6)%>%
  filter(abs(ParentMassErrorPPM)<=20)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_ppm_Xtandem<-Temp%>%
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
  filter(spectrum.id%in% UC_ThrMod_Fhit_ppm_Xtandem$spectrum.id)

Combo_AB<-left_join(UC_ThrMod_Fhit_ppm_Xtandem,In_A_B,by="spectrum.id")
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
Deamidated_ThrMod_NoMod_Xtandem<-Temp%>%
  filter(str_detect(modification_A, "0.9840"))
colnames(Deamidated_ThrMod_NoMod_Xtandem)

##
####Delta Score Filtration
##
##devide Deamidated_ThrMod_NoMod to same_backbone_AB and diff_backbone_AB
Same_Backbone_ThrMod_NoMod_Xtandem<-Deamidated_ThrMod_NoMod_Xtandem%>%
  filter(Deamidated_ThrMod_NoMod_Xtandem$pepSeq_A==Deamidated_ThrMod_NoMod_Xtandem$pepSeq_B)
##
Diff_Backbone_ThrMod_NoMod_Xtandem<-Deamidated_ThrMod_NoMod_Xtandem%>%
  filter(Deamidated_ThrMod_NoMod_Xtandem$pepSeq_A!=Deamidated_ThrMod_NoMod_Xtandem$pepSeq_B)

###MassError-Del_Score diagram
Temp<-Same_Backbone_ThrMod_NoMod_Xtandem%>%
  filter(XTandemEvalue_A<=10)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For same-seq spectra

#
Temp<-Diff_Backbone_ThrMod_NoMod_Xtandem%>%
  filter(XTandemEvalue_A<=10)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For diff-seq spectra

#pick out confident deamidated/citrullinated same_seq PSMs
Filter_Same_Backbone_ThrMod_NoMod_Xtandem<-Same_Backbone_ThrMod_NoMod_Xtandem%>%
  filter(Del_Score>0&
          XTandemEvalue_A<=0.5&
           abs(MassError_A)<=5
         |
           Del_Score>4&
           XTandemEvalue_A<=0.5&
           abs(MassError_A)<=10
  )

#pick out confident deamidated/citrullinated diff_seq PSMs
Filter_Diff_Backbone_ThrMod_NoMod_Xtandem<-Diff_Backbone_ThrMod_NoMod_Xtandem%>%
  filter(    Del_Score>2&
             XTandemEvalue_A<=0.5&
             abs(MassError_A)<=5)

##
#######For data with manual inspection annotation
##
####link dataset with Manual inspection annotation in MCP paper https://doi.org/10.1074/mcp.RA118.000696

#!Step: import Supplementary Table S1 by readr from paper: https://doi.org/10.1074/mcp.RA118.000696, title: Mining the Human Tissue Proteome for Protein Citrullination
colnames(Deamidated_ThrMod_NoMod_Xtandem)

#
Temp1<-Deamidated_ThrMod_NoMod_Xtandem%>%
  select(spectrum.id,
         spectrum.mh,
         pepSeq_A,
         type_A,
         modification_A,
         pepSeq_B,
         Del_Score,
         MassError_A,
         XTandemEvalue_A)
str(Temp1)
#!Step: import meta_data_ID by readr
meta_data_ID$mgf_spectrum_index<-as.numeric(meta_data_ID$mgf_spectrum_index)
str(meta_data_ID)
Temp1<-left_join(Temp1,meta_data_ID,by=c("spectrum.id"="mgf_spectrum_index"))

#
colnames(candidate_CitPep_MCP)
Temp2<-candidate_CitPep_MCP%>%
  select(`Scan number [raw file]`,MCP_PrecursorMZ=`Observed m/z`,
         `Manual inspection`,
         MCP_PepSeq=`Peptide sequence`,
         MCP_modifications=`Variable modifications identified by spectrum`)
#
colnames(Temp1)
colnames(Temp2)

#
Temp1$linking_id<-str_extract(Temp1$`spectrum title`,"[0-9]+")
Temp2$linking_id<-str_extract(Temp2$`Scan number [raw file]`,"[0-9]+")

Temp1$instrument_PrecursorMZ<-str_extract(Temp1$`spectrum title`,"[0-9]+[.][0-9]+@")
Temp1$instrument_PrecursorMZ<-str_replace(Temp1$instrument_PrecursorMZ,"@","")
Temp1$instrument_PrecursorMZ<-as.numeric(Temp1$instrument_PrecursorMZ)
class(Temp1$instrument_PrecursorMZ)

####
Temp<-left_join(Temp1,Temp2,by="linking_id")
Temp<-Temp%>%
  distinct(spectrum.id, .keep_all = T)
# Temp$PrecursorMZ<-round(Temp$PrecursorMZ,digits = 2)

Deamidated_ThrMod_NoMod_MI_Xtandem<-Temp
##


####comet##############
#####import results of comet
MCP_mod_comet #<-  import mod_comet by readr
MCP_ref_comet #<-  #import ref_comet by readr
colnames(MCP_mod_comet)
colnames(MCP_ref_comet)
##
#####Extract PSMs from Mod_search, mass error<=20ppm
##

####abs(DelM_ppm)<20
Temp<-MCP_mod_comet
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=((Temp$exp_neutral_mass-Temp$calc_neutral_mass)/Temp$calc_neutral_mass)*1e6)%>%
  filter(abs(ParentMassErrorPPM)<=20)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_ppm_comet<-Temp%>%
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
  filter(scan%in% UC_ThrMod_Fhit_ppm_comet$scan)

Combo_AB<-left_join(UC_ThrMod_Fhit_ppm_comet,In_A_B,by="scan")
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
Deamidated_ThrMod_NoMod_comet<-Temp%>%
  filter(str_detect(modification_A, "0.9840"))
colnames(Deamidated_ThrMod_NoMod_comet)

##
####Delta Score Filtration
##
##devide Deamidated_ThrMod_NoMod to same_backbone_AB and diff_backbone_AB
Same_Backbone_ThrMod_NoMod_comet<-Deamidated_ThrMod_NoMod_comet%>%
  filter(Deamidated_ThrMod_NoMod_comet$pepSeq_A==Deamidated_ThrMod_NoMod_comet$pepSeq_B)
##
Diff_Backbone_ThrMod_NoMod_comet<-Deamidated_ThrMod_NoMod_comet%>%
  filter(Deamidated_ThrMod_NoMod_comet$pepSeq_A!=Deamidated_ThrMod_NoMod_comet$pepSeq_B)

###MassError-Del_Score diagram
Temp<-Same_Backbone_ThrMod_NoMod_comet%>%
  filter(cometEvalue_A<=max(Same_Backbone_ThrMod_NoMod_comet$cometEvalue_A))
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For same-seq spectra

#
Temp<-Diff_Backbone_ThrMod_NoMod_comet%>%
  filter(cometEvalue_A<=max(Same_Backbone_ThrMod_NoMod_comet$cometEvalue_A))
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For diff-seq spectra

#pick out confident deamidated/citrullinated same_seq PSMs
Filter_Same_Backbone_ThrMod_NoMod_comet<-Same_Backbone_ThrMod_NoMod_comet%>%
  filter(Del_Score>0&
           cometEvalue_A<=0.5&
           abs(MassError_A)<=5
         |
           Del_Score>4&
           cometEvalue_A<=0.5&
           abs(MassError_A)<=10
  )

#pick out confident deamidated/citrullinated diff_seq PSMs
Filter_Diff_Backbone_ThrMod_NoMod_comet<-Diff_Backbone_ThrMod_NoMod_comet%>%
  filter(      str_count(modification_A,"0.984016")==1,
               Del_Score>2&
                 cometEvalue_A<=0.5&
                 abs(MassError_A)<=5)

##
#######For data with manual inspection annotation
##
####link dataset with Manual inspection annotation in MCP paper https://doi.org/10.1074/mcp.RA118.000696

#!Step: import Supplementary Table S1 by readr from paper: https://doi.org/10.1074/mcp.RA118.000696, title: Mining the Human Tissue Proteome for Protein Citrullination
colnames(Deamidated_ThrMod_NoMod_comet)

#
Temp1<-Deamidated_ThrMod_NoMod_comet%>%
  select(scan,
         pepSeq_A,
         peptide,
         modification_A,
         pepSeq_B,
         Del_Score,
         MassError_A,
         cometEvalue_A)
str(Temp1)
meta_data_ID$mgf_spectrum_index<-as.numeric(meta_data_ID$mgf_spectrum_index)
str(meta_data_ID)
Temp1<-left_join(Temp1,meta_data_ID,by=c("scan"="mgf_spectrum_index"))

#
colnames(candidate_CitPep_MCP)
Temp2<-candidate_CitPep_MCP%>%
  select(`Scan number [raw file]`,MCP_PrecursorMZ=`Observed m/z`,
         `Manual inspection`,
         MCP_PepSeq=`Peptide sequence`,
         MCP_modifications=`Variable modifications identified by spectrum`)
#
colnames(Temp1)
colnames(Temp2)

#
Temp1$linking_id<-str_extract(Temp1$`spectrum title`,"[0-9]+")
Temp2$linking_id<-str_extract(Temp2$`Scan number [raw file]`,"[0-9]+")

Temp1$instrument_PrecursorMZ<-str_extract(Temp1$`spectrum title`,"[0-9]+[.][0-9]+@")
Temp1$instrument_PrecursorMZ<-str_replace(Temp1$instrument_PrecursorMZ,"@","")
Temp1$instrument_PrecursorMZ<-as.numeric(Temp1$instrument_PrecursorMZ)
class(Temp1$instrument_PrecursorMZ)

####
Temp<-left_join(Temp1,Temp2,by="linking_id")
Temp<-Temp%>%
  distinct(scan, .keep_all = T)
# Temp$PrecursorMZ<-round(Temp$PrecursorMZ,digits = 2)

Deamidated_ThrMod_NoMod_MI_comet<-Temp
##


####maxquant##############
#####import results of comet
MCP_mod_maxquant #<-  import mod_maxquant
MCP_ref_maxquant #<-  #import ref_maxquant
colnames(MCP_mod_maxquant)
colnames(MCP_ref_maxquant)
##
#####Extract PSMs from Mod_search, mass error<=20ppm
##

####abs(DelM_ppm)<20
Temp<-MCP_mod_maxquant
Temp<-Temp%>%
  mutate(ParentMassErrorPPM=1e6*(Temp$`MS/MS m/z`-Temp$`m/z`)/Temp$`m/z`)%>%
  filter(abs(ParentMassErrorPPM)<=20)
colnames(Temp)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_ppm_maxquant<-Temp%>%
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
  filter(`MS/MS scan number`%in% UC_ThrMod_Fhit_ppm_maxquant$`MS/MS scan number`)

Combo_AB<-left_join(UC_ThrMod_Fhit_ppm_maxquant,In_A_B,by="MS/MS scan number")
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
    maxquentPEPValue_B=PEP.x,
    Del_Score,
    ParentMassErrorPPM_B=ParentMassErrorPPM.y,
    MassError_A
  )

#Extract candidate deamidated/citrullinated PSMs
Deamidated_ThrMod_NoMod_maxquant<-Temp%>%
  filter(str_detect(modification_A, "Deamidation"))
colnames(Deamidated_ThrMod_NoMod_maxquant)






