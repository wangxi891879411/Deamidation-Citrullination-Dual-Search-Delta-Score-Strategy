###Extract PSMs from mzid file of Mock_Mod_search, mass error<=20ppm#######
##
mzids<-file.choose() #Choose Results of MSGF+ search (.mzid files)
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
  arrange(PSM_Ori$`scan number(s)`) #sorting by scan number
####abs(DelM_ppm)<20
PSM_Ori_ppm<-PSM_Ori%>%
  filter(abs(ParentMassErrorPPM)<=20)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_ppm<-PSM_Ori_ppm %>%
  distinct(spectrumID, .keep_all = TRUE)


###Extract PSMs from mzid file of ref_search, keep all PSMs#######
##
mzids<-file.choose() #Choose Results of MSGF+ search (.mzid files)
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
  arrange(without_Ori$`scan number(s)`) #sorting by scan number
without_Ori

UC_NoMod_Fhit<-without_Ori%>%
  distinct(spectrumID, .keep_all = TRUE)


##############Linking PSMs from ref_search to mock_mod_search by scan number&Delta Score Calculation####

#calculation of systamic mass shift
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit%>%
  filter(abs(ParentMassErrorPPM)<=20&
           Temp$`MS-GF:SpecEValue`<=1e-10)
SMS<-median(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit%>%
  filter(`scan number(s)`%in% UC_ThrMod_Fhit_ppm$`scan number(s)`)

Combo_AB<-left_join(UC_ThrMod_Fhit_ppm,In_A_B,by="scan number(s)")

Temp<-Combo_AB %>%
  mutate(Del_Score=msmsScore.x-msmsScore.y,
         Del_Raw_Score=`MS-GF:RawScore.x`-`MS-GF:RawScore.y`,
         MassError_A=ParentMassErrorPPM.x-SMS,
         #,Del_Del_ppm=ParentMassErrorPPM.x-ParentMassErrorPPM.y
         #Experimental_Mass_A=experimentalMassToCharge.x*chargeState.x,
         #dataset="Nanowell_2DLC_islet_1_RZCol_270bar_55cm_022318"
         dataset_index=i
  )%>%
  select(`scan number(s)`,
         FragMethod=AssumedDissociationMethod.x,
         #Experimental_Mass_A,
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
         RawScore_A=`MS-GF:RawScore.x`,
         RawScore_B=`MS-GF:RawScore.y`,
         Del_Raw_Score,
         idFile=idFile.x,
         spectrumFile=spectrumFile.x,
         databaseFile=databaseFile.x,
         dataset_index,
         pepSeq_A=pepSeq.x,
         pepSeq_B=pepSeq.y,
         ParentMassErrorPPM_B=ParentMassErrorPPM.y,
         MassError_A,
         Del_Raw_Score,
         start_A=start.x,
         end_A=end.x
  )

#Extract candidate deamidated/citrullinated PSMs
Deamidated_ThrMod_NoMod_FDR<-Temp%>%
  filter(str_detect(modification_A, "0.984016")|str_detect(modification_A, "1.022694"))
colnames(Deamidated_ThrMod_NoMod_FDR)


####devide Deamidated_ThrMod_NoMod to same_backbone_AB and diff_backbone_AB#####

Temp<-Deamidated_ThrMod_NoMod_FDR
##
Same_Backbone_ThrMod_NoMod_FDR<-Temp%>%
  filter(Temp$pepSeq_A==Temp$pepSeq_B)
##
Diff_Backbone_ThrMod_NoMod_FDR<-Temp%>%
  filter(Temp$pepSeq_A!=Temp$pepSeq_B)
###

#########FDR_SameSeq###############################
library(dplyr)
library(stringr)
range(Same_Backbone_ThrMod_NoMod_FDR$Del_Score)
FDR_SameSeq<- data.frame(Del_score=double(),Target=integer(),Decoy=integer(),FDR=double(),PSMs=integer()) #creat an empty dataframe to store results
j<-c(-6)                                                      #!!!start point of Del_score
for (i in 1:(round(21/0.1)+1)){                                #!!!21 is the span of Del_score, 0.1 is the increasement
  Filter_Same_Backbone_FDR<-Same_Backbone_ThrMod_NoMod_FDR%>%
    filter(Del_Score>j&
             ParentMassErrorPPM_A>=c(-5+SMS)&
             ParentMassErrorPPM_A<=c(5+SMS)&
             msmsScore_A>=10
    )
  Target<-Filter_Same_Backbone_FDR%>%
    filter(str_detect(modification_A,"0.984016")
           &
             !str_detect(modification_A,"1.022694")
    )
  Decoy<-Filter_Same_Backbone_FDR%>%
    filter(str_detect(modification_A,"1.022694")
           #&
           #!str_detect(modification_A,"0.984016")
    )                  #!!!more strict for decoy
  Target<-length(Target$`scan number(s)`)
  Decoy<-length(Decoy$`scan number(s)`)
  FDR<-100*(2*Decoy/(Target+Decoy))
  PSMs<-length(Filter_Same_Backbone_FDR$`scan number(s)`)
  FDR_SameSeq[i,]<-c(round(j,digits = 2),Target,Decoy,round(FDR,digits = 2),PSMs)
  j=j+0.1                                                       #!!!increasement of Del_score
}
dir.create("Output_DScore_Benchmark")
write.csv(FDR_SameSeq,file.path("Output_DScore_Benchmark","FDR_SameSeq.csv"))

########FDR_DiffSeq##############
Temp<-Diff_Backbone_ThrMod_NoMod_FDR%>%       ##!!!update!!!
  filter(str_count(modification_A,"0.984016")+str_count(modification_A,"1.022694")==1)   # !NEW to estimate 1-site diff_backbone
range(Temp$Del_Score)                       ##!!!update!!!
FDR_DiffSeq<- data.frame(Del_score=double(),Target=integer(),Decoy=integer(),FDR=double(),PSMs=integer()) #creat an empty dataframe to store results
j<-c(-13)                                                      #!!!start point of Del_score
for (i in 1:(round(25/0.1)+1)){                                #!!!25 is the span of Del_score, 0.1 is the increasement
  Filter_Diff_Backbone_FDR<-Temp%>%
    filter(Del_Score>j&
             ParentMassErrorPPM_A>=c(-10+SMS)&
             ParentMassErrorPPM_A<=c(10+SMS)&
             msmsScore_A>=8
    )
  Target<-Filter_Diff_Backbone_FDR%>%
    filter(str_detect(modification_A,"0.984016")
           &!str_detect(modification_A,"1.022694"))
  Decoy<-Filter_Diff_Backbone_FDR%>%
    filter(str_detect(modification_A,"1.022694")
    )
  Target<-length(Target$`scan number(s)`)
  Decoy<-length(Decoy$`scan number(s)`)
  FDR<-100*(2*Decoy/(Target+Decoy))
  PSMs<-length(Filter_Diff_Backbone_FDR$`scan number(s)`)
  FDR_DiffSeq[i,]<-c(round(j,digits = 2),Target,Decoy,round(FDR,digits = 2),PSMs)
  j=j+0.1                                                       #!!!increasement of Del_score
}

range(FDR_DiffSeq$FDR)
##
dir.create("Output_DScore_Benchmark")
write.csv(FDR_DiffSeq,file.path("Output_DScore_Benchmark","FDR_DiffSeq_1site.csv")) #update name


