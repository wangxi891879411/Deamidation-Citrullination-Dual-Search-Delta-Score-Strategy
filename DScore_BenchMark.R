###Extract PSMs from mzid file of Mod_search, mass error<20ppm#######
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
####abs(DelM_ppm)<20
PSM_Ori_ppm<-PSM_Ori%>%
  filter(abs(ParentMassErrorPPM)<20)
#generate first-hit ppm-filtered PSMs for PSM_Ori
UC_ThrMod_Fhit_ppm<-PSM_Ori_ppm %>%
  distinct(spectrumID, .keep_all = TRUE)   #!Output_A



###Extract PSMs from mzid file of ref_search, keep all PSMs#######
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


##############Linking PSMs from ref_search to mod_search by scan number&Delta Score Calculation####

#calculation of systamic mass shift
library(stringr)
library(dplyr)
Temp<-UC_NoMod_Fhit%>%
  filter(abs(ParentMassErrorPPM)<20&
           UC_NoMod_Fhit$`MS-GF:SpecEValue`<1e-10)
SMS<-median(Temp$ParentMassErrorPPM)
hist(Temp$ParentMassErrorPPM)

#link datasets
In_A_B<-UC_NoMod_Fhit%>%
  filter(spectrumID%in% UC_ThrMod_Fhit_ppm$spectrumID)

Combo_AB<-left_join(UC_ThrMod_Fhit_ppm,In_A_B,by="spectrumID")

Temp<-Combo_AB %>%
  mutate(Del_Score=msmsScore.x-msmsScore.y,
         Del_Raw_Score=`MS-GF:RawScore.x`-`MS-GF:RawScore.y`,
         MassError_A=ParentMassErrorPPM.x-SMS#,
         #,Del_Del_ppm=ParentMassErrorPPM.x-ParentMassErrorPPM.y
         #Experimental_Mass_A=experimentalMassToCharge.x*chargeState.x,
         #dataset="Nanowell_2DLC_islet_1_RZCol_270bar_55cm_022318"
         #dataset_index=i
  )%>%
  select(spectrumID,
         `spectrum title`=`spectrum title.x`,
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
         #dataset_index,
         pepSeq_A=pepSeq.x,
         pepSeq_B=pepSeq.y,
         ParentMassErrorPPM_B=ParentMassErrorPPM.y,
         MassError_A,
         Del_Raw_Score,
         start_A=start.x,
         end_A=end.x
  )

#Extract candidate deamidated/citrullinated PSMs
Deamidated_ThrMod_NoMod<-Temp%>%
  filter(str_detect(modification_A, "0.984016"))
colnames(Deamidated_ThrMod_NoMod)


#####add "#" after NQR to indicate deamidation/citrullination#############

Temp<-Deamidated_ThrMod_NoMod%>%
  select(spectrumID,pepSeq_A,peptide_A,modification_A)
str(Temp)
Temp<-Temp%>%
  dplyr::mutate(Tidy_Modification=gsub(pattern = "\\)",replacement = "   )",x=Temp$modification_A)) #tidy strings in Modification.x

head(Temp$Tidy_Modification)

###extract modified sites
library(dplyr)
location_all<-str_locate_all(Temp$Tidy_Modification, "0.984016 \\(")
location_all

location_all<-lapply(location_all, function(x)
  as.data.frame(x))
location_all

location_all<-lapply(location_all, function(x)
  x[["start"]])
location_all

#count the number of modification sites
Temp<-Temp%>%
  dplyr::mutate(num_modification = sapply(location_all, function(x) ifelse(x[1] > 0, length(x), 0)))
#OR#Temp<-Temp%>%
#dplyr::mutate(num_modification=apply(Temp[,c("Modification_site_1","Modification_site_2","Modification_site_3")],1,function(x)
#length(x)-sum(is.na(x))))

P_1<-sapply(location_all, function(x)
  x[1])
#P_1[is.na(P_1)]<-0
P_1

P_2<-sapply(location_all, function(x)
  x[2])
#P_2[is.na(P_2)]<-0
#P_2
#or
#P_2<-replace(P_2,is.na(P_2),0)
#P_2
#or
#library(tidyr)
#P_2<-replace_na(P_2,0)
P_2


P_3<-sapply(location_all, function(x)
  x[3])
#P_3[is.na(P_3)]<-0
P_3

Modification_site_1<-substring(Temp$Tidy_Modification,P_1+nchar("0.984016 ("),P_1+nchar("0.984016 (")+1 )
Modification_site_1

Modification_site_1<-as.numeric(Modification_site_1) #coerce to numeric vector
Modification_site_1

#extract 2nd modified sites
Modification_site_2<-substring(Temp$Tidy_Modification,P_2+nchar("0.984016 ("),P_2+nchar("0.984016 (")+1 )
Modification_site_2

Modification_site_2<-as.numeric(Modification_site_2) #coerce to numeric vector
Modification_site_2

#extract 3rd modified sites
Modification_site_3<-substring(Temp$Tidy_Modification,P_3+nchar("0.984016 ("),P_3+nchar("0.984016 (")+1 )
Modification_site_3

Modification_site_3<-as.numeric(Modification_site_3) #coerce to numeric vector
Modification_site_3

################insert annotation
#library(tidyr)
Temp<-Temp%>%
  dplyr::mutate(Modification_site_1,Modification_site_2,Modification_site_3)
Temp<-Temp%>%
  dplyr::mutate(Last_Mod_site=apply(Temp[,c("Modification_site_1","Modification_site_2","Modification_site_3")],1,max,na.rm=TRUE))%>% #add Modification_site_vector
  dplyr::mutate(PepSeq_pre_modification_1=str_sub(Temp$pepSeq_A,start=1,end=Temp$Modification_site_1),
                PepSeq_inter_modification_12=str_sub(Temp$pepSeq_A,start=Temp$Modification_site_1+1,end=Temp$Modification_site_2),
                PepSeq_inter_modification_23=str_sub(Temp$pepSeq_A,start=Modification_site_2+1,end=Temp$Modification_site_3),
                PepSeq_after_modification_3 =str_sub(Temp$pepSeq_A,start=Last_Mod_site+1,end=-1)) #add  PepSeq_pre_modification


#annotate deamidation,citrulinnation with #,@ respectively
#Temp$pepSeq_annotated<-paste0(Temp$PepSeq_pre_modification_1,"#",Temp$PepSeq_inter_modification_12,
#  "#", Temp$PepSeq_inter_modification_23,
#   "#",Temp$PepSeq_after_modification_3)

Temp$pepSeq_annotated<-ifelse(Temp$num_modification==1,
                              paste0(Temp$PepSeq_pre_modification_1,"#",Temp$PepSeq_after_modification_3),
                              ifelse(Temp$num_modification==2,
                                     paste0(Temp$PepSeq_pre_modification_1,"#",
                                            Temp$PepSeq_inter_modification_12,"#",
                                            Temp$PepSeq_after_modification_3),
                                     paste0(Temp$PepSeq_pre_modification_1,"#",
                                            Temp$PepSeq_inter_modification_12,"#",
                                            Temp$PepSeq_inter_modification_23,"#",
                                            Temp$PepSeq_after_modification_3))
)

#change annotating symbol
#Temp$pepSeq_annotated<-str_replace_all(Temp$pepSeq_annotated,"R#","R@")



#add pre and post amino acid
library(tibble)
Temp$pre<-str_sub(Temp$peptide_A,start = 1,end = 2)
Temp$post<-str_sub(Temp$peptide_A,start = -2,end = -1)
Temp$peptide_annotated_A<-paste0(Temp$pre,
                                 Temp$pepSeq_annotated,
                                 Temp$post)

#Inserting the annotated column
#example
#library(tibble)
#dataset <- data.frame(a = 1:5, b = 2:6, c=3:7)
#dataset<-add_column(dataset, d = 4:8, .after = 2)
#add_column(dataset, d = 4:8, .after = "b")
#add_column(dataset, d = 4:8, .before = "c")

Deamidated_ThrMod_NoMod<-add_column(Deamidated_ThrMod_NoMod,
                                    peptide_annotated_A=Temp$peptide_annotated_A,
                                    .after = "peptide_A")
#length(unique(Deamidated_ThrMod_NoMod$dataset_index)) #test
##



#append A_score############
library(dplyr)
Temp<-Deamidated_ThrMod_NoMod   #!update!#
Temp<-Temp%>%
  mutate(Job=1,spectrumID_number=str_extract(Temp$spectrumID,"[0-9]+"))%>% #spectrumID means Number.X of spectrum in .mgf file here; start from 0.
  select(Job,
         spectrumID_number,
         Charge,
         PrecursorMZ,DelM_PPM=ParentMassErrorPPM_A,
         Peptide=peptide_annotated_A,MSGFDB_SpecEValue=SpecEvalue_A) #format file for A_score searching
Temp$spectrumID_number<-as.numeric(Temp$spectrumID_number)
Temp<-Temp%>%
  mutate(Job=1,Scan=spectrumID_number+1)%>% #spectrumID means Number.X of spectrum in .mgf file here; start from 0.
  select(Job,
         Scan,
         Charge,
         PrecursorMZ,DelM_PPM,
         Peptide,MSGFDB_SpecEValue) #format file for A_score searching
str(Temp)

#
dir.create("Output_DScore_Benchmark")
write.table(Temp,file.path("Output_DScore_Benchmark","Deamidated_ThrMod_NoMod_fht.txt"),
            quote = F,sep ="\t",row.names = F ) #output fht file for AScore searching

#!STEP: Ascore searching
#!STEP: import Deamidated_ThrMod_NoMod_fht_ascore.txt by readr, sep="\t"

####Append ascore
Temp1<-Deamidated_ThrMod_NoMod   #update#
#Temp1$dataset_index<-as.numeric(Temp1$dataset_index)
colnames(Temp1)
#colnames(Temp1)[1]<-c("Scan")
##
Temp2<-Deamidated_ThrMod_NoMod_fht_ascore #UC_Filter_Same_Backbone_fht_ascore#%>%    #update#
colnames(Temp2)
Temp2$spectrumID<-paste0("index=",(Temp2$Scan-1))
Temp2<-Temp2%>%
  group_by(Job,spectrumID)%>%
  select(spectrumID,BestSequence,AScore,Job)%>%
  dplyr::slice(which.max(AScore)) #to get unique scans form unique datasets
#append
Temp<-left_join(Temp1,
                Temp2,
                by="spectrumID")
Deamidated_ThrMod_NoMod_Ascore<-Temp
####################


####Delta Score Filtration#####
##devide Deamidated_ThrMod_NoMod to same_backbone_AB and diff_backbone_AB
Same_Backbone_ThrMod_NoMod<-Deamidated_ThrMod_NoMod_Ascore%>%
  filter(Deamidated_ThrMod_NoMod_Ascore$pepSeq_A==Deamidated_ThrMod_NoMod_Ascore$pepSeq_B)
##
Diff_Backbone_ThrMod_NoMod<-Deamidated_ThrMod_NoMod_Ascore%>%
  filter(Deamidated_ThrMod_NoMod_Ascore$pepSeq_A!=Deamidated_ThrMod_NoMod_Ascore$pepSeq_B)

###MassError-Del_Score diagram
Temp<-Same_Backbone_ThrMod_NoMod%>%
  filter(msmsScore_A>10)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For same-seq spectra

#
Temp<-Diff_Backbone_ThrMod_NoMod%>%
  filter(msmsScore_A>0)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A) #For diff-seq spectra

##
#pick out confident deamidated/citrullinated same_seq PSMs
##
Filter_Same_Backbone_ThrMod_NoMod<-Same_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==1)%>%
  filter(Del_Score>0.5&
           msmsScore_A>10&
           abs(MassError_A)<5
         |
           Del_Score>4&
           msmsScore_A>10&
           abs(MassError_A)<10
  )
##
#pick out confident deamidated/citrullinated diff_seq PSMs
##
Filter_Diff_Backbone_ThrMod_NoMod<-Diff_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==1&
             Del_Score>2&
             msmsScore_A>=10&
             abs(MassError_A)<=5)

###
############Output tables#########################
###
####Simplest output, SameSeq spectra, 1-site citrullination
###

Temp<-Deamidated_ThrMod_NoMod_Ascore%>%
  filter(str_detect(BestSequence,"R#"), #extract 1-site deamidated/citrullinated PSMs
         str_count(BestSequence,"R#")==1, #extract citrullinated PSMs
         pepSeq_A==pepSeq_B,            #extract same-seq PSMs
         SpecEvalue_A<1e-10,
         abs(MassError_A)<5,
         Del_Score>0.5
  )
#select and rename columns
Temp<-Temp%>%
  select(spectrumID,`spectrum title`,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         BestSequence,AScore,
         idFile)
length(unique(Temp$BestSequence))
length(unique(Temp$Protein_A))

#output
dir.create("Output_DScore_Benchmark") #creat a new folder for results of a new dataset #update
write.csv(Temp,file.path("Output_DScore_Benchmark", "PSMs_SameSeq_1Site_Citrullination.csv")) #output #update folder name

###
####All outputs, SameSeq&DiffSeq spectra, 1-site deamidation/citrullination
###
###unique_pep_same_seq
###Format identified PSMs
Temp1<-Filter_Same_Backbone_ThrMod_NoMod
Temp1$BestSequence<-ifelse(!is.na(Temp1$BestSequence),
                           Temp1$BestSequence,
                           Temp1$peptide_annotated_A)
Temp1<-Temp1%>%
  group_by(BestSequence)%>%
  select(Peptide=BestSequence,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
         Del_Score,#Del_Raw_Score,
         Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
         PrecursorMZ,#Charge,
         AScore,
         #dataset_index,
         idFile,
         spectrumID,
         `spectrum title`)%>%
  dplyr::slice(which.max(Del_Score)) #Peptides with highest Del_Score as surrogate

Temp2<-Filter_Same_Backbone_ThrMod_NoMod%>%
  group_by(BestSequence)%>%
  select(Peptide=BestSequence)%>%
  group_by(Peptide)%>%
  dplyr::summarise(PSM_count=n()) #Counts of PSMs

Pep_Identified<-left_join(Temp1,Temp2,by="Peptide")
Temp<-Pep_Identified[str_detect(Pep_Identified$Protein,"XXX",negate=T),] #to remove PSMs in decoy database
Temp<-Temp%>%
  filter(!is.na(Peptide))%>%#   #get rid of peptide without bestsequence
  select(Peptide,PepSeq,Protein,Description,Del_Score,
         Spectrum_Evalue,MassErrorPPM,PrecursorMZ,PSM_count,AScore,
         #dataset_index,
         spectrumID,
         `spectrum title`)

#output
dir.create("Output_DScore_Benchmark") #creat a new folder for results of a new dataset #update
write.csv(Temp,file.path("Output_DScore_Benchmark", "02unique_pep_SameSeq.csv")) #output #update folder name
#
#stat
PSMs<-length(Temp$Peptide)
Pep<-length(unique(Temp$Peptide))
Protein<-length(unique(Temp$Protein))
N_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"N#")])
Q_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"Q#")])
R_citrullinated<-length(Temp$Peptide[str_detect(Temp$Peptide,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)

####unique_pep_diffseq
###Diff_sequence_1site
Temp1<-Filter_Diff_Backbone_ThrMod_NoMod
Temp1$BestSequence<-ifelse(!is.na(Temp1$BestSequence),
                           Temp1$BestSequence,
                           Temp1$peptide_annotated_A)
Temp1<-Temp1%>%
  filter(str_count(BestSequence,"#")==1)%>% #extract 1-site peptides
  group_by(BestSequence)%>%
  select(Peptide=BestSequence,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
         Del_Score,
         #Del_Raw_Score,
         Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
         PrecursorMZ,
         AScore,
         #dataset_index,
         idFile,
         spectrumID
         #Charge
  )%>%
  dplyr::slice(which.max(Del_Score)) #Peptides with highest Del_Score as surrogate

Temp2<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  group_by(BestSequence)%>%
  select(Peptide=BestSequence)%>%
  group_by(Peptide)%>%
  dplyr::summarise(PSM_count=n()) #Counts of PSMs

Pep_Identified<-left_join(Temp1,Temp2,by="Peptide")
Temp<-Pep_Identified[str_detect(Pep_Identified$Protein,"XXX",negate=T),]
Temp<-Temp%>%
  filter(!is.na(Peptide))

#output
write.csv(Temp,file.path("Output_DScore_Benchmark", "03unique_pep_DiffSeq_1site.csv")) #output #update folder name
#stat
PSMs<-length(Temp$Peptide)
Pep<-length(unique(Temp$Peptide))
Protein<-length(unique(Temp$Protein))
N_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"N#")])
Q_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"Q#")])
R_citrullinated<-length(Temp$Peptide[str_detect(Temp$Peptide,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)

#####PSMs_sameseq
Filter_Same_Backbone_ThrMod_NoMod<-Same_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==1)%>% #extract 1-site peptides
  filter(Del_Score>0.5&
           msmsScore_A>10&
           abs(MassError_A)<5
         |
           Del_Score>4&
           msmsScore_A>10&
           abs(MassError_A)<10
  )
length(unique(Filter_Same_Backbone_ThrMod_NoMod$peptide_A))

#output
#Temp<-Temp%>%
#select(Peptide=peptide_annotated_A,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
#Del_Score,Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
#PrecursorMZ)
Temp<-Filter_Same_Backbone_ThrMod_NoMod%>%
  select(spectrumID,`spectrum title`,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         BestSequence,AScore,
         idFile)
write.csv(Temp,file.path("Output_DScore_Benchmark", "04PSMs_SameSeq.csv")) #output #update folder name
#stat
Temp<-Filter_Same_Backbone_ThrMod_NoMod #update
PSMs<-length(Temp$peptide_annotated_A)
Pep<-length(unique(Temp$peptide_annotated_A))
Protein<-length(unique(Temp$protein_A))
N_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"N#")])
Q_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"Q#")])
R_citrullinated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)
#


########PSMs_diffseq
#all
colnames(Diff_Backbone_ThrMod_NoMod)
Filter_Diff_Backbone_ThrMod_NoMod<-Diff_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==1)%>% #extract 1-site peptides
  filter(Del_Score>2&
           msmsScore_A>10&
           abs(MassError_A)<5)

#output

Temp<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  select(spectrumID,`spectrum title`,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         BestSequence,AScore,
         idFile)
write.csv(Temp,file.path("Output_DScore_Benchmark", "05PSMs_DiffSeq.csv")) #output #update folder name

#stat
PSMs<-length(Temp$BestSequence)
Pep<-length(unique(Temp$BestSequence))
Protein<-length(unique(Temp$Protein_A))
N_deamidated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"N#")])
Q_deamidated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"Q#")])
R_citrullinated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)


######unique_pep_combined
Temp1<-Filter_Same_Backbone_ThrMod_NoMod
Temp1$BestSequence<-ifelse(!is.na(Temp1$BestSequence),
                           Temp1$BestSequence,
                           Temp1$peptide_annotated_A)
#
Temp2<-Filter_Diff_Backbone_ThrMod_NoMod
Temp2$BestSequence<-ifelse(!is.na(Temp2$BestSequence),
                           Temp2$BestSequence,
                           Temp2$peptide_annotated_A)
Temp2<-Temp2%>%
  filter(str_count(BestSequence,"#")==1)

#
Temp<-rbind(Temp1,
            Temp2)

Temp1<-Temp%>%
  group_by(BestSequence)%>%
  select(Peptide=BestSequence,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
         Del_Score,#Del_Raw_Score,
         Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
         PrecursorMZ,Charge,
         AScore,
         idFile,
         spectrumID)%>%
  dplyr::slice(which.max(Del_Score)) #Peptides with highest Del_Score as surrogate

Temp2<-Temp%>%
  group_by(BestSequence)%>%
  select(Peptide=BestSequence)%>%
  group_by(Peptide)%>%
  dplyr::summarise(PSM_count=n()) #Counts of PSMs

Pep_Identified<-left_join(Temp1,Temp2,by="Peptide")
Temp<-Pep_Identified[str_detect(Pep_Identified$Protein,"XXX",negate=T),] #to remove PSMs in decoy database
Temp<-Temp%>%
  filter(!is.na(Peptide))%>%
  select(Peptide,PepSeq,Protein,Description,Del_Score,
         Spectrum_Evalue,MassErrorPPM,PrecursorMZ,PSM_count,
         AScore,
         idFile,
        spectrumID)
#output
write.csv(Temp,file.path("Output_DScore_Benchmark", "01unique_pep_combined.csv")) #output #update folder name
#stat
PSMs<-length(Temp$Peptide)
Pep<-length(unique(Temp$Peptide))
Protein<-length(unique(Temp$Protein))
N_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"N#")])
Q_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"Q#")])
R_citrullinated<-length(Temp$Peptide[str_detect(Temp$Peptide,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)



###Benchmarking dual-search Delta Score strategy by a mannually inspected dataset##########
###

####link dataset with Manual inspection annotation in MCP paper https://doi.org/10.1074/mcp.RA118.000696>>>>>>
###
colnames(Deamidated_ThrMod_NoMod_Ascore)
#!Step: import Supplementary Table S1 by readr from paper: https://doi.org/10.1074/mcp.RA118.000696, title: Mining the Human Tissue Proteome for Protein Citrullination
#name the imported data as candidate_CitPep_MCP
colnames(candidate_CitPep_MCP)
#
Temp1<-Deamidated_ThrMod_NoMod_Ascore%>%
  select(spectrumID,
         `spectrum title`,
         PrecursorMZ,
         Peptide=BestSequence,
         pepSeq_A,
         pepSeq_B,
         msmsScore_A,
         Del_Score,
         MassError_A,
         SpecEvalue_A,
         AScore,
         protein_A,
         description_A)
#
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
Temp$PrecursorMZ<-round(Temp$PrecursorMZ,digits = 2)
Temp$instrument_PrecursorMZ<-round(Temp$instrument_PrecursorMZ,digits = 2)
Temp$MCP_PrecursorMZ<-round(Temp$MCP_PrecursorMZ,digits = 2)
#str(Temp)
##
Temp<-Temp%>%
  filter(Temp$PrecursorMZ==Temp$MCP_PrecursorMZ
         |
           Temp$instrument_PrecursorMZ==Temp$MCP_PrecursorMZ
  )

##

Temp<-Temp%>%
  distinct(spectrumID,.keep_all = T)

Deamidated_ThrMod_NoMod_Ascore_MI<-Temp

#output
dir.create("Output_DScore_Benchmark") #creat a new folder for results of a new dataset #update
write.csv(Temp,file.path("Output_DScore_Benchmark", "Deamidated_ThrMod_NoMod_Ascore_MI.csv")) #output #update folder name

##
###The Delta_Score of "valid" spectra
###
Temp<-Deamidated_ThrMod_NoMod_Ascore_MI%>%
  filter(`Manual inspection`=="valid")
plot(Temp$Del_Score,Temp$MassError_A)
dir.create("Output_DScore_Benchmark") #creat a new folder for results of a new dataset #update
write.csv(Temp,file.path("Output_DScore_Benchmark", "Del_Score_valid_spectra.csv")) #output #update folder name

###
########Cadidate PSMs for manual inspection#########>>>
###
##PSMs annotated as "valid" while their Delta_socre<0
Temp<-Deamidated_ThrMod_NoMod_Ascore_MI_cit%>%
  filter(Del_Score<0,
         `Manual inspection`=="valid"
  )
#Deamidated_ThrMod_NoMod_Ascore_MI_cit
Temp<-Deamidated_ThrMod_NoMod_Ascore_MI%>%
  filter(Del_Score>0,
         abs(MassError_A)>10,
         `Manual inspection`=="valid"
  )

###
###PSMs annotated as "uninspected" or "ambigurious" while thery are PSMs passing Delta Score criteria,
#that is 1-site citrullination, DeltaScore >2, abs(masserror)<5ppm & SpecEvalue<1e-10 (msmsscore>10).
###
colnames(Deamidated_ThrMod_NoMod_Ascore_MI_cit)
#SameSeq
Temp1<-Deamidated_ThrMod_NoMod_Ascore_MI_cit%>%
  filter(Deamidated_ThrMod_NoMod_Ascore_MI_cit$pepSeq_A==Deamidated_ThrMod_NoMod_Ascore_MI_cit$pepSeq_B,
         str_count(Peptide,"#")==1&
           Del_Score>0,
         abs(MassError_A)<5,
         msmsScore_A>10,
         `Manual inspection`!="valid" #208 PSMs
  )
#DiffSeq
Temp2<-Deamidated_ThrMod_NoMod_Ascore_MI_cit%>%
  filter(Deamidated_ThrMod_NoMod_Ascore_MI_cit$pepSeq_A!=Deamidated_ThrMod_NoMod_Ascore_MI_cit$pepSeq_B,
         str_count(Deamidated_ThrMod_NoMod_Ascore_MI_cit$Peptide,"#")==1,
         Del_Score>2,
         abs(MassError_A)<5,
         msmsScore_A>10,
         `Manual inspection`!="valid" #11 PSMs
  )


