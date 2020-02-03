####Import and store ThrMod_Fhit mzid files in a list#############
#mzids <- file.choose() #choose searching result of msgf+, .mzid
#basename(mzids)
#mzids
# Vec<-c( #update
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_1_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_2_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_3_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_4_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_5_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_6_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_7_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_8_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_9_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_10_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_11_RZCol_270bar_55cm_022318_dta.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_12_RZCol_270bar_55cm_022318_dta.mzid"
# )



# UC_ThrMod_Fhit_ppm<-list()
# for (i in 1:12) {                #update
#   mzids<-Vec[i]
#   library("MSnID")
#   msnid <- MSnID(".")
#   msnid <- read_mzIDs(msnid, mzids) #tansfer .mzid file to S4 subject
#   #The following two function call create the new  numMisCleavages and numIrrCleabages columns in the MSnID object
#   msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
#   msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
#   #Trimming the data
#   #msnid <- apply_filter(msnid, "numIrregCleavages == 0")
#   #msnid <- apply_filter(msnid, "numMissCleavages <= 2")
#   msnid$ParentMassErrorPPM <- mass_measurement_error(msnid)
#   ##Filtering criteria
#   #-log10 transformed MS-GF+ Spectrum E-value,
#   #reflecting the goodness of match experimental and
#   #theoretical fragmentation patterns
#   msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
#   PSM_Ori<-msnid@psms          #retrive psms which data type is dataframe from msnid file.
#   #PSM_Mod<-msnid@psms
#
#   #######Manipulate exported dataframe file
#   library(dplyr)
#   library(stringr)
#   PSM_Ori<-as_tibble(PSM_Ori)
#   PSM_Ori<-PSM_Ori%>%
#     arrange(PSM_Ori$`scan number(s)`) #sorting by scan number
#   #PSM_Mod<-as_tibble(PSM_Mod)
#   #PSM_Ori_FHit<-PSM_Ori %>%
#   #distinct(spectrumID, .keep_all = TRUE)
#   #PSM_Ori_FHit
#   ####abs(DelM_ppm)<20
#   PSM_Ori_ppm<-PSM_Ori%>%
#     filter(abs(ParentMassErrorPPM)<=20)
#   #summary(PSM_Ori_ppm$ParentMassErrorPPM)
#
#   #generate first-hit ppm-filtered PSMs for PSM_Ori
#   UC_ThrMod_Fhit_ppm[[i]]<-PSM_Ori_ppm %>%
#     distinct(spectrumID, .keep_all = TRUE)
# }
#
# for (i in 1:12) { #update
#   print(unique(UC_ThrMod_Fhit_ppm[[i]]$idFile))
# }                                               #test to make sure every mzid file are stored


######Import and store NoMod_Fhit mzid files in a list########################
#mzids <- file.choose() #choose searching result of msgf+, .mzid
#basename(mzids)
#mzids
# Vec<-c( #update
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_1_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_2_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_3_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_4_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_5_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_6_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_7_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_8_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_9_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_10_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_11_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid",
#   "/Volumes/Xavier_Wang/Deamidation_R/data/Nanowell_files/mzid files/Nanowell_2DLC_islet_12_RZCol_270bar_55cm_022318_dta_without_deamidation.mzid"
# )
# UC_NoMod_Fhit<-list()
# for (i in 1:12) {    #update
#   mzids<-Vec[i]
#
#   library("MSnID")
#   msnid <- MSnID(".")
#   msnid <- read_mzIDs(msnid, mzids) #tansfer .mzid file to S4 subject
#   #The following two function call create the new  numMisCleavages and numIrrCleabages columns in the MSnID object
#   msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
#   msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
#   msnid$ParentMassErrorPPM <- mass_measurement_error(msnid) #add DelM_ppm
#   #-log10 transformed MS-GF+ Spectrum E-value,
#   #reflecting the goodness of match experimental and
#   #theoretical fragmentation patterns
#   msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
#   without_Ori<-msnid@psms          #retrive psms which data type is dataframe from msnid file.
#   #####tiding and sorting data
#   library(dplyr)
#   library(stringr)
#   without_Ori<-as_tibble(without_Ori)
#   without_Ori<-without_Ori%>%
#     arrange(without_Ori$`scan number(s)`) #sorting by scan number
#   without_Ori
#
#   UC_NoMod_Fhit[[i]]<-without_Ori%>%
#     distinct(spectrumID, .keep_all = TRUE)
#
# }
#
# for (i in 1:12) {   #update
#   print(unique(UC_NoMod_Fhit[[i]]$idFile))
# }                                               #test to make sure every mzid file are stored


##############Extract candidate deamidated/citrullinated PSMs######################
library(dplyr)
library(stringr)
#calculate Systamic Mass Shift
Temp<-rbind( #update
  UC_NoMod_Fhit[[1]],
  UC_NoMod_Fhit[[2]],
  UC_NoMod_Fhit[[3]],
  UC_NoMod_Fhit[[4]],
  UC_NoMod_Fhit[[5]],
  UC_NoMod_Fhit[[6]],
  UC_NoMod_Fhit[[7]],
  UC_NoMod_Fhit[[8]],
  UC_NoMod_Fhit[[9]],
  UC_NoMod_Fhit[[10]],
  UC_NoMod_Fhit[[11]],
  UC_NoMod_Fhit[[12]]
)
length(unique(Temp$idFile)) #test
colnames(Temp)
hist(Temp$ParentMassErrorPPM)

#calculation of systamic mass shift
library(stringr)
library(dplyr)
Temp<-Temp%>%
  filter(abs(ParentMassErrorPPM)<=20&
           Temp$`MS-GF:SpecEValue`<=1e-10)
SMS<-median(Temp$ParentMassErrorPPM)

#link datasets
Temp<-list()
for (i in 1:12) {                  #update
  In_A_B<-UC_NoMod_Fhit[[i]]%>%
    filter(`scan number(s)`%in% UC_ThrMod_Fhit_ppm[[i]]$`scan number(s)`)

  Combo_AB<-left_join(UC_ThrMod_Fhit_ppm[[i]],In_A_B,by="scan number(s)")

  Temp[[i]]<-Combo_AB %>%
    mutate(Del_Score=msmsScore.x-msmsScore.y,
           Del_Raw_Score=`MS-GF:RawScore.x`-`MS-GF:RawScore.y`,  #new 20181203
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
}

for (i in 1:12) {      #update
  print(unique(Temp[[i]]$dataset_index))
}  #to make sure every dataset is comboed

Temp<-rbind( #update
  Temp[[1]],
  Temp[[2]],
  Temp[[3]],
  Temp[[4]],
  Temp[[5]],
  Temp[[6]],
  Temp[[7]],
  Temp[[8]],
  Temp[[9]],
  Temp[[10]],
  Temp[[11]],
  Temp[[12]]
)
length(unique(Temp$dataset_index)) #test


#extracte candidate deamidated/citrullinated PSMs
Temp$modification_A<-gsub(pattern = "0.984015595000001",replacement = "0.984016",
                          Temp$modification_A) #for msgfplus_mzML method only
Deamidated_ThrMod_NoMod<-Temp%>%
  filter(str_detect(modification_A, "0.984016"))
colnames(Deamidated_ThrMod_NoMod)



#####go to annotation#############

Temp<-Deamidated_ThrMod_NoMod%>%
  select(`scan number(s)`,pepSeq_A,peptide_A,modification_A)
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
length(unique(Deamidated_ThrMod_NoMod$dataset_index)) #test
##########

#append A_score############
library(dplyr)
Temp<-Deamidated_ThrMod_NoMod   #!update!#
Temp<-Temp%>%
  select(Job=dataset_index,Scan =`scan number(s)`,Charge,
         PrecursorMZ,DelM_PPM=ParentMassErrorPPM_A,
         Peptide=peptide_annotated_A,MSGFDB_SpecEValue=SpecEvalue_A) #format file for A_score searching
View(Temp)

#
dir.create("Output_Nanowell_NoRefinary")
write.table(Temp,file.path("Output_Nanowell_NoRefinary","Deamidated_ThrMod_NoMod_fht.txt"),
            quote = F,sep ="\t",row.names = F ) #output fht file for AScore searching

#!STEP: Ascore searching
#!STEP: import Same_Backbone_ThrMod_NoMod_fht_ascore.txt by readr, sep="\t"


##
#Temp<-Filter_Same_Backbone
#colnames(Temp)[c(1,23)]<-c("Scan","Job")
#Add ascore
Temp1<-Deamidated_ThrMod_NoMod   #update#
#Temp1$dataset_index<-as.numeric(Temp1$dataset_index)
colnames(Temp1)
colnames(Temp1)[1]<-c("Scan")
Temp2<-Deamidated_ThrMod_NoMod_fht_ascore#UC_Filter_Same_Backbone_fht_ascore#%>%    #update#

Temp2<-Temp2%>%
  group_by(Job,Scan)%>%
  select(Scan,BestSequence,AScore,Job)%>%
  dplyr::slice(which.max(AScore)) #to get unique scans form unique datasets

#Temp2$Job<-as.numeric(Temp2$Job)

length(Temp1[Temp1$dataset_index==1,]$Scan)
length(Temp2[Temp2$Job==1,]$Scan)             #to see if Temp1 and Temp2 have same length

L<-list()   ##empty box to store
for (i in 1:12) {   #update number
  L[[i]]<-left_join(Temp1[Temp1$dataset_index==i,],
                    Temp2[Temp2$Job==i,],
                    by="Scan")
}
Deamidated_ThrMod_NoMod_Ascore<-rbind(L[[1]],L[[2]],L[[3]],L[[4]],L[[5]],L[[6]],L[[7]],L[[8]],L[[9]],L[[10]],
                                      L[[11]],L[[12]]
                                      # ,L[[13]],L[[14]],L[[15]],L[[16]],L[[17]],L[[18]],L[[19]],L[[20]],
                                      # L[[21]],L[[22]],L[[23]],L[[24]],L[[25]],L[[26]],L[[27]],L[[28]],L[[29]],L[[30]]
)

length(unique(Deamidated_ThrMod_NoMod_Ascore$dataset_index)) ##test
####################

####devide Deamidated_ThrMod_NoMod to same_backbone_AB and diff_backbone_AB#####

##
Same_Backbone_ThrMod_NoMod<-Deamidated_ThrMod_NoMod_Ascore%>%
  filter(Deamidated_ThrMod_NoMod_Ascore$pepSeq_A==Deamidated_ThrMod_NoMod_Ascore$pepSeq_B)

# Temp<-Same_Backbone_ThrMod_NoMod
#
# #stat
# PSMs<-length(Temp$peptide_annotated_A)
# Pep<-length(unique(Temp$peptide_annotated_A))
# Protein<-length(unique(Temp$protein_A))
# Stat<-cbind(PSMs,Pep,Protein)
# #output
# write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "PSMs_SameSeq_20ppm.csv")) #output


##
Diff_Backbone_ThrMod_NoMod<-Deamidated_ThrMod_NoMod_Ascore%>%
  filter(Deamidated_ThrMod_NoMod_Ascore$pepSeq_A!=Deamidated_ThrMod_NoMod_Ascore$pepSeq_B)

# Temp<-Deamidated_ThrMod_NoMod_Ascore
# #stat
# PSMs<-length(Temp$peptide_annotated_A)
# Pep<-length(unique(Temp$peptide_annotated_A))
# Protein<-length(unique(Temp$protein_A))
# Stat<-cbind(PSMs,Pep,Protein)
#
# #output
# write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "PSMs_DiffSeq_20ppm.csv")) #output


#####plot########
Temp<-Same_Backbone_ThrMod_NoMod%>%
  filter(msmsScore_A>=10)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A)

#
Temp<-Diff_Backbone_ThrMod_NoMod%>%
  filter(msmsScore_A>=10)
hist(Temp$MassError_A,
     breaks = 36)

plot(Temp$Del_Score,Temp$MassError_A)
#################


####Filtration#################
#same_Backbone
Filter_Same_Backbone_ThrMod_NoMod<-Same_Backbone_ThrMod_NoMod%>%
  filter(Del_Score>0&
           msmsScore_A>=10&
           abs(MassError_A)<=5
         |
           Del_Score>4&
           msmsScore_A>=10&
           abs(MassError_A)<=10
  )#%>%
#filter(Del_Raw_Score>13)
length(unique(Filter_Same_Backbone_ThrMod_NoMod$peptide_A))


# #stat
# Temp<-Filter_Same_Backbone_ThrMod_NoMod #update
# PSMs<-length(Temp$peptide_annotated_A)
# Pep<-length(unique(Temp$peptide_annotated_A))
# Protein<-length(unique(Temp$protein_A))
# N_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"N#")])
# Q_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"Q#")])
# R_citrullinated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"R#")])
# Stat<-cbind(PSMs,Pep,Protein,
#             N_deamidated,Q_deamidated,R_citrullinated)
# #output
# Temp<-Filter_Same_Backbone_ThrMod_NoMod
# #Temp<-Temp%>%
#        #select(Peptide=peptide_annotated_A,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
#        #Del_Score,Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
#        #PrecursorMZ)
# Temp<-Temp%>%
#   select(Scan,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
#          PSMs_A=peptide_annotated_A,PSMs_B=peptide_B,protein_A,protein_B,
#          SpecEvalue_A,SpecEvalue_B,Del_Score,BestSequence,AScore)
# write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "PSMs_identified_SameSeq.csv")) #output #update folder name
#
# #


# #####diff-backbone
#all
colnames(Diff_Backbone_ThrMod_NoMod)
Filter_Diff_Backbone_ThrMod_NoMod<-Diff_Backbone_ThrMod_NoMod%>%
  filter(Del_Score>2&
           msmsScore_A>=10&
           abs(MassError_A)<=5)
#
# #stat
# Temp<-Filter_Diff_Backbone_ThrMod_NoMod #update
# PSMs<-length(Temp$peptide_annotated_A)
# Pep<-length(unique(Temp$peptide_annotated_A))
# Protein<-length(unique(Temp$protein_A))
# N_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"N#")])
# Q_deamidated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"Q#")])
# R_citrullinated<-length(Temp$peptide_annotated_A[str_detect(Temp$peptide_annotated_A,"R#")])
# Stat<-cbind(PSMs,Pep,Protein,
#             N_deamidated,Q_deamidated,R_citrullinated)
# #output
# Temp<-Filter_Diff_Backbone_ThrMod_NoMod%>%
#   select(Peptide=peptide_annotated_A,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
#          Del_Score,Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
#          PrecursorMZ)
# write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "PSMs_identified_DiffSeq.csv")) #output #update folder name


#######################################

#unique_pep_same_seq#################
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
         dataset_index,
         idFile,
         Scan)%>%
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
         dataset_index,
         idFile,
         Scan)

#output

dir.create("Output_Nanowell_NoRefinary") #creat a new folder for results of a new dataset #update
write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "02unique_pep_SameSeq.csv")) #output #update folder name


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

####unique_pep_diffseq######################

###Diff_sequence_1site#####
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
         dataset_index,
         idFile,
         Scan
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
write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "03unique_pep_DiffSeq_1site.csv")) #output #update folder name


#stat
PSMs<-length(Temp$Peptide)
Pep<-length(unique(Temp$Peptide))
Protein<-length(unique(Temp$Protein))
N_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"N#")])
Q_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"Q#")])
R_citrullinated<-length(Temp$Peptide[str_detect(Temp$Peptide,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)

#####PSMs_sameseq#######################
Filter_Same_Backbone_ThrMod_NoMod<-Same_Backbone_ThrMod_NoMod%>%
  filter(Del_Score>0&
           msmsScore_A>=10&
           abs(MassError_A)<=5
         |
           Del_Score>4&
           msmsScore_A>=10&
           abs(MassError_A)<=10
  )#%>%
#filter(Del_Raw_Score>13)
length(unique(Filter_Same_Backbone_ThrMod_NoMod$peptide_A))

#output
#Temp<-Temp%>%
#select(Peptide=peptide_annotated_A,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
#Del_Score,Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
#PrecursorMZ)
Temp<-Filter_Same_Backbone_ThrMod_NoMod%>%
  select(Scan,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         dataset_index,BestSequence,AScore,
         dataset_index,
         idFile,
         Scan)
write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "04PSMs_SameSeq.csv")) #output #update folder name



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


########PSMs_diffseq##############
#all
colnames(Diff_Backbone_ThrMod_NoMod)
Filter_Diff_Backbone_ThrMod_NoMod<-Diff_Backbone_ThrMod_NoMod%>%
  filter(Del_Score>2&
           msmsScore_A>=10&
           abs(MassError_A)<=5)

#output

Temp<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  select(Scan,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         dataset_index,
         BestSequence,AScore,
         dataset_index,
         idFile,
         Scan)
write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "05PSMs_DiffSeq.csv")) #output #update folder name

#stat
Temp<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  select(Scan,Charge,PrecursorMZ,MassErrorPPM_A=MassError_A,
         Peptide_A=peptide_annotated_A,Peptide_B=peptide_B,
         Protein_A=protein_A,Protein_B=protein_B,
         SpecEvalue_A,SpecEvalue_B,Del_Score,
         dataset_index,
         BestSequence,AScore) #update

PSMs<-length(Temp$BestSequence)
Pep<-length(unique(Temp$BestSequence))
Protein<-length(unique(Temp$Protein_A))
N_deamidated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"N#")])
Q_deamidated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"Q#")])
R_citrullinated<-length(Temp$BestSequence[str_detect(Temp$BestSequence,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)


######unique_pep_combined################
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
         dataset_index,
         idFile,
         Scan)%>%
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
         dataset_index,
         idFile,
         Scan)

#output

write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "01unique_pep_combined.csv")) #output #update folder name


#stat
#Temp<-X01unique_pep_combined
PSMs<-length(Temp$Peptide)
Pep<-length(unique(Temp$Peptide))
Protein<-length(unique(Temp$Protein))
N_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"N#")])
Q_deamidated<-length(Temp$Peptide[str_detect(Temp$Peptide,"Q#")])
R_citrullinated<-length(Temp$Peptide[str_detect(Temp$Peptide,"R#")])
Stat<-cbind(PSMs,Pep,Protein,
            N_deamidated,Q_deamidated,R_citrullinated)


###unique_pep_diff_seq_many_sites_modified#######

#2-site
Temp1<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==2)%>% #extract 2-site peptides
  group_by(BestSequence)%>%
  select(Peptide=BestSequence,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
         Del_Score,
         #Del_Raw_Score,
         Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
         PrecursorMZ,
         AScore,
         dataset_index,
         idFile,
         Scan
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

length(unique(Temp$Protein))

#3-site

Temp1<-Filter_Diff_Backbone_ThrMod_NoMod%>%
  filter(str_count(BestSequence,"#")==3)%>% #extract 2-site peptides
  group_by(BestSequence)%>%
  select(Peptide=BestSequence,PepSeq=pepSeq_A,Protein=protein_A,Description=description_A,
         Del_Score,
         #Del_Raw_Score,
         Spectrum_Evalue=SpecEvalue_A,MassErrorPPM=MassError_A,
         PrecursorMZ,
         AScore,
         dataset_index,
         idFile,
         Scan
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

##stat
stat<-data.frame(PSM_1site=dim(filter(Filter_Diff_Backbone_ThrMod_NoMod,str_count(BestSequence,"#")==1))[1],
                 PSM_2site=dim(filter(Filter_Diff_Backbone_ThrMod_NoMod,str_count(BestSequence,"#")==2))[1],
                 PSM_3site=dim(filter(Filter_Diff_Backbone_ThrMod_NoMod,str_count(BestSequence,"#")==3))[1])

############################################################




location<-str_locate(Temp4$sequence,Temp4$PepSeq_DB_clean)
location<-as.data.frame(location)
Temp4$pep_DB_on_protein<-location$start

#ModSite_DB on protein
location_all<-str_extract_all(Temp4$PepSeq_DB_ModSite,"\\d+")

Temp4$P_1<-sapply(location_all, function(x)
  x[1])

Temp4$P_2<-sapply(location_all, function(x)
  x[2])

Temp4$P_3<-sapply(location_all, function(x)
  x[3])

Temp4[,c("pep_DB_on_protein","P_1","P_2","P_3")]<-apply(Temp4[,c("pep_DB_on_protein","P_1","P_2","P_3")],2,as.numeric)
str(Temp4) #check class

Temp4$P_1<-Temp4$pep_DB_on_protein+Temp4$P_1-1
Temp4$P_2<-Temp4$pep_DB_on_protein+Temp4$P_2-1
Temp4$P_3<-Temp4$pep_DB_on_protein+Temp4$P_3-1


Temp4$ModSite_DB_on_protein<-paste0(Temp4$P_1," ",Temp4$P_2," ",Temp4$P_3)
Temp4$ModSite_DB_on_protein<-str_remove_all(Temp4$ModSite_DB_on_protein,"NA")

###

#identified site in DB
str(Temp4)
Temp4$Identified_site_in_DB<-ifelse(  (Temp4$P_1-4)<Temp4$ModSite_on_protein&
                                        Temp4$ModSite_on_protein<(Temp4$P_1+4)|
                                        (Temp4$P_2-4)<Temp4$ModSite_on_protein&
                                        Temp4$ModSite_on_protein<(Temp4$P_2+4)|
                                        (Temp4$P_3-4)<Temp4$ModSite_on_protein&
                                        Temp4$ModSite_on_protein<(Temp4$P_3+4),
                                      paste0("Yes"),
                                      paste0("No"))
Temp4$Identified_site_in_DB[is.na(Temp4$Identified_site_in_DB)]<-"No" #replace NA with No

Temp4$seq_overlap<-"Yes"

#PepSeq_ModSite
location<-as.data.frame(str_locate(Temp4$Peptide,"[A-Z]#"))
AA<-str_sub(Temp4$Peptide,location$start,location$start)

location<-location$start-2

Temp4$PepSeq_ModSite<-paste0(AA,location)

#rearange column order
colnames(Temp4)
Temp<-Temp4%>%
  select(Peptide,
         PepSeq,
         Epitope_identified=PepSeq_ModSite,
         Antigen_iedb=PepSeq_DB_clean,
         Epitope_iedb=PepSeq_DB_ModSite,
         overlap,
         UniprotID,
         Description,
         `Organism Name`,
         Spectrum_Evalue,
         Del_Score,
         AScore,
         PSM_count,
         ProSeq=sequence,
         ProSite_epitope_identified=ModSite_on_protein,
         ProSite_epitope_P1_iedb=P_1,
         ProSite_epitope_P2_iedb=P_2,
         ProSite_epitope_P3_iedb=P_3,
         Epitope_identified_in_iedb=Identified_site_in_DB,
         seq_overlap,
         dataset_index,
         idFile,
         Scan,
         `Epitope IRI`
  )

write.csv(Temp,file.path("Output_Nanowell_NoRefinary", "12peptide_in_dea_cit_database.csv"))
peptide_in_dea_cit_database<-Temp

####


#####simple statistics

length(unique(unique_pep_combined$Protein)) #354
colnames(Protein_known_NQR)

unique(Protein_known_NQR$UniprotID) #129
unique(peptide_in_dea_cit_database$UniprotID) #39 proteins sequence overlap
unique(peptide_in_dea_cit_database$Peptide) #106

###