######mod_search#########
###load package "rTANDEM"
library(rTANDEM)
##fasta file
fasta<-file.choose()
taxonomy <- rTTaxo(taxon="Homo sapiens",format="peptide",URL=fasta)
taxonomy
param <- rTParam()
param <- setParamValue(param, 'protein', 'taxon', value="Homo sapiens")
param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
##parameters
param <- setParamValue(param, 'list path', 'default parameters',value=system.file("extdata/default_input.xml", package="rTANDEM"))
param <- setParamOrbitrap(param)
param<-setParamValue(param, 'residue', 'potential modification mass',
                     value=c("15.994915@M,
                             0.984016@N,0.984016@N,0.984016@Q,0.984016@R"))#set dynamic modification in norm mode
# param<-setParamValue(param, 'residue', 'potential modification mass',
#                      value=NULL)#set dynamic modification in normal mode
param<-setParamValue(param, 'residue', 'modification mass',
                     value=NULL)#set fixed modification in normal mode
param<-setParamValue(param, 'spectrum', 'fragment monoisotopic mass error',
                     value=0.02)#set fixed modification in normal mode #set ms/ms error tolerance
print.rTParam(param)
##peak list
mgf<-file.choose()

param <- setParamValue(param, 'spectrum', 'path',
                       value=mgf)
##output
param<-setParamValue(param, 'output', 'results',
                     value="all")#set fixed modification in normal mode #set ms/ms error tolerance
param<-setParamValue(param, 'output', 'path hashing',
                     value="no")#set fixed modification in normal mode #set ms/ms error tolerance
param <- setParamValue(param, 'output', 'xsl path',
                       value=system.file("extdata/tandem-input-style.xsl", package="rTANDEM"))
dir.create("output_Xtandem")

param <- setParamValue(param, 'output', 'path',
                       value=file.path("output_Xtandem","MCP_cit_concatenated_mod.xml")) #update the name of output file
print.rTParam(param)

##database search
#result.path <- tandem(param)
tandem(param)

######ref_search#########
#parameter changed
param<-setParamValue(param, 'residue', 'potential modification mass',
                     value=c("15.994915@M")
                             )#set dynamic modification in norm mode
#output changed
param <- setParamValue(param, 'output', 'path',
                       value=file.path("output_Xtandem","MCP_cit_concatenated_ref.xml")) #update the name of output file
print.rTParam(param)
#Xtandem search
tandem(param)






