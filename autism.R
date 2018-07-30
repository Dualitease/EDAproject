library(VariantAnnotation)
library(vcfR)
library(ggpubr)
library(ggplot2)
library(MASS)

#definitions of color format
whitecolor='pink1'
blackcolor='black'
hispaniccolor='tan3'
apicolor='yellow1'
mytol<-var(c(0,0,0,0.001))

#autism rate by group: white,black,hispanic,api (asian-pacific islander). It should be noted that all api samples are east asian, none pacific islander.
y<-c(12,10.2,7.9, 9.7) #retrieved from https://www.cdc.gov/mmwr/preview/mmwrhtml/ss6103a1.htm#Tab2



setwd('~/data/udacity/projects/1')

# Gene files obtained at https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/
# Method: search by gene name, click downloads, click download data for this region (change data to aggregates by population)
#         change file name to gene name, move file to working directory
#
# Selection: Genes selected for hypothesis autism and williams are opposites on a similar axis (overmentation->oversociality:systematize-empathize)
# these disorders overlap at ion metabolism, especially calcium, which is crucial in all vesicle release. Williams associated with too much
# calcium (over or pre release, proposed epigenetic-caused upregulated oxytocin receptors), autism associated with underfiring, especially of underdeveloped
# GABAergic networks and oxytocin activity. Autism is known to affect males more than females.
#
# relevant sites: https://www.nature.com/articles/mp201477
#                 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5667738/
#                 https://www.proteinatlas.org/ENSG00000188467-SLC24A5/tissue (slc expresses in male tissues, oxtr/avp extresses in female tissues)
#                 https://www.researchgate.net/publication/273148753_Oxytocin_Vasopressin_and_Williams_syndrome_Epigenetic_effects_on_Abnormal_Social_Behavior
#    
#                 http://www.wired.co.uk/video/how-music-can-help-with-strokes-autism-and-parkinsons-disease (music)
#                 https://pdfs.semanticscholar.org/58b6/0a01dd60cf45577fa39d23baa823f260608b.pdf (autism/music)
#                 https://en.wikipedia.org/wiki/Williams_syndrome (music: multiple references, they are talented)


oxtrfile<-read.vcfR('oxtr.popvcf') # oxytocin is a pro social hormone at behavioral level
avpr1afile<-read.vcfR('avpr1a.popvcf') # AVP intricately involved with social behavior, involved in osmotic pressure avoidance
                                       # in c elegans, as per Temple Grandin, autistics may prefer a hug machine (relates to sensory pressure)
slc24a5file<-read.vcfR('slc24a5.popvcf') #solute transporter (na,k,ca). At cellular level oxy/avp mediate solute concentration
syvn1file<-read.vcfR('syvn1.popvcf')     #ER misfolded protein gene. Protein blockage/vesicle inteference associated with dementia (and autism)

# Other genes of interest - cacna1c and vdr (cellular calcium), htr3a and slc6a4 (serotonin), CAMK2A, SNAP25 (intracellular calcium)
#                          caps2, INSR (BDNF release, relates to development of -among others- GABAergic networks, correlates to waist size)
# Note: plots do not necessarily represent n-terminus to c-terminus along x-axis. oxtr and avpr1a transcripts are read from
# right to left by position on chromosome.


#-----------------------------------OXTR GENE---------------------------------------------------
# 3 letter code is for distinct population, for example GBR is great britain, JPT is Tokyo, Japan
# vcf format contains actual data in gt. Character formating issue with , commas substituted to periods.

oxtr<-data.frame(rows=seq(1,length(oxtrfile@gt[,'GBR'])))
oxtr$GBR<-gsub(',','.',oxtrfile@gt[,'GBR'])
oxtr$FIN<-gsub(',','.',oxtrfile@gt[,'FIN'])
oxtr$CEU<-gsub(',','.',oxtrfile@gt[,'CEU'])

oxtr$YRI<-gsub(',','.',oxtrfile@gt[,'YRI'])
oxtr$MSL<-gsub(',','.',oxtrfile@gt[,'MSL'])
oxtr$GWD<-gsub(',','.',oxtrfile@gt[,'GWD'])

oxtr$CLM<-gsub(',','.',oxtrfile@gt[,'CLM'])
oxtr$PUR<-gsub(',','.',oxtrfile@gt[,'PUR'])
oxtr$PEL<-gsub(',','.',oxtrfile@gt[,'PEL'])

oxtr$CHB<-gsub(',','.',oxtrfile@gt[,'CHB'])
oxtr$CHS<-gsub(',','.',oxtrfile@gt[,'CHS'])
oxtr$JPT<-gsub(',','.',oxtrfile@gt[,'JPT'])

# creating numerator for division. numerator = count of alternate alleles per pop at SNP, aggregating by groups
numer<-as.numeric(unlist(strsplit(oxtr$GBR,':'))[seq(2,4*length(oxtr$GBR),4)])+
            as.numeric(unlist(strsplit(oxtr$FIN,':'))[seq(2,4*length(oxtr$FIN),4)])+
            as.numeric(unlist(strsplit(oxtr$CEU,':'))[seq(2,4*length(oxtr$CEU),4)])
# denom is quantity of polymorphism sites (all SNP) per gene, replicated to be vector of equal length as numer
denom<-(rep((as.numeric(unlist(strsplit(oxtr$GBR,':'))[1])+as.numeric(unlist(strsplit(oxtr$FIN,':'))[1])+
               as.numeric(unlist(strsplit(oxtr$CEU,':'))[1])),length(oxtr$GBR)))

oxtr$white<-numer/denom
##################################################################################


#rinse, repeat per group
numer<-as.numeric(unlist(strsplit(oxtr$YRI,':'))[seq(2,4*length(oxtr$YRI),4)])+
  as.numeric(unlist(strsplit(oxtr$MSL,':'))[seq(2,4*length(oxtr$MSL),4)])+
  as.numeric(unlist(strsplit(oxtr$GWD,':'))[seq(2,4*length(oxtr$GWD),4)])
denom<-(rep((as.numeric(unlist(strsplit(oxtr$YRI,':'))[1])+as.numeric(unlist(strsplit(oxtr$MSL,':'))[1])+
               as.numeric(unlist(strsplit(oxtr$GWD,':'))[1])),length(oxtr$YRI)))

oxtr$black<-numer/denom
##################################################################################


numer<-as.numeric(unlist(strsplit(oxtr$CLM,':'))[seq(2,4*length(oxtr$CLM),4)])+
  as.numeric(unlist(strsplit(oxtr$PUR,':'))[seq(2,4*length(oxtr$PUR),4)])+
  as.numeric(unlist(strsplit(oxtr$PEL,':'))[seq(2,4*length(oxtr$PEL),4)])
denom<-(rep((as.numeric(unlist(strsplit(oxtr$CLM,':'))[1])+as.numeric(unlist(strsplit(oxtr$PUR,':'))[1])+
               as.numeric(unlist(strsplit(oxtr$PEL,':'))[1])),length(oxtr$CLM)))

oxtr$hispanic<-numer/denom
##################################################################################




numer<-as.numeric(unlist(strsplit(oxtr$CHB,':'))[seq(2,4*length(oxtr$CHB),4)])+
  as.numeric(unlist(strsplit(oxtr$CHS,':'))[seq(2,4*length(oxtr$CHS),4)])+
  as.numeric(unlist(strsplit(oxtr$JPT,':'))[seq(2,4*length(oxtr$JPT),4)])
denom<-(rep((as.numeric(unlist(strsplit(oxtr$CHB,':'))[1])+as.numeric(unlist(strsplit(oxtr$CHS,':'))[1])+
               as.numeric(unlist(strsplit(oxtr$JPT,':'))[1])),length(oxtr$CHB)))

oxtr$api<-numer/denom
##################################################################################

#plot has polymorphism loci as x axis, and rate of variation on y axis. Smoothened by loess.

p1<-ggplot(data=oxtr, aes(x=seq(1,length(white))))+
  stat_smooth(method='loess',geom='line',aes(y=white),color=whitecolor)+
  stat_smooth(method='loess',geom='line',aes(y=black),color=blackcolor)+
  stat_smooth(method='loess',geom='line',aes(y=hispanic),color=hispaniccolor)+
  stat_smooth(method='loess',geom='line',aes(y=api),color=apicolor)+
  labs(y='alternate allele frequency', x='position on chromosome',title='OXTR gene')
#------------------------end OXTR GENE---------------------------------------------
##########################################################################################
##########################################################################################

#-----------------------------------avpr1a GENE---------------------------------------------------
avpr1a<-data.frame(rows=seq(1,length(avpr1afile@gt[,'GBR'])))
avpr1a$GBR<-gsub(',','.',avpr1afile@gt[,'GBR'])
avpr1a$FIN<-gsub(',','.',avpr1afile@gt[,'FIN'])
avpr1a$CEU<-gsub(',','.',avpr1afile@gt[,'CEU'])

avpr1a$YRI<-gsub(',','.',avpr1afile@gt[,'YRI'])
avpr1a$MSL<-gsub(',','.',avpr1afile@gt[,'MSL'])
avpr1a$GWD<-gsub(',','.',avpr1afile@gt[,'GWD'])

avpr1a$CLM<-gsub(',','.',avpr1afile@gt[,'CLM'])
avpr1a$PUR<-gsub(',','.',avpr1afile@gt[,'PUR'])
avpr1a$PEL<-gsub(',','.',avpr1afile@gt[,'PEL'])

avpr1a$CHB<-gsub(',','.',avpr1afile@gt[,'CHB'])
avpr1a$CHS<-gsub(',','.',avpr1afile@gt[,'CHS'])
avpr1a$JPT<-gsub(',','.',avpr1afile@gt[,'JPT'])

numer<-as.numeric(unlist(strsplit(avpr1a$GBR,':'))[seq(2,4*length(avpr1a$GBR),4)])+
  as.numeric(unlist(strsplit(avpr1a$FIN,':'))[seq(2,4*length(avpr1a$FIN),4)])+
  as.numeric(unlist(strsplit(avpr1a$CEU,':'))[seq(2,4*length(avpr1a$CEU),4)])
denom<-(rep((as.numeric(unlist(strsplit(avpr1a$GBR,':'))[1])+as.numeric(unlist(strsplit(avpr1a$FIN,':'))[1])+
               as.numeric(unlist(strsplit(avpr1a$CEU,':'))[1])),length(avpr1a$GBR)))

avpr1a$white<-numer/denom
##################################################################################



numer<-as.numeric(unlist(strsplit(avpr1a$YRI,':'))[seq(2,4*length(avpr1a$YRI),4)])+
  as.numeric(unlist(strsplit(avpr1a$MSL,':'))[seq(2,4*length(avpr1a$MSL),4)])+
  as.numeric(unlist(strsplit(avpr1a$GWD,':'))[seq(2,4*length(avpr1a$GWD),4)])
denom<-(rep((as.numeric(unlist(strsplit(avpr1a$YRI,':'))[1])+as.numeric(unlist(strsplit(avpr1a$MSL,':'))[1])+
               as.numeric(unlist(strsplit(avpr1a$GWD,':'))[1])),length(avpr1a$YRI)))

avpr1a$black<-numer/denom
##################################################################################


numer<-as.numeric(unlist(strsplit(avpr1a$CLM,':'))[seq(2,4*length(avpr1a$CLM),4)])+
  as.numeric(unlist(strsplit(avpr1a$PUR,':'))[seq(2,4*length(avpr1a$PUR),4)])+
  as.numeric(unlist(strsplit(avpr1a$PEL,':'))[seq(2,4*length(avpr1a$PEL),4)])
denom<-(rep((as.numeric(unlist(strsplit(avpr1a$CLM,':'))[1])+as.numeric(unlist(strsplit(avpr1a$PUR,':'))[1])+
               as.numeric(unlist(strsplit(avpr1a$PEL,':'))[1])),length(avpr1a$CLM)))

avpr1a$hispanic<-numer/denom
##################################################################################




numer<-as.numeric(unlist(strsplit(avpr1a$CHB,':'))[seq(2,4*length(avpr1a$CHB),4)])+
  as.numeric(unlist(strsplit(avpr1a$CHS,':'))[seq(2,4*length(avpr1a$CHS),4)])+
  as.numeric(unlist(strsplit(avpr1a$JPT,':'))[seq(2,4*length(avpr1a$JPT),4)])
denom<-(rep((as.numeric(unlist(strsplit(avpr1a$CHB,':'))[1])+as.numeric(unlist(strsplit(avpr1a$CHS,':'))[1])+
               as.numeric(unlist(strsplit(avpr1a$JPT,':'))[1])),length(avpr1a$CHB)))

avpr1a$api<-numer/denom
##################################################################################



p2<-ggplot(data=avpr1a, aes(x=seq(1,length(white))))+
  stat_smooth(method='loess',geom='line',aes(y=white),color=whitecolor)+
  stat_smooth(method='loess',geom='line',aes(y=black),color=blackcolor)+
  stat_smooth(method='loess',geom='line',aes(y=hispanic),color=hispaniccolor)+
  stat_smooth(method='loess',geom='line',aes(y=api),color=apicolor)+
  labs(y='alternate allele frequency', x='position on chromosome',title='avpr1a gene')

#------------------------end avpr1a GENE---------------------------------------------
##########################################################################################
##########################################################################################




#-----------------------------------slc24a5 GENE---------------------------------------------------
slc24a5<-data.frame(rows=seq(1,length(slc24a5file@gt[,'GBR'])))
slc24a5$GBR<-gsub(',','.',slc24a5file@gt[,'GBR'])
slc24a5$FIN<-gsub(',','.',slc24a5file@gt[,'FIN'])
slc24a5$CEU<-gsub(',','.',slc24a5file@gt[,'CEU'])

slc24a5$YRI<-gsub(',','.',slc24a5file@gt[,'YRI'])
slc24a5$MSL<-gsub(',','.',slc24a5file@gt[,'MSL'])
slc24a5$GWD<-gsub(',','.',slc24a5file@gt[,'GWD'])

slc24a5$CLM<-gsub(',','.',slc24a5file@gt[,'CLM'])
slc24a5$PUR<-gsub(',','.',slc24a5file@gt[,'PUR'])
slc24a5$PEL<-gsub(',','.',slc24a5file@gt[,'PEL'])

slc24a5$CHB<-gsub(',','.',slc24a5file@gt[,'CHB'])
slc24a5$CHS<-gsub(',','.',slc24a5file@gt[,'CHS'])
slc24a5$JPT<-gsub(',','.',slc24a5file@gt[,'JPT'])

numer<-as.numeric(unlist(strsplit(slc24a5$GBR,':'))[seq(2,4*length(slc24a5$GBR),4)])+
  as.numeric(unlist(strsplit(slc24a5$FIN,':'))[seq(2,4*length(slc24a5$FIN),4)])+
  as.numeric(unlist(strsplit(slc24a5$CEU,':'))[seq(2,4*length(slc24a5$CEU),4)])
denom<-(rep((as.numeric(unlist(strsplit(slc24a5$GBR,':'))[1])+as.numeric(unlist(strsplit(slc24a5$FIN,':'))[1])+
               as.numeric(unlist(strsplit(slc24a5$CEU,':'))[1])),length(slc24a5$GBR)))

slc24a5$white<-numer/denom
##################################################################################



numer<-as.numeric(unlist(strsplit(slc24a5$YRI,':'))[seq(2,4*length(slc24a5$YRI),4)])+
  as.numeric(unlist(strsplit(slc24a5$MSL,':'))[seq(2,4*length(slc24a5$MSL),4)])+
  as.numeric(unlist(strsplit(slc24a5$GWD,':'))[seq(2,4*length(slc24a5$GWD),4)])
denom<-(rep((as.numeric(unlist(strsplit(slc24a5$YRI,':'))[1])+as.numeric(unlist(strsplit(slc24a5$MSL,':'))[1])+
               as.numeric(unlist(strsplit(slc24a5$GWD,':'))[1])),length(slc24a5$YRI)))

slc24a5$black<-numer/denom
##################################################################################


numer<-as.numeric(unlist(strsplit(slc24a5$CLM,':'))[seq(2,4*length(slc24a5$CLM),4)])+
  as.numeric(unlist(strsplit(slc24a5$PUR,':'))[seq(2,4*length(slc24a5$PUR),4)])+
  as.numeric(unlist(strsplit(slc24a5$PEL,':'))[seq(2,4*length(slc24a5$PEL),4)])
denom<-(rep((as.numeric(unlist(strsplit(slc24a5$CLM,':'))[1])+as.numeric(unlist(strsplit(slc24a5$PUR,':'))[1])+
               as.numeric(unlist(strsplit(slc24a5$PEL,':'))[1])),length(slc24a5$CLM)))

slc24a5$hispanic<-numer/denom
##################################################################################




numer<-as.numeric(unlist(strsplit(slc24a5$CHB,':'))[seq(2,4*length(slc24a5$CHB),4)])+
  as.numeric(unlist(strsplit(slc24a5$CHS,':'))[seq(2,4*length(slc24a5$CHS),4)])+
  as.numeric(unlist(strsplit(slc24a5$JPT,':'))[seq(2,4*length(slc24a5$JPT),4)])
denom<-(rep((as.numeric(unlist(strsplit(slc24a5$CHB,':'))[1])+as.numeric(unlist(strsplit(slc24a5$CHS,':'))[1])+
               as.numeric(unlist(strsplit(slc24a5$JPT,':'))[1])),length(slc24a5$CHB)))

slc24a5$api<-numer/denom
##################################################################################



p3<-ggplot(data=slc24a5, aes(x=seq(1,length(white))))+
  stat_smooth(method='loess',geom='line',aes(y=white),color=whitecolor)+
  stat_smooth(method='loess',geom='line',aes(y=black),color=blackcolor)+
  stat_smooth(method='loess',geom='line',aes(y=hispanic),color=hispaniccolor)+
  stat_smooth(method='loess',geom='line',aes(y=api),color=apicolor)+
  labs(y='alternate allele frequency', x='position on chromosome',title='slc24a5 gene')
#------------------------end slc24a5 GENE---------------------------------------------
##########################################################################################
##########################################################################################
#-----------------------------------syvn1 GENE---------------------------------------------------
syvn1<-data.frame(rows=seq(1,length(syvn1file@gt[,'GBR'])))
syvn1$GBR<-gsub(',','.',syvn1file@gt[,'GBR'])
syvn1$FIN<-gsub(',','.',syvn1file@gt[,'FIN'])
syvn1$CEU<-gsub(',','.',syvn1file@gt[,'CEU'])

syvn1$YRI<-gsub(',','.',syvn1file@gt[,'YRI'])
syvn1$MSL<-gsub(',','.',syvn1file@gt[,'MSL'])
syvn1$GWD<-gsub(',','.',syvn1file@gt[,'GWD'])

syvn1$CLM<-gsub(',','.',syvn1file@gt[,'CLM'])
syvn1$PUR<-gsub(',','.',syvn1file@gt[,'PUR'])
syvn1$PEL<-gsub(',','.',syvn1file@gt[,'PEL'])

syvn1$CHB<-gsub(',','.',syvn1file@gt[,'CHB'])
syvn1$CHS<-gsub(',','.',syvn1file@gt[,'CHS'])
syvn1$JPT<-gsub(',','.',syvn1file@gt[,'JPT'])

numer<-as.numeric(unlist(strsplit(syvn1$GBR,':'))[seq(2,4*length(syvn1$GBR),4)])+
  as.numeric(unlist(strsplit(syvn1$FIN,':'))[seq(2,4*length(syvn1$FIN),4)])+
  as.numeric(unlist(strsplit(syvn1$CEU,':'))[seq(2,4*length(syvn1$CEU),4)])
denom<-(rep((as.numeric(unlist(strsplit(syvn1$GBR,':'))[1])+as.numeric(unlist(strsplit(syvn1$FIN,':'))[1])+
               as.numeric(unlist(strsplit(syvn1$CEU,':'))[1])),length(syvn1$GBR)))

syvn1$white<-numer/denom
##################################################################################



numer<-as.numeric(unlist(strsplit(syvn1$YRI,':'))[seq(2,4*length(syvn1$YRI),4)])+
  as.numeric(unlist(strsplit(syvn1$MSL,':'))[seq(2,4*length(syvn1$MSL),4)])+
  as.numeric(unlist(strsplit(syvn1$GWD,':'))[seq(2,4*length(syvn1$GWD),4)])
denom<-(rep((as.numeric(unlist(strsplit(syvn1$YRI,':'))[1])+as.numeric(unlist(strsplit(syvn1$MSL,':'))[1])+
               as.numeric(unlist(strsplit(syvn1$GWD,':'))[1])),length(syvn1$YRI)))

syvn1$black<-numer/denom
##################################################################################


numer<-as.numeric(unlist(strsplit(syvn1$CLM,':'))[seq(2,4*length(syvn1$CLM),4)])+
  as.numeric(unlist(strsplit(syvn1$PUR,':'))[seq(2,4*length(syvn1$PUR),4)])+
  as.numeric(unlist(strsplit(syvn1$PEL,':'))[seq(2,4*length(syvn1$PEL),4)])
denom<-(rep((as.numeric(unlist(strsplit(syvn1$CLM,':'))[1])+as.numeric(unlist(strsplit(syvn1$PUR,':'))[1])+
               as.numeric(unlist(strsplit(syvn1$PEL,':'))[1])),length(syvn1$CLM)))

syvn1$hispanic<-numer/denom
##################################################################################




numer<-as.numeric(unlist(strsplit(syvn1$CHB,':'))[seq(2,4*length(syvn1$CHB),4)])+
  as.numeric(unlist(strsplit(syvn1$CHS,':'))[seq(2,4*length(syvn1$CHS),4)])+
  as.numeric(unlist(strsplit(syvn1$JPT,':'))[seq(2,4*length(syvn1$JPT),4)])
denom<-(rep((as.numeric(unlist(strsplit(syvn1$CHB,':'))[1])+as.numeric(unlist(strsplit(syvn1$CHS,':'))[1])+
               as.numeric(unlist(strsplit(syvn1$JPT,':'))[1])),length(syvn1$CHB)))

syvn1$api<-numer/denom
##################################################################################



p4<-ggplot(data=syvn1, aes(x=seq(1,length(white))))+
  stat_smooth(method='loess',geom='line',aes(y=white),color=whitecolor)+
  stat_smooth(method='loess',geom='line',aes(y=black),color=blackcolor)+
  stat_smooth(method='loess',geom='line',aes(y=hispanic),color=hispaniccolor)+
  stat_smooth(method='loess',geom='line',aes(y=api),color=apicolor)+
  labs(y='alternate allele frequency', x='position on chromosome',title='syvn1 gene')
#------------------------end syvn1 GENE---------------------------------------------
##########################################################################################
##########################################################################################
p1<-p1+  theme(panel.background = element_rect(fill='lightskyblue',color='lightskyblue'),panel.grid.major=element_line(color='lightskyblue'),panel.grid.minor=element_line(color='lightskyblue'))
 
p2<-p2+  theme(panel.background = element_rect(fill='lightskyblue',color='lightskyblue'),panel.grid.major=element_line(color='lightskyblue'),panel.grid.minor=element_line(color='lightskyblue'))
  
p3<-p3+  theme(panel.background = element_rect(fill='lightskyblue',color='lightskyblue'),panel.grid.major=element_line(color='lightskyblue'),panel.grid.minor=element_line(color='lightskyblue'))
  
p4<-p4+  theme(panel.background = element_rect(fill='lightskyblue',color='lightskyblue'),panel.grid.major=element_line(color='lightskyblue'),panel.grid.minor=element_line(color='lightskyblue'))

ggarrange(p1,p3,p2,p4,ncol=2,nrow=2) #changed arrangement to put avp under oxy, since I think they compare better


#common autism SNP:rs53576, homozygous A is worst, homozygous G is best
# reference for this loci is A (autistic) per console command: oxtrfile@fix[,4][which(oxtrfile@fix[,3]=='rs53576')[1]]
#increased variance at this loci is negatively correlated to autism
p5<-p1 + geom_vline(show.legend=TRUE,xintercept=which(oxtrfile@fix[,3]=='rs53576')[1],color="darkolivegreen",linetype=2)
p5



###############################################
# Looking to compare variance on genes to autism rate by group membership
# Intended formula is autism rate~polymorphism 1+polymorphism 2+...+polymorphism 1943   (~ represents equals sign)
# This formula looks like a 4 row x 1 column answer (y), set equal 4 rows of sum((weighted computed differentials times) 1943 columns)
# lda() should return a 1 row x 1943 column of differentials.
#
# Known issue: lda() doesn't like columns having same data across any rows. It only wants unique variable inputs. Many of the loci have zero variants
#              per row.
# Obstacle:    Forumla won't accept a dataframe or matrix of 4x1943. It wants format of var1+var2+var3 (each representing a column of the 1943).
#              Not sure how to get the data into that format.


#aggregating the data into 4 rows x 1943 columns
combined<-data.frame(white=unlist(c(oxtr$white,avpr1a$white,slc24a5$white,syvn1$white)), #row 1
                     black=unlist(c(oxtr$black,avpr1a$black,slc24a5$black,syvn1$black)), #row 2
                     hispanic=unlist(c(oxtr$hispanic,avpr1a$hispanic,slc24a5$hispanic,syvn1$hispanic)), #row 3
                     api=unlist(c(oxtr$api,avpr1a$api,slc24a5$api,syvn1$api))) #row 4

mymatrix<-matrix(c(white=unlist(c(oxtr$white,avpr1a$white,slc24a5$white,syvn1$white)), #row 1
                   black=unlist(c(oxtr$black,avpr1a$black,slc24a5$black,syvn1$black)), #row 2
                   hispanic=unlist(c(oxtr$hispanic,avpr1a$hispanic,slc24a5$hispanic,syvn1$hispanic)), #row 3
                   api=unlist(c(oxtr$api,avpr1a$api,slc24a5$api,syvn1$api))),nrow=4,ncol=1943)

combined<-data.frame(mymatrix,row.names = c("white","black","hipanic","api"))

#preparing loci as column names in format of chr.position, eg chr3 position 87654321 = 3.87654321
oxtr$factors<-as.character(as.numeric(oxtrfile@fix[,1])+(as.numeric(oxtrfile@fix[,2])/1e7))
avpr1a$factors<-as.character(as.numeric(avpr1afile@fix[,1])+(as.numeric(avpr1afile@fix[,2])/1e8))
slc24a5$factors<-as.character(as.numeric(slc24a5file@fix[,1])+(as.numeric(slc24a5file@fix[,2])/1e8))
syvn1$factors<-as.character(as.numeric(syvn1file@fix[,1])+(as.numeric(syvn1file@fix[,2])/1e8))

mylist=list()
tempfactor=character()
#adding trailing zeroes
for (each in oxtr$factors){
  tempfactor<-each
  while (nchar(tempfactor)<9){
  tempfactor<-paste(tempfactor,"0",sep="")}
  mylist<-c(mylist,tempfactor)
}
oxtr$factors<-unlist(mylist)
mylist<-NULL

for (each in avpr1a$factors){
  tempfactor<-each
  while (nchar(tempfactor)<10){
    tempfactor<-paste(tempfactor,"0",sep="")}
  mylist<-c(mylist,tempfactor)
}
avpr1a$factors<-unlist(mylist)
mylist<-NULL

for (each in slc24a5$factors){
  tempfactor<-each
  while (nchar(tempfactor)<10){
    tempfactor<-paste(tempfactor,"0",sep="")}
  mylist<-c(mylist,tempfactor)
}
slc24a5$factors<-unlist(mylist)
mylist<-NULL

for (each in syvn1$factors){
  tempfactor<-each
  while (nchar(tempfactor)<10){
    tempfactor<-paste(tempfactor,"0",sep="")}
  mylist<-c(mylist,tempfactor)
}
syvn1$factors<-unlist(mylist)



colnames(combined)<-c(as.character(oxtr$factors),avpr1a$factors,slc24a5$factors,syvn1$factors)
combined<-combined[,!apply(combined, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))] #reference https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
#DATA prepared###############################################################################

combined<-combined/.001

ans<-lda(y~rownames(combined), data=combined)#,tol=mytol)



# combined<-data.frame(row=c("white","black","hispanic","api"))
# combined<-NULL
# combined<-rbind(combined,white=unlist(c(oxtr$white,avpr1a$white,slc24a5$white,syvn1$white)))
# combined<-rbind(combined,black=unlist(c(oxtr$black,avpr1a$black,slc24a5$black,syvn1$black)))
# combined<-rbind(combined,hispanic=unlist(c(oxtr$hispanic,avpr1a$hispanic,slc24a5$hispanic,syvn1$hispanic)))
# combined<-rbind(combined,api=unlist(c(oxtr$api,avpr1a$api,slc24a5$api,syvn1$api)))
# 
# combined<-combined[,1:4]
# #stuff that doesn't work
# ans<-lm(y~combined)
# 
# 
# #combined<-t(combined) 
# ans<-lda(combined[1,]+combined[2,]+combined[3,]+combined[4,],grouping=rows)
# 
#                     
                     