library(dplyr)
library(maftools)
library(paletteer)



####Load MAF files
Driver.genes <- c("RB1", "TP53", "NFE2L2", "CDKN2A", "PIK3CA")


TargetFile="../data/mutation_data/InputMutation.maf"
clinicalInput="../data/mutation_data/Clinical.tsv"
LungData = read.maf(maf = TargetFile,
                    clinicalData=clinicalInput,
                    vc_nonSyn=c('Missense',
                                'Nonsense',
                                'Nonstop',
                                'Frame_Shift',
                                "Inframe_INDEL",
                                "Splice_Site"))


ClinicalInfo<-data.frame(read.table(clinicalInput,header=T))
names(ClinicalInfo)[2]<-"Center"

# Set Color for Mutation Categories 
vc_cols = RColorBrewer::brewer.pal(n = 6, name = 'Paired')
names(vc_cols) = c(
  'Missense',
  'Nonsense',
  'Nonstop',
  'Frame_Shift',
  "Inframe_INDEL",
  "Splice_Site"
)

vc_cols[1]<-"#006DDBFF"
vc_cols[2]<-"#3B1B53FF"
vc_cols[3]<-"#79CC3DFF"
vc_cols[4]<-"#FF7F00FF"
vc_cols[5]<-"#117733FF"
vc_cols[6]<-"#7F00FFFF"



# Calculate 
ToCalculate.Frequency<-subsetMaf(maf = LungData, genes = Driver.genes, 
                                 mafObj = FALSE,
                                 fields="Center")[subsetMaf(maf = LungData, 
                                                            genes = Driver.genes, 
                                                            mafObj = FALSE,
                    fields="Center")$Variant_Classification%in%names(vc_cols),]

total_per_center <- ClinicalInfo %>%
  group_by(Center) %>%
  summarise(Total_Barcode = n_distinct(Tumor_Sample_Barcode))

df.frequency<-ToCalculate.Frequency %>%group_by(Center,Hugo_Symbol) %>% 
  summarise(Sample_Barcode_Frequency = n_distinct(Tumor_Sample_Barcode))
df.Percentage<-df.frequency%>% left_join(total_per_center,by="Center") %>%
  mutate(Percentage = (Sample_Barcode_Frequency / Total_Barcode) * 100)

Mut.Percentage.LUSC<-df.Percentage[which(df.Percentage$Center=="LUSC"),c(2,5)]
names(Mut.Percentage.LUSC)[2]<-"Percentage in LUSC"
Mut.Percentage.SCLC<-df.Percentage[which(df.Percentage$Center=="SCLC"),c(2,5)]
names(Mut.Percentage.SCLC)[2]<-"Percentage in SCLC"



# Sort Samples based on Mutation Order
# Taking All Sample IDs
SubsetLUSC<-subsetMaf(maf = LungData, genes = Driver.genes, mafObj = FALSE,fields="Center")
Total.LUSC<-ClinicalInfo[which(ClinicalInfo$Center=="LUSC"),]$Tumor_Sample_Barcode
Total.SCLC<-ClinicalInfo[which(ClinicalInfo$Center=="SCLC"),]$Tumor_Sample_Barcode
# For SCLC
# RB1 Mutated SCLC Samples
RB1.SCLC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="SCLC"&
                                                       SubsetLUSC$Hugo_Symbol=="RB1"&
                                                       SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# TP53 Mutated SCLC Samples
TP53.SCLC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="SCLC"&
                                                        SubsetLUSC$Hugo_Symbol=="TP53"&SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# RB1&TP53 Mutated SCLC Samples
Both.SCLC.MutID<-intersect(RB1.SCLC.MutID,TP53.SCLC.MutID)
# RB1 Mutated but TP53 Not mutated SCLC Samples
RB1.NonTP53.SCLC<-setdiff(RB1.SCLC.MutID,TP53.SCLC.MutID)
# TP53 Mutated but RB1 Not mutated SCLC Samples
OnlyTP53.SCLC<-setdiff(TP53.SCLC.MutID,RB1.SCLC.MutID)
# RB1 Not mutated SCLC Samples
RB1.SCLC.nonMutID<-setdiff(Total.SCLC,RB1.SCLC.MutID)
# RB1&TP53 Not mutated SCLC Samples
TP53.RB1.SCLC.nonMutID<-setdiff(RB1.SCLC.nonMutID,TP53.SCLC.MutID)

# For LUSC
# RB1 Mutated LUSC Samples
RB1.LUSC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="LUSC"&
                                                       SubsetLUSC$Hugo_Symbol=="RB1"&SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# TP53 Mutated LUSC Samples
TP53.LUSC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="LUSC"&
                                                        SubsetLUSC$Hugo_Symbol=="TP53"&
                                                        SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# NFE2L2 Mutated LUSC Samples
NFE2L2.LUSC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="LUSC"&
                                                          SubsetLUSC$Hugo_Symbol=="NFE2L2"&SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# NOTCH1 Mutated LUSC Samples
NOTCH1.LUSC.MutID<-as.character(unique(SubsetLUSC[which(SubsetLUSC$Center=="LUSC"&
                                                          SubsetLUSC$Hugo_Symbol=="NOTCH1"&
                                                          SubsetLUSC$Variant_Classification%in%names(vc_cols)),]$Tumor_Sample_Barcode))
# RB1&TP53 Mutated LUSC Samples
Both.LUSC.MutID<-intersect(RB1.LUSC.MutID,TP53.LUSC.MutID)
# RB1 Mutated but TP53 Not mutated LUSC Samples
RB1.NonTP53.LUSC<-setdiff(RB1.LUSC.MutID,TP53.LUSC.MutID)
# Only TP53 mutated LUSC Samples
OnlyTP53.LUSC<-setdiff(TP53.LUSC.MutID,RB1.LUSC.MutID)
# RB1 Not mutated LUSC Samples
RB1.LUSC.nonMutID<-setdiff(Total.LUSC,RB1.LUSC.MutID)
# RB1&TP53 Not mutated LUSC Samples
TP53.RB1.LUSC.nonMutID<-setdiff(RB1.LUSC.nonMutID,TP53.LUSC.MutID)
# Only NFE2L2 mutated LUSC Samples
OnlyNFE2L2.MutId<-intersect(TP53.RB1.LUSC.nonMutID,NFE2L2.LUSC.MutID)
# Notch1 Mutated among RB1 Not mutated LUSC Samples
Notch1.Among.Non.RB1.TP53<-intersect(TP53.RB1.LUSC.nonMutID,NOTCH1.LUSC.MutID)
# Only Notch1 Mutated LUSC Samples
OnlyNotch1.MutID<-setdiff(Notch1.Among.Non.RB1.TP53,OnlyNFE2L2.MutId)
# NFE2L2&TP53&RB1 Not mutated LUSC Samples
NFE2L2.TP53.RB1Non.MutID<-setdiff(TP53.RB1.LUSC.nonMutID,OnlyNFE2L2.MutId)
# Notech1&NFE2L2&TP53&RB1 Not mutated LUSC Samples
Notch1.NFE2L2.TP53.RB1Non.MutID<-setdiff(NFE2L2.TP53.RB1Non.MutID,OnlyNotch1.MutID)
# Assemble
SampleOrder=c(Both.SCLC.MutID,RB1.NonTP53.SCLC,OnlyTP53.SCLC,TP53.RB1.SCLC.nonMutID, #Total should be 109
              Both.LUSC.MutID,RB1.NonTP53.LUSC,OnlyTP53.LUSC,
              OnlyNFE2L2.MutId,OnlyNotch1.MutID,Notch1.NFE2L2.TP53.RB1Non.MutID) #Total should be 47



# Plot
pdf("../figures/oncoplot.pdf")
oncoplot(maf = LungData, colors = vc_cols, genes=Driver.genes,clinicalFeatures = c('Cancer_Type'),
         sortByAnnotation = F,draw_titv = F,showTitle=F,
         leftBarData = Mut.Percentage.SCLC,
         leftBarLims = c(0, 100),
         rightBarData = Mut.Percentage.LUSC,
         rightBarLims = c(0,100),GeneOrderSort = F,
         keepGeneOrder=T,removeNonMutated = F,sampleOrder=SampleOrder,
         sepwd_samples = 0,sepwd_genes = 0,legendFontSize = 1,annotationFontSize =1,fontSize = 0.7,anno_height=0.5
)
dev.off()