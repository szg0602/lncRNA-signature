options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
library(BiocManager)
#####################1 import data from file with downloaded from TCGA
rm(list = ls())
##gene document
metadata<-jsonlite::fromJSON("metadata.cart.2018-10-29.json")
if(!require(dplyr)){
  require(dplyr)
}

#Correspondence between files and samples
metadata<-metadata%>%
  select(c(file_name,associated_entities))

naid_df <- data.frame()
for (i in 1:dim(metadata)[1]){
  naid_df[i,1] <- substr(metadata$file_name[i],1,nchar(metadata$file_name[i])-3)
  naid_df[i,2] <- metadata$associated_entities[i][[1]]$entity_submitter_id
}
DIR<-getwd()
if(!dir.exists("files")){
  dir.create("files")
}
##copy count file 
DIR1<-paste0(DIR,"/counts")
for(dirnames in dir(DIR1)[1:length(dir(DIR1))]){
  dirnames<-paste0(DIR1,"/",dirnames)
  file<-list.files(dirnames,pattern = "*.counts.gz$")
  file.copy(paste0(dirnames,"/",file),"files")
}



##read count files
countFile<-list.files("files",pattern = "*.counts$",recursive = T)

if(!require(data.table))
{
  require(data.table)
}

location <- which(naid_df==countFile[1],arr.ind = TRUE)
TCGA_id <- as.character(naid_df[location[1],2])
expr_df<- fread(paste0("files/",countFile[1]))
names(expr_df) <- c("gene_id",TCGA_id)

for(i in 2:length(countFile)){
  
  location<-which(naid_df==countFile[i],arr.ind = TRUE)
  
  TCGA_id<-as.character(naid_df[location[1],2])
  dfnew<-fread(paste0("files/",countFile[i]))
  names(dfnew)<-c("gene_id",TCGA_id)
  expr_df<-inner_join(expr_df,dfnew,by="gene_id")
}
#####save expression

save(expr_df,file="exprSet.Rdata")

##import miRNA expression data
metadata<-jsonlite::fromJSON("metadata.cart.miRNA.2019-01-04.json")
metadata<-metadata%>%
  select(c(file_name,associated_entities))

naid_df <- data.frame()
for (i in 1:dim(metadata)[1]){
  naid_df[i,1] <- metadata$file_name[i]
  naid_df[i,2] <- metadata$associated_entities[i][[1]]$entity_submitter_id
}

DIR<-getwd()
if(!dir.exists("miRNA_files")){
  dir.create("miRNA_files")
}
DIR1<-paste0(DIR,"/miRNAcounts")
for(dirnames in dir(DIR1)[1:length(dir(DIR1))]){
  dirnames<-paste0(DIR1,"/",dirnames)
  file<-list.files(dirnames,pattern = "*mirnas.quantification.txt$")
  file.copy(paste0(dirnames,"/",file),"miRNA_files")
}
countFile<-list.files("miRNA_files",pattern = "*.txt$",recursive = T)
location <- which(naid_df==countFile[1],arr.ind = TRUE)
TCGA_id <- as.character(naid_df[location[1],2])
miRNA_expr_df<- fread(paste0("miRNA_files/",countFile[1]))
miRNA_expr_df<-miRNA_expr_df[,1:2]
names(miRNA_expr_df) <- c("gene_id",TCGA_id)

for(i in 2:length(countFile)){
  
  location<-which(naid_df==countFile[i],arr.ind = TRUE)
  
  TCGA_id<-as.character(naid_df[location[1],2])
  dfnew<-fread(paste0("miRNA_files/",countFile[i]))
  dfnew<-dfnew[,1:2]
  names(dfnew)<-c("gene_id",TCGA_id)
  miRNA_expr_df<-inner_join(miRNA_expr_df,dfnew,by="gene_id")
}

miRNA_group<-data.frame(names(miRNA_expr_df)[-1])
for(i in 1:length(miRNA_group[,1])){
  num<-as.numeric(substring(miRNA_group[i,1],14,15))
  if(num %in% seq(1,9)){
    miRNA_group[i,2]<-"tumor"
  }
  if(num %in% seq(10,29)){
    miRNA_group[i,2]<-"normal"
  }
}

names(miRNA_group)<-c("TGCA_id","sample")
save(miRNA_expr_df,miRNA_group,file = "miRNA_exprSet.Rdata")

#import miRNA expression data end


###############cinical data################

rm(list = ls())
Dir<-getwd()
if(!dir.exists("clinical_files")){
  dir.create("clinical_files")
}
Dir<-paste0(Dir,"/clinical")
for(dirnames in dir(Dir)[1:length(dir(Dir))]){
  dirnames<-paste0(Dir,"/",dirnames)
  file<-list.files(dirnames,pattern = "*.xml$")
  file.copy(paste0(dirnames,"/",file),"clinical_files")
}
Dir<-getwd()
Dir<-paste0(Dir,"/clinical_files")
all_fiels=list.files(path = Dir ,pattern='*.xml$',recursive=T)
all_fiels
length(unique(all_fiels))



library("XML")
library("methods")
cl = lapply(all_fiels
            , function(x){
              result <- xmlParse(file = file.path(Dir,x)) 
              rootnode <- xmlRoot(result)  
              xmldataframe <- xmlToDataFrame( rootnode[2] ) 
              return(t(xmldataframe))
            })


j<-1
k<-1
cl_cc<-list()
cl_cl<-list()

############Exclude invalid information#######
for(i in cl){
  if(length(i)==55){
    cl_cc[[j]]<-i
    j<-j+1
  }else {
    cl_cl[[k]]<-i
    k<-k+1
  }
 }

cl_df <- t(do.call(cbind,cl_cc))
cl_df1 <- t(do.call(cbind,cl_cl))
dim(cl_df)
cl_df[1:4,1:4]
save(cl_df,file = "clinical_info.Rdata")


#################2 handle data to get gene exprssion
#requeire package
if(!require("rtracklayer")){
  BiocManager::install("rtracklayer")
}
if(!require("SummarizedExperiment")){
  BiocManager::install("SummarizedExperiment")
}
###read gene annotation file to annotate gene
rm(list = ls())
load(file = "exprSet.Rdata")
#gtf <- rtracklayer::import('Homo_sapiens.GRCh38.90.chr.gtf')
#gtf_df <- as.data.frame(gtf)

#save(gtf_df,file = "gtf_df.Rdata")
###load gene anntotion data
load(file = "gtf_df.Rdata")

library(dplyr)
library(tidyr)

expr_df_nopoint <- expr_df %>% 
  separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>% 
  select(-drop)

save(expr_df_nopoint,file = "exprSet_nopoint.Rdata")

#############lncRNA##########

ncRNA <- c("sense_overlapping","lincRNA","3prime_overlapping_ncRNA",
           "antisense_RNA","sense_intronic","bidirectional_promoter_lncRNA","non_coding")


LncRNA_exprSet<-gtf_df%>%
  filter(type=="transcript",transcript_biotype%in% ncRNA) %>%
  select(c(gene_name,gene_id,transcript_biotype)) %>%
  distinct()%>%
  inner_join(expr_df_nopoint,by ="gene_id") %>% 
  unite(gene_id,gene_name,gene_id,transcript_biotype,sep = " | ")

LncRNA_exprSet<-na.omit(LncRNA_exprSet)

###divided into tumor or normal group by TCGA_id
group<-data.frame(names(LncRNA_exprSet)[-1])
for(i in 1:length(group[,1])){
  num<-as.numeric(substring(group[i,1],14,15))
  if(num %in% seq(1,9)){
    group[i,2]<-"tumor"
  }
  if(num %in% seq(10,29)){
    group[i,2]<-"normal"
  }
}

names(group)<-c("TGCA_id","sample")
save(LncRNA_exprSet,group,file = "LncRNA_exprSet.Rdata")


##############lncRNA end################

#####mRNA###############
mRNA_exprSet<-gtf_df%>%
  filter(type=="gene",gene_biotype=="protein_coding") %>%
  select(c(gene_name,gene_id,gene_biotype)) %>%
  inner_join(expr_df_nopoint,by ="gene_id") %>% 
  unite(gene_id,gene_name,gene_id,gene_biotype,sep = " | ")

mRNA_exprSet<-na.omit(mRNA_exprSet)

###divided into tumor or normal group by TCGA_id
group<-data.frame(names(mRNA_exprSet)[-1])
for(i in 1:length(group[,1])){
  num<-as.numeric(substring(group[i,1],14,15))
  if(num %in% seq(1,9)){
    group[i,2]<-"tumor"
  }
  if(num %in% seq(10,29)){
    group[i,2]<-"normal"
  }
}

names(group)<-c("TGCA_id","sample")
save(mRNA_exprSet,group,file = "mRNA_exprSet.Rdata")

#####mRNA end#######

#########start differential analysis###########
#require packages
if(!require("DESeq2")){
  install("DESeq2")
}
if(!require("edgeR")){
  install("edgeR")
}

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))

####load lncRNA expression data
rm(list = ls())
load(file ="LncRNA_exprSet.Rdata" )
library(dplyr)
library(tidyr)

group$sample<-as.factor(group$sample)

group%>%group_by(sample)%>%summarise(n())


mycounts <- LncRNA_exprSet1
##Exclude certain genes whose expression levels are low in all samples by edgeR
keepGene=rowSums(edgeR::cpm(mycounts[-1])>0) >=2
table(keepGene);dim(mycounts)
mycounts<-mycounts[keepGene,]

#differential analysis by DESeq2
dds<-DESeqDataSetFromMatrix(countData = mycounts,
                            colData = group2,
                            design = ~sample,
                            tidy=TRUE)

vsd<-vst(dds,blind = FALSE)
#Standardized data
LncRNA_exprSet_vst<-as.data.frame(assay(vsd))
save(LncRNA_exprSet_vst,file = "LncRNA_exprSet_vst.Rdata")
load(file = "LncRNA_exprSet_vst.Rdata")

###get result of differential analysis
dds<-DESeq(dds)


res<-results(dds,tidy = T)
res<-as.data.frame(res)
res<-res%>%
 separate(row,c("symbol","ensemble","genetype"),sep="\\|")%>%
  select(- c(ensemble,genetype))%>%
  arrange(abs(padj))%>%
  distinct(symbol,.keep_all = TRUE)%>%
  arrange(padj)
save(res,file = "result.Rdata")

####load mRNA data and start differential analysis
rm(list = ls())
load(file ="mRNA_exprSet.Rdata" )
library(dplyr)
library(tidyr)

group$sample<-as.factor(group$sample)

group%>%group_by(sample)%>%summarise(n())

##Exclude certain genes whose expression levels are low in all samples by edgeR
mycounts <- mRNA_exprSet
keepGene=rowSums(edgeR::cpm(mycounts[-1])>0) >=2
table(keepGene);dim(mycounts)

mycounts<-mycounts[keepGene,]


dds<-DESeqDataSetFromMatrix(countData = mycounts,
                            colData = group,
                            design = ~sample,
                            tidy=TRUE)
vsd<-vst(dds,blind = FALSE)
mRNA_exprSet_vst<-as.data.frame(assay(vsd))
save(mRNA_exprSet_vst,file = "mRNA_exprSet_vst.Rdata")
load(file = "mRNA_exprSet_vst.Rdata")
###get result of differential analysis
mRNA_dds<-DESeq(dds)
save(mRNA_dds,file = "mRNA_dds_DEseq.Rdata")
load(file = "mRNA_dds_DEseq.Rdata")
mRNA_res<-results(mRNA_dds,tidy = TRUE)
mRNA_res<-as.data.frame(mRNA_res)

index<-!is.na(mRNA_res$padj)
mRNA_res<-mRNA_res[index,]


mRNA_res<-mRNA_res%>%
  separate(row,c("symbol","ensemble","genetype"),sep="\\|")%>%
  select(- c(ensemble,genetype))%>%
  arrange(desc(abs(log2FoldChange)))%>%
  distinct(symbol,.keep_all = TRUE)%>%
  arrange(desc(log2FoldChange))

save(mRNA_res,file = "mRNA_result.Rdata")


###miRNA differential analysis
miRNA_group$sample<-as.factor(miRNA_group$sample)

miRNA_group%>%group_by(sample)%>%summarise(n())


miRNA_mycounts <- miRNA_expr_df
keepGene=rowSums(edgeR::cpm(miRNA_mycounts[-1])>0) >=2
table(keepGene);dim(miRNA_mycounts)
miRNA_mycounts<-miRNA_mycounts[keepGene,]

dds<-DESeqDataSetFromMatrix(countData = miRNA_mycounts,
                            colData = miRNA_group,
                            design = ~sample,
                            tidy=TRUE)
dds<-DESeq(dds)
save(dds,file = "miRNA_dds_DEseq.Rdata")
load(file = "miRNA_dds_DEseq.Rdata")
miRNA_res<-results(dds,tidy = TRUE)

save(miRNA_res,file = "miRNA_result.Rdata")
index<-!is.na(miRNA_res$padj)
miRNA_res<-miRNA_res[index,]
save(miRNA_res,file = "miRNA_result.Rdata")
###miRNA differential analysis end



