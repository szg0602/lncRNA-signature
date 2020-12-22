options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
library(BiocManager)
load(file = "LncRNA_exprSet.Rdata")
load(file = "result.Rdata")
diff_lncRNA_lab<-subset(res,(padj<0.01)&(abs(log2FoldChange))>2)

library(dplyr)
library(tidyr)
LncRNA_exprSet<-LncRNA_exprSet%>% 
  separate(gene_id,c("symbol","ensemble","genetype"),sep="\\|")%>%##separate the tcga_id
  select(- c(ensemble,genetype))%>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(LncRNA_exprSet))]))%>%
  arrange(desc(rowMean)) %>% #arrange by median of expression
  distinct(symbol,.keep_all = T) %>% # distinct  symbol
  select(-rowMean)


##differentially epressed lncRNAs#################################
index<-data.frame()
j<-1
for(i in diff_lncRNA_lab$symbol){
  location<-which(i==LncRNA_exprSet$symbol)
  index[j,1]<-location
  j<-j+1
  
}
diff_lncRNA_exprSet<-LncRNA_exprSet[index$V1,]
rownames(diff_lncRNA_exprSet)<-diff_lncRNA_exprSet$symbol
diff_lncRNA_exprSet<-select(diff_lncRNA_exprSet,-symbol)


mycounts <- diff_lncRNA_exprSet
keepGene=rowSums(edgeR::cpm(mycounts[-1])>0) >=2
table(keepGene);dim(mycounts)
mycounts<-mycounts[keepGene,]
mycounts_tumor<-mycounts[,group$sample=="tumor"]

########################################################

#2.2####################################################
load(file = "clinical_info.Rdata")  
cl_df<-as.data.frame(cl_df)

###Filter clinical information
clinical<-cl_df%>%
  select(c(bcr_patient_barcode,vital_status,gender,age_at_initial_pathologic_diagnosis,days_to_last_followup
           ,days_to_death,stage_event,race_list))

clinical$days_to_last_followup<-as.numeric(as.character(clinical$days_to_last_followup))
clinical$days_to_death<-as.numeric(as.character(clinical$days_to_death))

index<-which(!is.na(clinical$days_to_last_followup)&!is.na(clinical$days_to_death))
clinical[index,"days_to_last_followup"]<-NA

index<-which(is.na(clinical$days_to_last_followup)&is.na(clinical$days_to_death))
clinical<-clinical[-index,]
rownames(clinical)<-seq(1,nrow(clinical))

index<-which(is.na(clinical$days_to_last_followup))
clinical$days_to_last_followup[index]<-0

index<-which(is.na(clinical$days_to_death))
clinical$days_to_death[index]<-0

clinical$days<-clinical$days_to_last_followup+clinical$days_to_death

clinical<-clinical%>%
  select(-c(days_to_last_followup,days_to_death))

colnames(clinical)=c('ID','event','gender','age','stage','race',"days")

clinical$time=clinical$days/30

clinical$event_nu=ifelse(as.character(clinical$event)=='Alive',0,1)

clinical$age<-as.numeric(as.character(clinical$age))
clinical<-clinical[!(is.na(clinical$age)),]
clinical$age_group=ifelse(clinical$age>median(clinical$age),'older','younger')

library(stringr)
clinical$stage<-as.character(clinical$stage)
clinical$stage<-str_trim(clinical$stage,"both")


phe=clinical
head(phe)
phe$ID<-str_trim(phe$ID,"both")
phe$ID=toupper(phe$ID) 

#Match expression samples and clinical information samples by TCGA ID
phe=phe[match(substr(colnames(mycounts_tumor),1,12),phe$ID),]
list<-which(is.na(phe$ID))
mycounts_tumor<-mycounts_tumor[,-list]
phe<-phe[-list,]
head(phe)

##Remove genes whose expression levels are exactly the same in the sample
list<-apply(mycounts_tumor,1,function(gene){
  return(all(gene == gene[1]))
})
mycounts_tumor_final<-mycounts_tumor[!list,]
##################################################

#Univariate cox regression

library(survival)
library(survminer)
mySurv<-with(phe,Surv(time, event_nu))
xx<-log2(mycounts_tumor_final+1)
cox_results <-apply(xx , 1 , function(gene){
  gene<-as.numeric(gene)
  high_or_low=(ifelse(gene>median(gene),'high','alow') )
  survival_dat <- data.frame(high_or_low=high_or_low,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ high_or_low, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  # summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['high_or_lowhigh',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results_0.05<-as.data.frame(cox_results[cox_results[,4]<0.05,])
cox_results_0.05<-cbind(rownames(cox_results_0.05),cox_results_0.05)
colnames(cox_results_0.05)[1]<-"symbol"
cox_results_0.05_sort<-arrange(cox_results_0.05,p)

###################20 lncRNAs################################
index<-data.frame()
j<-1
ll<-c(1:12,14:20,23)
for(i in cox_results_0.05_sort$symbol[ll]){
  location<-which(i==rownames(mycounts_tumor_final))
  index[j,1]<-location
  j<-j+1
  
}
diff_lncRNA_exprSet_surv<-mycounts_tumor_final[index$V1,]

################################

#The data is divided into two groups randomly and equally
set.seed(23333)

k=sample(1:368,184)
diff_lncRNA_exprSet_surv_t<-diff_lncRNA_exprSet_surv[,k]
diff_lncRNA_exprSet_surv_v<-diff_lncRNA_exprSet_surv[,-k]

phe_t<-phe[k,]
phe_v<-phe[-k,]

#################################################################
#Further screening by rbsurv

library(rbsurv)
fit2 <- rbsurv(time=phe_t$time,status=phe_t$event_nu, x=log2(diff_lncRNA_exprSet_surv_t+1), method="efron", max.n.genes=12,n.iter = 1000,n.fold = 3,
              gene.ID = rownames(diff_lncRNA_exprSet_surv_t),seed = 123456)
save(fit2,file="rbsurv_fit2.Rdata")
load(file = "rbsurv_fit2.Rdata")
e=t(diff_lncRNA_exprSet_surv_t[fit2$gene.list[c(1:10)],])#top 10

gene<-c("AC0934261","LINC01446","AC0053071","OVAAL","NOVA1AS1","LINC02293","AL5920431","LINC01929","LINC00973","LINC00601")
colnames(e)<- gene

e<-log2(e+1)
########Multivariate cox regression##########################
dat=cbind(phe_t,e)
s=Surv(time, event_nu) ~AC0934261 + LINC01446 + AC0053071+OVAAL+
  NOVA1AS1 +LINC02293+AL5920431+LINC01929+LINC00973+LINC00601
model_coxph <- coxph(s, data = dat )
summary(model_coxph,data=dat)
########################################


############Model evaluation in training dataset################
library(ggfortify)
library(ggplot2)

#the risk score
fp_risk <- predict(model_coxph,dat,type="risk")
fp_dat<-data.frame(s=1:length(fp_risk),risk=as.numeric(sort(fp_risk )))
fp_dat$risk_group<-ifelse(fp_dat$risk>median(fp_dat$risk),"high","low")
ggplot(fp_dat,aes(x=s,y=risk))+geom_point(aes(col=risk_group))+
  labs(title="", x="Patien numbers",y="Risk Score")+
  geom_vline(xintercept = 92.5,lty=4,lwd=0.6,alpha=0.8)


sur_dat=data.frame(n=1:length(fp_risk),
                   t=phe[names(sort(fp_risk )),'time'] ,
                   e=phe[names(sort(fp_risk )),'event_nu']  ) 
sur_dat$event=ifelse(sur_dat$e==0,'alive','death')
ggplot(sur_dat,aes(x=n,y=t))+geom_point(aes(col=event))+
  labs(title="", x="Patien numbers",y="Survial Time (month)")+
  geom_vline(xintercept = 92.5,lty=4,lwd=0.6,alpha=0.8)

new_dat<-cbind(fp_dat,sur_dat)

sfit <- survfit(Surv(t, e)~risk_group, data=new_dat)
#K-M
survminer::ggsurvplot(sfit,alette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           xlab ="Time in months", 
           ggtheme = theme_gray())
#############################################

##########ROC curve##############
library(ROCR)
fp<-predict(model_coxph,dat)
re<-cbind(dat ,fp)

###########################
library(timeROC)
with(re,
     ROC <<- timeROC(T=time,
                     delta=event_nu,#survial status
                     marker=fp,
                     cause=1,
                     weighting="marginal",
                     times=c(36,60),
                     ROC = TRUE,
                     iid = TRUE)
)
# plot ROC curve
plot(ROC,time=36,col = "black",add =FALSE)
plot(ROC,time=60,col = "red",add = F)


##########for validate dataset#############

e=t(diff_lncRNA_exprSet_surv_v[ fit2$gene.list[c(1:10)],])
colnames(e)<- gene

dat_v=cbind(phe_v,e)
#the risk score
fp_risk_v <- predict(model_coxph,dat_v,type="risk")
fp_dat_v<-data.frame(s=1:length(fp_risk_v),risk=as.numeric(sort(fp_risk_v )))
fp_dat_v$risk_group<-ifelse(fp_dat_v$risk>median(fp_dat$risk),"high","low")
ggplot(fp_dat_v,aes(x=s,y=risk))+geom_point(aes(col=risk_group))+
  labs(title="", x="Patien numbers",y="Risk Score")+
  geom_vline(xintercept = 95.5,lty=4,lwd=0.6,alpha=0.8)
sur_dat_v=data.frame(n=1:length(fp_risk_v),
                     t=phe[names(sort(fp_risk_v )),'time'] ,
                     e=phe[names(sort(fp_risk_v )),'event_nu']  ) 
sur_dat_v$event=ifelse(sur_dat_v$e==0,'alive','death')
ggplot(sur_dat_v,aes(x=n,y=t))+geom_point(aes(col=event))+
  labs(title="", x="Patien numbers",y="Survial Time (month)")+
  geom_vline(xintercept = 95.5,lty=4,lwd=0.6,alpha=0.8)
new_dat_v<-cbind(fp_dat_v,sur_dat_v)

sfit <- survfit(Surv(t, e)~risk_group, data=new_dat_v)
#K-M
survminer::ggsurvplot(sfit,alette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           xlab ="Time in months", 
           ggtheme = theme_gray())

#####ROC of validate dataset
dat=cbind(phe_v,e)
fp<-predict(model_coxph,dat)
re<-cbind(dat ,fp)
###########################
library(timeROC)
with(re,
     ROC <<- timeROC(T=time,
                     delta=event_nu,#survial status
                     marker=fp,
                     cause=1,
                     weighting="marginal",
                     times=c(36,60),
                     ROC = TRUE,
                     iid = TRUE)
)
# plot roc curve
plot(ROC,time=36,col = "black",add =FALSE)

plot(ROC,time=60,col = "red",add = F)


########################################
########for entire dataset###########################
e=t(diff_lncRNA_exprSet_surv[fit2$gene.list[c(1:10)],])
colnames(e)<-gene
e<-log2(e+1)
dat_all=cbind(phe,e)
#the risk score
fp_risk_all <- predict(model_coxph,dat_all,type="risk")

fp_dat_all<-data.frame(s=1:length(fp_risk_all),risk=as.numeric(sort(fp_risk_all )))
fp_dat_all$risk_group<-ifelse(fp_dat_all$risk>median(fp_dat$risk),"high","low")
ggplot(fp_dat_all,aes(x=s,y=risk))+geom_point(aes(col=risk_group))+
  labs(title="", x="Patien numbers",y="Risk Score")+
  geom_vline(xintercept = 187.5,lty=4,lwd=0.6,alpha=0.8)

sur_dat_all=data.frame(n=1:length(fp_risk_all),
                       t=phe[names(sort(fp_risk_all )),'time'] ,
                       e=phe[names(sort(fp_risk_all )),'event_nu']  ) 
sur_dat_all$event=ifelse(sur_dat_all$e==0,'alive','death')
ggplot(sur_dat_all,aes(x=n,y=t))+geom_point(aes(col=event))+
  labs(title="", x="Patien numbers",y="Survial Time (month)")+
  geom_vline(xintercept = 187.5,lty=4,lwd=0.6,alpha=0.8)

new_dat_all<-cbind(fp_dat_all,sur_dat_all)

sfit <- survfit(Surv(t, e)~risk_group, data=new_dat_all)
#K-M
survminer::ggsurvplot(sfit,alette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           xlab ="Time in months", 
           ggtheme = theme_gray())

#################################################
##########entire dataset ROC
dat=cbind(phe,e)
fp<-predict(model_coxph,dat)
re<-cbind(dat ,fp)

###########################
library(timeROC)
with(re,
     ROC <<- timeROC(T=time, 
                     delta=event_nu,#survival status
                     marker=fp,
                     cause=1,
                     weighting="marginal",#
                     times=c(36,60),
                     ROC = TRUE,
                     iid = TRUE)
)
# plot ROC curve
plot(ROC,time=36,col = "black",add =FALSE)
plot(ROC,time=60,col = "red",add = F)



#####HeatMap and Volcano Plots 

#1.import lncRNA diferential analysis result
load(file = "LncRNA_exprSet.Rdata")
load(file = "result.Rdata")
diff_lncRNA_lab<-subset(res,(padj<0.01)&(abs(log2FoldChange))>2)

library(dplyr)
library(tidyr)
LncRNA_exprSet<-LncRNA_exprSet%>%
  separate(gene_id,c("symbol","ensemble","genetype"),sep="\\|")%>%
  select(- c(ensemble,genetype))%>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(LncRNA_exprSet))]))%>%
  arrange(desc(rowMean)) %>% #arrange by rowmean of expression
  distinct(symbol,.keep_all = T) %>% # Keep the first one
  select(-rowMean)


##the result#################################
index<-data.frame()
j<-1
for(i in diff_lncRNA_lab$symbol){
  location<-which(i==LncRNA_exprSet$symbol)
  index[j,1]<-location
  j<-j+1
  
}
diff_lncRNA_exprSet<-LncRNA_exprSet[index$V1,]

#################################################
#导import mRNAFinial diferential analysis result
load(file = "mRNA_exprSet.Rdata")
load(file = "mRNA_result.Rdata")
diff_mRNA_lab<-subset(mRNA_res,(padj<0.01)&(abs(log2FoldChange))>2)

mRNA_exprSet<-mRNA_exprSet%>%
  separate(gene_id,c("symbol","ensemble","genetype"),sep="\\|")%>%
  select(- c(ensemble,genetype))%>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(mRNA_exprSet))]))%>%
  arrange(desc(rowMean)) %>% #arrange by rowmean of expression
  distinct(symbol,.keep_all = T) %>% # Keep the first one
  select(-rowMean)


##the result#################################
index<-data.frame()
j<-1
for(i in diff_mRNA_lab$symbol){
  location<-which(i==mRNA_exprSet$symbol)
  index[j,1]<-location
  j<-j+1
  
}
diff_mRNA_exprSet<-mRNA_exprSet[index$V1,]
rownames(diff_mRNA_exprSet)<-diff_mRNA_exprSet$symbol
diff_mRNA_exprSet<-select(diff_mRNA_exprSet,-symbol)
#################################################
#import miRNA diferential analysis result
load(file = "miRNA_exprSet.Rdata")
load(file = "miRNA_result.Rdata")
diff_miRNA_lab<-subset(miRNA_res,(padj<0.01)&(abs(log2FoldChange))>2)

##the result#################################
index<-data.frame()
j<-1
for(i in diff_miRNA_lab$row){
  location<-which(i==miRNA_expr_df$gene_id)
  index[j,1]<-location
  j<-j+1
  
}
diff_miRNA_exprSet<-miRNA_expr_df[index$V1,]
rownames(diff_miRNA_exprSet)<-diff_miRNA_exprSet$gene_id
diff_miRNA_exprSet<-select(diff_miRNA_exprSet,-gene_id)
#################################################

########lncRNA HeatMap and Volcano Plots ########
####heatmap

if(! require("pheatmap")) install.packages("pheatmap")
#Make group information
res<-res%>%arrange(desc(log2FoldChange))
diffLab <- subset(res,(padj<0.01)&(abs(log2FoldChange))>2)

heatdata <- diff_lncRNA_exprSet[diffLab$symbol,]
annotation_col <- data.frame(group_list=group$sample)
rownames(annotation_col) <- colnames(heatdata)

group1<-group%>%arrange(sample)
annotation_col1 <- data.frame(group_list=group1$sample)
heatdata1<-heatdata%>%select(as.character(group1$TGCA_id))
rownames(annotation_col1) <- colnames(heatdata1)
colnames(annotation_col1)<-" "

bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
pheatmap(log2(heatdata1+1), #the data
         cluster_rows = T,
         cluster_cols = T,
          annotation_col =annotation_col1, 
          annotation_legend=T, 
         show_rownames = F,
         show_colnames = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/5*3),colorRampPalette(colors = c("white","#FF0000"))(length(bk)/5*2)),
         legend_breaks=seq(-5,5,1),## the legend of color
         breaks=bk,
         sepcolor="black",
         cellwidth   = 0.8,

)

library(ggplot2)
library(ggrepel)
library(dplyr)

data <- res

data$significant<-ifelse(data$padj<0.01 & abs(data$log2FoldChange) > 2,(ifelse(data$log2FoldChange>0,"UP(627)","DOWN(220)")),"NO(12883)")
list<-is.na(data$significant)
data$significant[list]<-"NO(12883)"
library(stringr)
data$symbol<-str_trim(data$symbol,"both")
genes<-c("AC093426.1","LINC01446","AC005307.1","OVAAL","NOVA1-AS1",
         "LINC02293","AL592043.1","LINC01929","LINC00973","LINC00601")
data$significant[which(data$symbol %in% genes)]<-"UP(627) target"

# Volcano Plots
ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("blue","black","red","green"))+
  labs(x="log2 (fold change)",y="-log10 (p-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.01),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(2,-2),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(#panel.border = element_blank(),
       panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
geom_text_repel(data=data[which(data$symbol %in% genes),], aes(label=symbol),col="black",alpha = 0.8, nudge_y = 25,
                segment.size = 0.5, 
                arrow = arrow(length=unit(0.01, "npc")),force = 100)


################################

########mRNA HeatMap and Volcano Plots########
####heatmap

#Make group information
mRNA_res<-mRNA_res%>%arrange(desc(log2FoldChange))
diffLab <- subset(mRNA_res,(padj<0.01)&(abs(log2FoldChange))>2)

heatdata <- diff_mRNA_exprSet[diffLab$symbol,]
annotation_col <- data.frame(group_list=group$sample)
rownames(annotation_col) <- colnames(heatdata)

group1<-group%>%arrange(sample)
annotation_col1 <- data.frame(group_list=group1$sample)
heatdata1<-heatdata%>%select(as.character(group1$TGCA_id))
rownames(annotation_col1) <- colnames(heatdata1)
colnames(annotation_col1)<-" "

bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
pheatmap(log2(heatdata1+1), #the data
         cluster_rows = T,
         cluster_cols = T,
         annotation_col =annotation_col1, 
         annotation_legend=T, 
         show_rownames = F,
         show_colnames = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/5*3),colorRampPalette(colors = c("white","#FF0000"))(length(bk)/5*2)),
         legend_breaks=seq(-5,5,1),#the legend of color
         breaks=bk,
         sepcolor="black",
         cellwidth   = 0.8,
)
library(ggplot2)
library(ggrepel)
library(dplyr)

data <- mRNA_res

data$significant<-ifelse(data$padj<0.01 & abs(data$log2FoldChange) > 2,(ifelse(data$log2FoldChange>0,"UP(656)","DOWN(816)")),"NO(18088)")
data$symbol<-str_trim(data$symbol,"both")
genes<-read.csv(file = "mRNAFinial.csv")

data$significant[which(data$symbol %in% genes$symbol)]<-"UP(656) target"

# Volcano Plots
ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("blue","black","red","green"))+
  labs(x="log2 (fold change)",y="-log10 (p-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.01),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(2,-2),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(#panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  geom_text_repel(data=data[which(data$symbol %in% genes$symbol),], aes(label=symbol),col="black",alpha = 0.8, nudge_y = 40,
                  segment.size = 0.5, 
                  arrow = arrow(length=unit(0.01, "npc")),force = 50)


################################

########miRNAHeatMap and Volcano Plots########
####heatmap
#Make group information
miRNA_res<-miRNA_res%>%arrange(desc(log2FoldChange))
diffLab <- subset(miRNA_res,(padj<0.01)&(abs(log2FoldChange))>2)

heatdata <- diff_miRNA_exprSet[diffLab$row,]
annotation_col <- data.frame(group_list=miRNA_group$sample)
rownames(annotation_col) <- colnames(heatdata)

miRNA_group1<-miRNA_group%>%arrange(sample)
annotation_col1 <- data.frame(group_list=miRNA_group1$sample)
heatdata1<-heatdata%>%select(as.character(miRNA_group1$TGCA_id))
rownames(annotation_col1) <- colnames(heatdata1)
colnames(annotation_col1)<-" "

bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
pheatmap(log2(heatdata1+1), #the data
         cluster_rows = T,
         cluster_cols = T,
         annotation_col =annotation_col1, 
         annotation_legend=T, 
         show_rownames = F,
         show_colnames = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/5*3),colorRampPalette(colors = c("white","#FF0000"))(length(bk)/5*2)),
         legend_breaks=seq(-5,5,1),# the legend of color
         breaks=bk,
         sepcolor="black",
         cellwidth   = 0.8,
)
library(ggplot2)
library(ggrepel)
library(dplyr)

data <- miRNA_res

data$significant<-ifelse(data$padj<0.01 & abs(data$log2FoldChange) > 2,(ifelse(data$log2FoldChange>0,"UP(44)","DOWN(20)")),"NO(789)")
data$row<-str_trim(data$row,"both")
genes<-read.csv(file = "miFinial.csv")

data$significant[which(data$row %in% tolower(genes$symbol))]<-"DOWN(20) target"

# Volcano Plots
ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("blue","green","black","red"))+
  labs( x="log2 (fold change)",y="-log10 (p-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.01),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(2,-2),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(#panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  geom_text_repel(data=data[which(data$row %in% tolower(genes$symbol)),], aes(label=row),col="black",alpha = 0.8, nudge_y = 30,nudge_x=4,
                  segment.size = 0.5, 
                  arrow = arrow(length=unit(0.01, "npc")),force = 50)

