library(tidyr)
library(dplyr)
library(data.table)
setwd("/work/PD_Proteomics/PrimaryDATA/")
df=as.data.frame(fread("V4_Proteomic_incident_PD_1112.csv")) #蛋白加covariate加AD
df=df[-1]
ProteinSuppInfo <- as.data.frame(fread("ProteinSuppInfo.csv")) 
Protein <- left_join(Protein,ProteinSuppInfo,by="eid")
Protein <- rename(Protein,"HLAA_SampAge"="HLA-A_SampAge","HLADRA_SampAge"="HLA-DRA_SampAge",
                  "HLAE_SampAge"="HLA-E_SampAge","HLAA"="HLA.A","HLAE"="HLA.E","HLADRA"="HLA.DRA",
                  "ERVV1"="ERVV.1","ERVV1_SampAge"="ERVV-1_SampAge")

#随机分组（这里按照YJ已经随机分好的组）#
## 随机切割 分2/3 1/2
{
  set.seed(1314)
  library(caret)
  Datasplit<-createDataPartition(y=df$eid,p=2/3,list=FALSE)  
  training_dataset<-df[Datasplit,]
  testing_dataset<-df[-Datasplit,]
  rm(df)
}

split_info <- read.csv("/work/PD_Proteomics/PrimaryDATA/Train_test_split_info.csv")
training_eid <- filter(split_info,split_id==0)
testing_eid <- filter(split_info,split_id==1)
training_dataset <- left_join(training_eid,Protein,by="eid")
training_dataset <- training_dataset[,c(1,3:5856)]
#write.csv(training_dataset,"/work/PD_Proteomics/PrimaryDATA/trainingYJ.csv")
testing_dataset <- left_join(testing_eid,Protein,by="eid")
testing_dataset <- testing_dataset[,c(1,3:5856)]
#write.csv(testing_dataset,"/work/PD_Proteomics/PrimaryDATA/testingYJ.csv")
training_dataset$Female <- as.factor(training_dataset$Female) 
training_dataset$smoking.status <- as.factor(training_dataset$smoking.status) 
training_dataset$Ethnicity <- as.factor(training_dataset$Ethnicity) 
training_dataset$alcohol.intake.frequency <- as.factor(training_dataset$alcohol.intake.frequency) 

for (i in 2:2921) {
  training_dataset[[i]]=(training_dataset[[i]]-mean(training_dataset[[i]],na.rm=T))/sd(training_dataset[[i]],na.rm = T)
}

## 分层
### 4 5Y内
eid5=Protein$eid[which(Protein$PD_y==1&Protein$follow_up<5*365)]# 更新后follow up是days
Protein2=Protein
Protein2$PD_y[which(!Protein2$eid%in%eid5)]=0
### 5 5Y外
Protein3=Protein
Protein3=Protein[-which(Protein$eid%in%c(eid5)),]

## 循环单因素cox(改一些变量名称)
library(survival)
Uni_glm_model=
  function(x){
    FML=as.formula(paste0("Surv(follow_up,PD_y==1)~",x,"+TDI+BMI+Female+Age+alcohol.intake.frequency+Ethnicity+smoking.status+fastingtime+season",
                          "+",paste0(x,"_SampAge"))) 
    Proname=x
    glm1=coxph(FML,data=Protein10) # cox回归方程
    GSum<- summary(glm1)
    coef=GSum$coefficients[,1]
    HR<- round(GSum$coefficients[,2],4)
    se=GSum$coefficients[,3]
    CI5=round(exp(coef-1.96*se),4)
    CI95=round(exp(coef+1.96*se),4)
    Pvalue<- GSum$coefficients[,5]
    Z=GSum$coefficients[,4]
    
    Uni_glm_model=cbind(Proname,HR,CI5,CI95,Pvalue,Z)#,,P_value_for_Schoenfeld_residuals,P_value_for_global_test
    Uni_glm_model=as.data.frame(Uni_glm_model)
    dimnames(Uni_glm_model)[[2]]=c("Protein","HR","CI5","CI95","P","Z")#,,,"P-value for Schoenfeld residuals","P-value for global test"
    
    return(Uni_glm_model) 
  } 

variable.names=colnames(Protein10)[2:2921] #ncol(predictors_data)
Uni_glm=vector(mode="list",length=length(variable.names))
Uni_glm=lapply(variable.names,Uni_glm_model)      

## 合并结果
result=data.frame()       
Uni_glm2=vector(mode="list",length=length(variable.names))
for (i in 1:length(variable.names)) {
  Uni_glm2[[i]]=Uni_glm[[i]][1,]
  result=rbind(result,Uni_glm2[[i]]) 
}

### Bonferroni矫正
result$P=as.numeric(result$P)
result=result[order(result$P),]
result$P_Bonferroni_adjust=p.adjust(
  result$P,  # P值列表
  method ="bonferroni"     # bonferroni校正的方法
)


full=result
write.csv(full,"/work/PD_Proteomics/Results/cox_resultse/train_full.csv")
less5=result
write.csv(less5,"/work/PD_Proteomics/Results/cox_results/train_less5.csv")
more5=result
write.csv(more5,"/work/PD_Proteomics/Results/cox_results/train_more5.csv")

train_full <- rename(train_full,"HR_full"="HR","CI5_full"="CI5","CI95_full"="CI95","P_full"="P",
                     "Z_full"="Z","P_Bonferroni_adjust_full"="P_Bonferroni_adjust")
train_less5Y <- rename(train_less5Y,"HR_less5"="HR","CI5_less5"="CI5","CI95_less5"="CI95","P_less5"="P",
                       "Z_less5"="Z","P_Bonferroni_adjust_less5"="P_Bonferroni_adjust")
train_more5Y <- rename(train_more5Y,"HR_more5"="HR","CI5_more5"="CI5","CI95_more5"="CI95","P_more5"="P",
                       "Z_more5"="Z","P_Bonferroni_adjust_more5"="P_Bonferroni_adjust")
coxresult_train <- left_join(train_full[,c(2:8)],train_less5Y[,c(2:8)],by="Protein")
coxresult_train <- left_join(coxresult_train,train_more5Y[,c(2:8)],by="Protein")

###internal validation###
pro_result <- as.data.frame(fread("/work/PD_Proteomics/Results/cox_results/coxresult_train"))
pro=data.frame("Assay"=unique(c(pro_result$Pro_code[which(pro_result$P_Bonferroni_full<0.05)],
                                pro_result$Pro_code[which(pro_result$P_Bonferroni_less5<0.05)],
                                pro_result$Pro_code[which(pro_result$P_Bonferroni_more5<0.05)])))# 需要富集的蛋白
name <- pro[,c(1)]
test_signif <- testing_dataset[,colnames(testing_dataset) %in% name]
test_signif <- cbind(testing_dataset$eid,test_signif)
test_signif <- rename(test_signif,"eid"="testing_dataset$eid")
Covariance <- testing_dataset[,c(2,2926:2935)]
test_signif <- left_join(Covariance,test_signif,by="eid")

## 分层
### 4 5Y内
eid10=test_signif$eid[which(test_signif$PD_y==1&test_signif$follow_up<5*365)]# 更新后follow up是days
Protein2=test_signif
Protein2$PD_y[which(!Protein2$eid%in%eid10)]=0
### 5 5Y外
Protein3=test_signif
Protein3=test_signif[-which(test_signif$eid%in%c(eid10)),]

## 循环单因素cox(改一些变量名称)
Uni_glm_model=
  function(x){
    FML=as.formula(paste0("Surv(follow_up,PD_y==1)~",x,"+TDI+BMI+Female+Age+alcohol.intake.frequency+Ethnicity+smoking.status")) # +APOE4 #+X21022  "+I(",x,"^2)",+Cholesterol+HDL_Cholesterol+SBP
    Proname=x
    glm1=coxph(FML,data=Protein3) # cox回归方程
    GSum<- summary(glm1)
    coef=GSum$coefficients[,1]
    HR<- round(GSum$coefficients[,2],4)
    se=GSum$coefficients[,3]
    CI5=round(exp(coef-1.96*se),4)
    CI95=round(exp(coef+1.96*se),4)
    Pvalue<- GSum$coefficients[,5]
    Z=GSum$coefficients[,4]
    
    Uni_glm_model=cbind(Proname,HR,CI5,CI95,Pvalue,Z)#,,P_value_for_Schoenfeld_residuals,P_value_for_global_test
    Uni_glm_model=as.data.frame(Uni_glm_model)
    dimnames(Uni_glm_model)[[2]]=c("Protein","HR","CI5","CI95","P","Z")#,,,"P-value for Schoenfeld residuals","P-value for global test"
    
    return(Uni_glm_model) 
  } 

variable.names=colnames(Protein3)[12:73] #ncol(predictors_data)
Uni_glm=vector(mode="list",length=length(variable.names))
Uni_glm=lapply(variable.names,Uni_glm_model)      

## 合并结果
result=data.frame()       
Uni_glm2=vector(mode="list",length=length(variable.names))
for (i in 1:length(variable.names)) {
  Uni_glm2[[i]]=Uni_glm[[i]][1,]
  result=rbind(result,Uni_glm2[[i]]) 
}

more5Y_test=result
write.csv(more5Y_test,"/work/PD_Proteomics/Results/random_New/testing_dataset_more5Y.csv")
less5Y_test=result
write.csv(less5Y_test,"/work/PD_Proteomics/Results/random_New/testing_dataset_less5Y.csv")
full_test=result
write.csv(full_test,"/work/PD_Proteomics/Results/random_New/testing_dataset_full.csv")


###功能分簇###
library(BiocManager)
library(ggplot2) 
library(stringr)
library(dplyr)
library(tidyr)
library(STRINGdb)
library(igraph)
library(data.table)

## load the protein ##
setwd("/work/PD_Proteomics/PrimaryDATA/")
pro_result=read.csv("/work/PD_Proteomics/Results/cox_results/random_yj/coxresult_train_YJ_1228.csv")
pro=data.frame("Assay"=unique(c(pro_result$Pro_code[which(pro_result$P_Bonferroni_full<0.05)],
                                pro_result$Pro_code[which(pro_result$P_Bonferroni_less5<0.05)],
                                pro_result$Pro_code[which(pro_result$P_Bonferroni_more5<0.05)])))# 需要富集的蛋白
symbol=read.csv("Olink_3072_panel_published_1127.csv") # 蛋白和gene名称不完全相同，可以按照这个nature提供的文件对应
pro=left_join(pro,symbol[,c(1,10)],by="Assay") # 写成标准symbol
rm(symbol)

pro_signif <- filter(pro_result,P_Bonferroni_full<0.05|P_Bonferroni_less5<0.05|P_Bonferroni_more5<0.05)
pro_signif$full <- 0
pro_signif$less5 <- 0
pro_signif$more5 <- 0
pro_signif$less5[which(pro_signif$P_Bonferroni_less5<0.05)] <- 1
pro_signif$more5[which(pro_signif$P_Bonferroni_more5<0.05)] <- 1
pro_signif$full[which(pro_signif$P_Bonferroni_full<0.05)] <- 1
write.csv(pro_signif,"/work/PD_Proteomics/Results/PPI/Pro_signif_group_1229.csv")


######################### 选择STRINGdb类型
string_db <- STRINGdb$new( version="12.0", #数据库版本。
                           species=9606,   #人9606，小鼠10090 
                           score_threshold=500, #蛋白互作的得分 默认400, 低150，高700，极高900 可以在400-900里选择连续数，比如650
                           input_directory="") 
######################### 获取STRING_id 
getOption('timeout')
#[1] 60
#设定timeout时间
options(timeout=100000)
##确认一下
getOption('timeout')
dat_map <- string_db$map(my_data_frame=pro, 
                         my_data_frame_id_col_names="HGNC.symbol", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = FALSE )
hits <- unique(dat_map$STRING_id)
#data_links <- string_db$get_interactions(hits)
######################### clustering分簇 
## iGraph clustering 互作网络分簇
#algorithm: fastgreedy(默认), walktrap, edge.betweenness
clustersList <- string_db$get_clusters(string_ids = hits ,
                                       algorithm  = "fastgreedy" ) 

Pro_clust=data.frame(matrix(NA, nrow = 4, ncol = 26)) # nrow:类别里最多有多少蛋白；ncol：有多少类别 
for (i in 1:26) {
  temp=clustersList[[i]]
  temp=unique(dat_map$HGNC.symbol[which(dat_map$STRING_id%in%temp)])
  length(temp)=4
  Pro_clust[i]=temp
} # 提取出每一类的gene symbol名称
write.csv(Pro_clust,"/work/PD_Proteomics/Results/PPI/Pro_cluster_38_600_0115.csv")


# 4 每个类别的富集，加上语义相似度计算
library(GOSemSim)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)

# 读入功能分簇文件生成的类别============
result=data.frame("ID"=NA,"Description"=NA,"GeneRatio"=NA,"BgRatio"=NA,"pvalue" =NA,
                  "p.adjust"=NA,"qvalue"=NA,"geneID"=NA,
                  "Count"=NA,"Redundant_with_term"=NA,"Category"=NA )


# 循环富集=================
for (i in c(1:5)) { # 1:category数量
  genes_full <- bitr(na.omit(Pro_clust[[i]]), fromType ="SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db) ## 转换差异表达基因
  genes=na.omit(genes_full) # 去掉未转化成功的gene_id
  genes=genes$ENTREZID  
  ego <- enrichGO(gene          = genes,
                  keyType = "ENTREZID", 
                  #universe      = geneList,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH", # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                  pvalueCutoff  =0.05, 
                  qvalueCutoff  = 0.05,
                  minGSSize     = 1,
                  readable      = TRUE)
  GO_data=data.frame(ego)
  GO_data=GO_data[order(GO_data$Count,decreasing = T),] #按count排序
  ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min) # 改成0.7
  GO_data2=data.frame(ego2)#去掉冗余的terms
  GO_data2=GO_data2[order(GO_data2$Count,decreasing = T),] #按count排序
  # 计算语义相似度:
  hsGO <- godata('org.Hs.eg.db', ont="BP") 
  SemSim=mgoSim(GO_data2$ID, GO_data2$ID, semData=hsGO, measure="Wang", combine=NULL)
  GO_data2$Redundant_with_term=NA
  for (j in 1:nrow(GO_data2)) {
    v=SemSim[j,][order(SemSim[j,],decreasing = T)]
    name=rownames(as.data.frame(which(v>=0.25&v!=1))) # 相似度评分0.25我都显示了，也可以调整
    name=GO_data2$Description[grep(paste0(name,collapse = "|"),GO_data2$ID)]
    GO_data2$Redundant_with_term[j]=paste0(name,collapse = "|")
  }
  GO_data2$GeneRatio=gsub("/","_", GO_data2$GeneRatio)
  GO_data2$BgRatio=gsub("/","_", GO_data2$BgRatio)
  GO_data2$Category=i
  result=rbind(result,GO_data2)
}
result <- result[-1,]
write.csv(result,"/work/PD_Proteomics/Results/PPI/Pro_enrich_38_650_1228.csv")

### Figure2 GO富集 ###
###按时间分组###
library(data.table)
library(dplyr)
library(ggplot2)
setwd("/work/PD_Proteomics/PrimaryDATA/")
full_result=as.data.frame(fread("/work/PD_Proteomics/Results/cox_results/random_yj/Derivation_full.csv"))
less5=as.data.frame(fread("/work/PD_Proteomics/Results/cox_results/random_yj/Derivation_5yrs.csv"))
more5=as.data.frame(fread("/work/PD_Proteomics/Results/cox_results/random_yj/Derivation_over5yrs.csv"))

pro=data.frame("Assay"=unique(c(full_result$Pro_code[which(full_result$pval_bfi<0.05)])))# 需要富集的蛋白
pro=data.frame("Assay"=unique(c(less5$Pro_code[which(less5$pval_bfi<0.05)])))# 需要富集的蛋白
pro=data.frame("Assay"=unique(c(more5$Pro_code[which(more5$pval_bfi<0.05)])))# 需要富集的蛋白

symbol=read.csv("Olink_3072_panel_published_1127.csv") # 蛋白和gene名称不完全相同，可以按照这个nature提供的文件对应
pro=left_join(pro,symbol[,c(1,10)],by="Assay") # 写成标准symbol
rm(symbol)

library(GOSemSim)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)

# 读入功能分簇文件生成的类别============
result=data.frame("ID"=NA,"Description"=NA,"GeneRatio"=NA,"BgRatio"=NA,"pvalue" =NA,
                  "p.adjust"=NA,"qvalue"=NA,"geneID"=NA,
                  "Count"=NA,"Redundant_with_term"=NA)
genes_full <- bitr(na.omit(pro$HGNC.symbol), fromType ="SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db)
genes=na.omit(genes_full) # 去掉未转化成功的gene_id
genes=genes$ENTREZID  
ego <- enrichGO(gene          = genes,
                keyType = "ENTREZID", 
                #universe      = geneList,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH", # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff  =0.05, 
                qvalueCutoff  = 0.05,
                minGSSize     = 1,
                readable      = TRUE)
GO_data=data.frame(ego)
GO_data=GO_data[order(GO_data$Count,decreasing = T),] #按count排序
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min) # 改成0.7
GO_data2=data.frame(ego2)#去掉冗余的terms
GO_data2=GO_data2[order(GO_data2$Count,decreasing = T),] #按count排序
# 计算语义相似度:
hsGO <- godata('org.Hs.eg.db', ont="BP") 
SemSim=mgoSim(GO_data2$ID, GO_data2$ID, semData=hsGO, measure="Wang", combine=NULL)
GO_data2$Redundant_with_term=NA
for (j in 1:nrow(GO_data2)) {
  v=SemSim[j,][order(SemSim[j,],decreasing = T)]
  name=rownames(as.data.frame(which(v>=0.25&v!=1))) # 相似度评分0.25我都显示了，也可以调整
  name=GO_data2$Description[grep(paste0(name,collapse = "|"),GO_data2$ID)]
  GO_data2$Redundant_with_term[j]=paste0(name,collapse = "|")
}
GO_data2$GeneRatio=gsub("/","_", GO_data2$GeneRatio)
GO_data2$BgRatio=gsub("/","_", GO_data2$BgRatio)
result=rbind(result,GO_data2)
result <- result[-1,]
write.csv(result,"/work/PD_Proteomics/Results/PPI/Full_enrich_36_600_1229.csv")
write.csv(result,"/work/PD_Proteomics/Results/PPI/Less5_enrich_17_600_1229.csv")
write.csv(result,"/work/PD_Proteomics/Results/PPI/More5_enrich_5_600_1229.csv")

### 可视化 ###
library(data.table)
library(dplyr)
result <- as.data.frame(fread("/work/PD_Proteomics/Results/PPI/Less5_enrich_17_600_1229.csv"))
result_order <- result[order(result$p.adjust),]
result_order <- rename(result_order,"padjust"="p.adjust")
result_signif <- filter(result_order,padjust < 0.05)

result_plot <- result_order[1:15,]
result_plot$term <- factor(result_plot$Description, levels = rev(result_plot$Description))
library(RColorBrewer)
display.brewer.all()
library(ggplot2)
pdf("/work/PD_Proteomics/Results/PPI/Barplot_less5(YJ)_1229.pdf",width = 10,height = 6)
ggplot(data = result_plot, aes(x = -log10(padjust), y = term, fill = -log10(padjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "-log10(adjusted P)", y = "") +
  scale_fill_distiller(palette = "Blues",direction = 1)+
  theme(
    panel.grid.major = element_line(colour=NA),
    panel.grid.minor = element_line(colour=NA),
    panel.background = element_rect(fill = "white"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.line = element_line(colour="black",linewidth = 0.6)
  )
dev.off()


###Trajectory###
library(data.table)
library(dplyr)
library(ggplot2)
df_loess <- as.data.frame(fread("/work/PD_Proteomics/Results/Trajectory/38pr/df_loess_38protein.csv"))
df_loess$Protein=as.factor(df_loess$Protein)
category <- as.data.frame(fread("/work/PD_Proteomics/Results/Trajectory/38pr/joined_clusters_protein.csv"))
library(data.table)
library(tidyr)
df_loess <- as.data.frame(fread("/work/PD_Proteomics/Results/Trajectory/38pr/df_loess_38protein.csv"))
df_loess$Protein=as.factor(df_loess$Protein)

df_heat=df_loess
max(df_heat$Estimate_loess,na.rm = T) # 0.4720419
min(df_heat$Estimate_loess,na.rm = T) # -0.9444244
df_heat$Estimate_loess[which(abs(df_heat$Estimate_loess)<0.3)]=NA #自定义一个“异常”的阈值，可以看看成图再选择
df_heat=pivot_wider(df_heat,names_from = Year,values_from = Estimate_loess)
for (i in ncol(df_heat):2) {
  df_heat=df_heat[order(abs(df_heat[[i]]),decreasing = T),]
} 
write.csv(df_heat,"/work/PD_Proteomics/Results/Trajectory/df_heat_0127.csv")
heat=df_heat[1:38,2:(ncol(df_heat))]
rownames(heat)=df_heat$Protein[1:38]

library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
cellwidth = 0.1
cellheight = 0.7
cn = dim(heat)[2]
rn = dim(heat)[1]
w=cellwidth*cn
h=cellheight*rn
col_fun1 = colorRamp2(c(-1.0, 0, 0.5), c("#0f86a9","#C6DBEF","#C0151B"))
col_fun1 = colorRamp2(c(-1.0, 0, 0.5), c("#0f86a9","white","#C0151B"))

setwd("/work/PD_Proteomics/Results/Trajectory/38pr/")
pdf("V1_full_heatmap_0121_015.pdf", width = 10, height = 15)
Heatmap(
  heat,  
  name = "r",  
  col = col_fun1,  
  width = unit(w, "cm"),  
  height = unit(h, "cm"),  
  rect_gp = gpar(col = "black", lwd = 0.05),  
  border_gp = gpar(col = "black", lty = 1, lwd = 1.2),  
  na_col = "white",  
  cluster_rows = F,  
  cluster_columns = F,
  row_title = NULL,  
  column_title = NULL, 
  column_names_gp = gpar(fontsize = 11),  
  row_names_gp = gpar(fontsize = 11), 
  heatmap_legend_param = list(  
    legend_height = unit(3, "cm"),
    grid_width = unit(0.4, "cm"),
    labels_gp = gpar(col = "gray20", fontsize = 8),
    title = "Zscore",
    direction = "horizontal"
  ))

dev.off()