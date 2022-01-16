ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/Pan-cancer-graphs"

setwd(ROOT_DIR)

data <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_data.txt", header=T, sep='\t', as.is=T)

head(data$Genome_doublings)
TYPES<-unique(data$Type)

library("stringr")

#subsetting data based on WGD status
CA <- data
CA_WGD0<-CA[CA$Genome_doublings==0,]
CA_WGD1<-CA[CA$Genome_doublings==1,]
CA_WGD2<-CA[CA$Genome_doublings==2,]
CA_WGD<-CA[CA$Genome_doublings>0,]

#all tumor types with >=20 events in each group
g_20_types <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC")


  
WGD_fractions <- c()
WGD0_fractions <- c()

for (t in 1:length(g_20_types)){

  library("MASS")
  
  g_20_types = c("ACC","BLCA","CESC", "BRCA", "COAD", "ESCA","GBM" ,"HNSC", "KIRC","LGG" ,"LIHC", "LUAD", "LUSC","OV", "PAAD", "PCPG", "PRAD", "READ","SARC", "SKCM", "STAD", "UCEC")
  
  
  df <- data.frame(frac_top3 = c(), wgd_status = as.character(c()), type = as.character(c()))
  
  
  ty <- g_20_types[t]
  
  CA<-data[data$Type==ty,]
  CA_WGD0<-CA[CA$Genome_doublings==0,]
  CA_WGD1<-CA[CA$Genome_doublings==1,]
  CA_WGD2<-CA[CA$Genome_doublings==2,]
  CA_WGD<-CA[CA$Genome_doublings>0,]
  
  data_plus_minus <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_plus_minus_data.txt", header=T, sep='\t', as.is=T)
  

  CA_PM<-data_plus_minus[data_plus_minus$Type==ty,]
  CA_WGD0_plus_minus<-CA_PM[CA_PM$Genome_doublings==0,]
  CA_WGD1_plus_minus<-CA_PM[CA_PM$Genome_doublings==1,]
  CA_WGD2_plus_minus<-CA_PM[CA_PM$Genome_doublings==2,]
  CA_WGD_plus_minus<-CA_PM[CA_PM$Genome_doublings>0,]
  
  
  CA_WGD_mat <- CA_WGD[,14:(ncol(CA_WGD)-17)]
  CA_WGD0_mat <- CA_WGD0[,14:(ncol(CA_WGD0)-17)]
  
  CA_PM_WGD0<-c()
  CA_PM_WGD<-c()
  
  
  for (i in 14:91){
    WGD0_s <- sum(CA_WGD0_plus_minus[,i]==1, na.rm=T)
    WGD_s <- sum(CA_WGD_plus_minus[,i]==1, na.rm=T)
    CA_PM_WGD0[i]<- WGD0_s
    CA_PM_WGD[i]<- WGD_s
  }

  
  WGD_sum_desc <- sort(CA_PM_WGD, decreasing = T)
  WGD0_sum_desc <- sort(CA_PM_WGD0, decreasing = T)
  
  WGD_sum_top_3 <- WGD_sum_desc[1:3]
  WGD0_sum_top_3 <- WGD0_sum_desc[1:3]
  WGD_sum_total <- sum(WGD_sum_top_3, na.rm=T)
  WGD0_sum_total <- sum(WGD0_sum_top_3, na.rm=T)
  WGD_total <- sum(CA_PM_WGD, na.rm=T)
  WGD0_total <- sum(CA_PM_WGD0, na.rm=T)
  WGD_frac <- WGD_sum_total/WGD_total
  WGD0_frac <- WGD0_sum_total/WGD0_total
  WGD_fractions[t] <- WGD_frac
  WGD0_fractions[t] <- WGD0_frac

  
}
  
frac_top3_t <- c(WGD0_fractions, WGD_fractions)
wgd_status_t <- c(rep(as.character("non-WGD"), length(WGD0_fractions)),rep("WGD", length(WGD_fractions)))
type_t <- c(rep(as.character(g_20_types), 2))

new_data <- cbind(frac_top3_t, wgd_status_t, type_t)


df <- rbind(df, new_data)

t.test(WGD0_fractions,WGD_fractions, paired = T)