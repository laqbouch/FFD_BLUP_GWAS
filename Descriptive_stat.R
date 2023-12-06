#Compare phenotypic distribution per genetic cluster
library("AMR")
library("tidyr") 
library("ggplot2")
phenot <- read_excel("Pheno_7yrs_FFD_CP.xlsx")
phenot <- phenot[,c(1,5)]
clusers <- read_csv("SNMF_Matrix_K3.csv")

cluster_blup_data <- merge(clusers, phenot, by.x = "...1", 
                           by.y = "Name_on_cluster", all.x = T, all.y = F)
cluster_blup_data <- cluster_blup_data[!duplicated(cluster_blup_data), ]
q_df <- cluster_blup_data
q_df %>%
  add_column(cluster_SNP_70p = NA)
q_df$cluster_SNP_70p[q_df$P1>=0.7] <- "Eastern"
q_df$cluster_SNP_70p[q_df$P2>=0.7] <- "Central"
q_df$cluster_SNP_70p[q_df$P3>=0.7] <- "Western"
q_df$cluster_SNP_70p[q_df$P1<0.7 & q_df$P2<0.7 & q_df$P3<0.7] <- "M"

###violin plot per genetic cluster
violin_data<-q_df[!(q_df$cluster_SNP_70p=="M"),]
png(file = paste0(getwd(),"/Distribution_FFD_BLUP_per_cluster.png"),width = 570, height = 630)
 ggplot(violin_data, aes(x=violin_data$cluster_SNP_70p,y=violin_data$BLUP_FFD_7yrs,
                             color=violin_data$cluster_SNP_70p))+
  geom_violin(trim = FALSE)+scale_color_manual(values= c("#228B22","#004080","#FF6347"))+
geom_jitter(shape=16, position=position_jitter(0.2))+ geom_boxplot(width=0.1)+
  theme(legend.position="none") + stat_summary(fun.data="mean_sdl", mult=1, width=0.2 )+
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black")+
  labs(x="Genetic cluster", y = "FFD BLUP")+
  scale_x_discrete(limits=c("Eastern","Central","Western"),labels = c("Eastern","Central","Western"))
#ggsave(filename = paste0(getwd(),"/Distribution_FFD_BLUP_per_cluster.png"),
  #     plot = last_plot(), width = 570, height = 630, units = "px")
 dev.off()
 
#Wilcoxon test
wilcox_test_result <- pairwise.wilcox.test(violin_data$BLUP_FFD_7yrs, violin_data$cluster_SNP_70p, p.adjust.method = "none")
wilcox_test_result
#mean per cluster
mean_per_class <- violin_data %>%
   group_by(cluster_SNP_70p) %>%
   summarise(mean_BLUP_FFD_7yrs = round(mean(BLUP_FFD_7yrs), 2))
print(mean_per_class)
 
 