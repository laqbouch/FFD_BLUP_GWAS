#Title: "FFD_BLUP_GWAS_MM4LMM"
#Author: "Laila Aqbouch"
#contact: "aqbouchlaila@gmail.com"
#date: 05/12/2023

library(anyLib)
anyLib(c("devtools", "emma","apercu","corpcor","data.table","MM4LMM","qqman",
         "tibble","readxl","LEA","readr", "hierfstat",
         "MASS","PerformanceAnalytics","tidyverse","ggplot2","vioplot","ggpattern"))

###Import genotypic data
SNP <- read_delim("positions.txt", 
                  delim = "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE)
names <- read_csv("sample_names.txt", 
                  col_names = FALSE)

###Import imputed genotypic data by LFMM
imputed_geno_LFMM <- read_table("Oe9_CHR_scaff_WOGBM_NA_SNP_90_NA_indiv_25_MAC1_MAF_5_nuclear.lfmm_imputed.lfmm", 
                                col_names = FALSE)
genot.imp <- imputed_geno_LFMM
SNP$positions <- paste(SNP$X1,"_",SNP$X2, sep = "")
ap(genot.imp)

##add cultivar names
genot.mat <- as.matrix(genot.imp)
rownames(genot.mat) <- names$X1
colnames(genot.mat) <- SNP$positions

#MAF allele should be coded 012 to use this command line
p <- colMeans(genot.mat) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
ap(maf)

#visual support (facltatif)
#png(file="MAF_before_filter_imputed_SNMF.png")
#hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))
#dev.off()

##filter MAF at 0.05
sum(maf < 0.05)
genot.ok <- genot.mat[, maf >= 0.05]
dim(genot.ok)
##verif
p <- colMeans(genot.ok) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
#png(file="MAF_after_filter_0.05_imputed_SNMF.png")
#hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))
#dev.off()
dim(genot.ok)
ap(genot.ok)
#clean environement
rm(genot.mat, genot.imp, maf, p, q)

####Carte
map <- SNP[, c(3, 1, 2)]
colnames(map) <- c('SNP','Chr','Pos')
# Make the file simple to use (Replace s of scaff by 9)
map$Chr <- sub('9','',map$Chr)
map$Chr <- sub('s','9',map$Chr)
map$Chr <- gsub('[Oe_LG]','',map$Chr)
ap(map)
tail(map)
#write.csv(map , "map.csv", row.names=T)

# filter conserved SNP
map <- map[map$SNP %in% colnames(genot.ok), ]
map$Chr <- as.numeric(map$Chr)
#order by chromosome and position 
map <- map[order(map$Pos), ]
map <- map[order(map$Chr), ]
head(map)
tail(map)
#write.csv(map, "ordered_map.csv", row.names=T)

####Import phenotypic data
phenot <- read_excel("Pheno_7yrs_FFD_CP.xlsx")
head(phenot)
phenot <- phenot[,c(1,5)]
##FFD data
y_FFD <- phenot$BLUP_FFD_7yrs
names(y_FFD) <- phenot$Name_on_cluster

#Phenotypic data distribution
y <- y_FFD
length(y)
ap(y)
#summary stat
summary(y)
###test of normality
shapiro.test(y)
# plot phenotype distribution
png(file="hist_phenotype_distribution_FFD_7yrs.png",width =575 ,height = 733,res = 120)
vioplot(y, horizontal = FALSE, col = "pink", names = "FFD_BLUP", main = "Distribution of BLUP FFD")
dev.off()

#####VerifY data sorted 
#Order geno to fit pheno file
genot.ok <- genot.ok[match(names(y), rownames(genot.ok)), ]
#png(file="verif_genot_names.png")
plot(match(rownames(genot.ok), names(y)))
#dev.off()
#png(file="verif_genot_SNP.png")
plot(match(map$SNP, colnames(genot.ok)))
#dev.off()

#####create Kinship martix using Weir and Goudet (2017) beta estimator 
ap(genot.ok)
dim(genot.ok)
K <- beta.dosage(genot.ok)
K <- make.positive.definite(K)
is.positive.definite(K)
ap(K)
#write.csv(K, file = "K_values.csv", row.names = FALSE)

###structure matrix
Q.mat <- sapply(SNMF_Matrix_K3[,2:ncol(SNMF_Matrix_K3)], as.numeric)
rownames(Q.mat) <- SNMF_Matrix_K3$...1
colnames(Q.mat)
dim(Q.mat)
ap(Q.mat)
ap(genot.ok)

#Compare three models: model with kinship only, model with structure only, and model with both
#calcul of the AIC and BIC based on model without geno data
mmest_K_test <- MMEst(Y = y,VarList = list(Additive = K, Error = diag(length(y))))
ap(mmest_K_test)
AIC_K = -2*(mmest_K_test$NullModel$`LogLik (Reml)`) +2*1 
BIC_K = -2 * (mmest_K_test$NullModel$`LogLik (Reml)`)+ log(318)*1 

mmest_QK_test <- MMEst(Y = y, Cofactor = Q.mat[, 1:2],
                       VarList = list(Additive = K, Error = diag(length(y))))
ap(mmest_QK_test)
AIC_QK = -2*(mmest_QK_test$NullModel$`LogLik (Reml)`) +2*3 
BIC_QK = -2 * (mmest_QK_test$NullModel$`LogLik (Reml)`)+ log(318)*3 

mmest_Q_test <- MMEst(Y = y, Cofactor = Q.mat[, 1:2],
                      VarList = list(Error = diag(length(y))))
ap(mmest_Q_test)
AIC_Q = -2*(mmest_Q_test$NullModel$`LogLik (Reml)`) +2*2
BIC_Q = -2 * (mmest_Q_test$NullModel$`LogLik (Reml)`)+ log(318)*2 

###GWAS of FFD BLUP using MM4LMM, model with Kinship only is the best based on AIC and BIC
mmest_K <- MMEst(Y = y, X = genot.ok, VarList = list(Additive = K, Error = diag(length(y))), Verbose = TRUE,Cofactor = NULL)
#ap(mmest_K)

# Calcul de la covariance des estimateurs de variance
#Estimation of effects
covariance_K <- MMVcov(ResMM = mmest_K, Y = y, VarList = list(Additive = K, Error = diag(length(y))))

#test of the significance of effects
out.test_K <- AnovaTest(mmest_K)
ap(out.test_K)
res.mm4lmm_K <- cbind(map, P = sapply(out.test_K, function(x){x["Xeffect", "pval"]}))
ap(res.mm4lmm_K)

#Explort results (facultatif)
write.csv(res.mm4lmm_K, file = paste0(getwd(),"/Result_mm4lmm_p_values.csv"))

#QQ_plot with -log(p_values)
png(file=paste0(getwd(),"/QQ_plot_mm4lmm_FFD_K.png"))
qq(res.mm4lmm_K$P, main = "QQ_plot_GWAS_MM4LMM_FFD_K")
dev.off()
#QQ_plot with real values of p_values
u <- runif(118948,0,1)
png(file=paste0(getwd(),"/QQ__true_p_values_mm4lmm_FFD_K.png"))
qqplot(u,res.mm4lmm_K$P, main = "QQ_p_values_mm4lmm_FFD_7yrs_K")
dev.off()

#Manhattan plot with only SNPS on chromosome
res.mm4lmm_K_chr <- subset(res.mm4lmm_K,Chr%in%1:23)
tail(res.mm4lmm_K_chr)
png(file=paste0(getwd(),"/GWAS_plot_mm4lmm_FFD_K_FDR5_CHR.png"),res = 100,width = 1000, height = 400)
manhattan(res.mm4lmm_K_chr, chr = "Chr", bp = "Pos", p = "P",
          suggestiveline = -log10(9.60E-06),
          col = c("gray0", "chartreuse4"),  
          genomewideline = -log10(9.60E-06),
          main = "GWAS_MM4LMM_FFD_WG_SNMFimp_K_FDR_5",
          ylim = c(0, 7.5),
          annotatePval = 0.00001)
dev.off()

###FDR correction
p_values <- res.mm4lmm_K$P
p_adjusted <- p.adjust(p_values, method = "fdr")
p_adjusted_coor <- cbind(res.mm4lmm_K$SNP,res.mm4lmm_K$P,p_adjusted)
colnames(p_adjusted_coor)= c("SNP","P","p_adjusted")
p_adjusted_coor <- as.data.frame(p_adjusted_coor)
#ap(p_adjusted_coor)

###Extract significant (Threshold at 5% FDR)
signifSNPs <- p_adjusted_coor[p_adjusted_coor$p_adjusted < 0.05,]
signifSNPs
write.csv(signifSNPs, paste0(getwd(),"/K/result_mm4lmm_K_sig_SNP_FDR_05.csv"))

###Significant SNP effets 
effects_sig_SNP <- sapply(mmest_K[signifSNPs$SNP], function(x) x$Beta)[2, ]
# Convert to data frame
effects_sig_SNP_df <- data.frame(SNP_name = names(effects_sig_SNP), effects_sig_SNP = effects_sig_SNP, row.names = NULL)
# Write the data frame to a file in your working directory
write.table(effects_sig_SNP_df, file = "effects_sig_SNP_FDR5.csv", sep = ",", col.names = TRUE, row.names = TRUE)

#Allelic repartition per genetic cluster
SNMF_Matrix_K3 <- read_csv("SNMF_Matrix_K3.csv")
SNMF_Matrix_K3 <- SNMF_Matrix_K3[match(rownames(genot.ok), SNMF_Matrix_K3$...1),]
#ap(SNMF_Matrix_K3)
#verify order
plot(match(rownames(genot.ok), SNMF_Matrix_K3$...1))

#classify individuals by genetic cluster
q_df <- SNMF_Matrix_K3
q_df <- q_df %>%
  add_column(cluster_SNP_70p = NA)
q_df$cluster_SNP_70p[q_df$P1>=0.7] <- "Eastern"
q_df$cluster_SNP_70p[q_df$P2>=0.7] <- "Central"
q_df$cluster_SNP_70p[q_df$P3>=0.7] <- "Western"
q_df$cluster_SNP_70p[q_df$P1<0.7 & q_df$P2<0.7 & q_df$P3<0.7] <- "M"
cluster_info <- q_df[,c(1,5)]
# Iterate over each SNP
for (i in signifSNPs$SNP) {
    # Combining the specific SNP data with cluster_info
  combined_data <- data.frame(cluster_info, genot.ok[, i])
  combined_data <- combined_data[combined_data$cluster_SNP_70p != "M", ]
  names(combined_data)[3] <- "geno"
  
  # Count the number of genotypes per cluster
  genotype_counts <- table(combined_data$cluster_SNP_70p)
  #beta value for the SNP
  beta_value <- effects_sig_SNP_df[effects_sig_SNP_df$SNP_name == i, "effects_sig_SNP"]
  
  # Create the bar plot for the specific SNP with annotations
  #png(file = paste0(getwd(),"/Allelic_repartition_per_cluster_for",i,".png"))
  ggplot(combined_data, aes(x = factor(cluster_SNP_70p, levels = c("Eastern", "Central", "Western")))) +
    geom_bar(aes(fill = factor(geno)), position = "stack") +
    scale_fill_manual(values = c("#FF5733", "#5DADE2", "#2ECC71")) +  
    labs(x = "Genetic cluster", y = "Count", title = paste("Allelic repartition for", i)) +
    geom_text(stat='count', aes(label = ..count.., group = geno), position = position_stack(vjust = 0.5), show.legend = FALSE) +
    annotate("text", x = Inf, y = Inf, label = paste("Beta = ", round(beta_value, 2)), hjust = 1, vjust = 1, size = 4, fontface = "bold") +
    guides(fill = guide_legend(title = "Genotype"))
  
  ggsave(filename = paste0(getwd(), "/Allelic_repartition_per_cluster_for", i, ".png"),
         plot = last_plot(), width = 1200, height = 1000, units = "px")
    # Save the current plot
  #dev.off()
}