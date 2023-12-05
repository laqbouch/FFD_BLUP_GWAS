#Title: "FFD_BLUP_GWAS_MLMM"
#Author: "Laila Aqbouch"
#contact: "aqbouchlaila@gmail.com"
#date: 05/12/2023

library(anyLib)
anyLib(c("devtools","apercu","data.table","mlmm","readr","readxl","hierfstat","corpcor","tibble","emma"))

###Import geno data
SNP <- read_delim("positions.txt", 
                  delim = "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE)

ap(SNP)
names <- read_csv("sample_names.txt", 
                  col_names = FALSE)
ap(names)
###Genotype
imputed_geno_LFMM <- read_table("Oe9_CHR_scaff_WOGBM_NA_SNP_90_NA_indiv_25_MAC1_MAF_5_nuclear.lfmm_imputed.lfmm", 
                                col_names = FALSE)
genot.imp <- imputed_geno_LFMM
SNP$positions <- paste(SNP$X1,"_",SNP$X2, sep = "")
ap(genot.imp)

##add cultivar names
dim(genot.imp)
ap(genot.imp) 
genot.mat <- as.matrix(genot.imp)
ap(genot.mat)
rownames(genot.mat) <- names$X1
colnames(genot.mat) <- SNP$positions
class(genot.mat)
dim(genot.mat)
genot.mat[1:20,1:3]
#MAF allele should be coded 012
p <- colMeans(genot.mat) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
#plot MAF
png(file="MAF_before_filter_imputed_SNMF.png")
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))
dev.off()
##filter MAF at 0.05
sum(maf < 0.05)
genot.ok <- genot.mat[, maf >= 0.05]
dim(genot.ok)
##verif
p <- colMeans(genot.ok) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
png(file="MAF_after_filter_0.05_imputed_SNMF.png")
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))
dev.off()
dim(genot.ok)
ap(genot.ok)
#clean environement
rm(genot.mat, genot.imp, maf, p, q)

####Carte
map <- SNP[, c(3, 1, 2)]
colnames(map) <- c('SNP','Chr','Pos')
class(map)
dim(map)
ap(map)
# Replace String with Another Stirng, we replmace s of scaff by 9
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

####Import pheno data
phenot <- read_excel("Pheno_7yrs_FFD_CP.xlsx")
head(phenot)
phenot <- phenot[,c(1,5)]
##FFD data
y_FFD <- phenot$BLUP_FFD_7yrs
names(y_FFD) <- phenot$Name_on_cluster
#verifier la distribution
y <- y_FFD
length(y)
ap(y)
#####verif data sorted
#order geno to fit pheno
genot.ok <- genot.ok[match(names(y), rownames(genot.ok)), ]
plot(match(rownames(genot.ok), names(y)))
plot(match(map$SNP, colnames(genot.ok)))
#####Matrice of genetic relationship 
K <- beta.dosage(genot.ok)
K <- make.positive.definite(K)
is.positive.definite(K)

###MLMM with Kinship in the model and the same treshold as MM4LMM GWAS
mygwas <- mlmm(Y = y, X = genot.ok, K = K, maxsteps = 7, nbchunks = 2,thresh = 0.0000096)
#summary of all steps
step_table <- mygwas$step_table
#variance component
variance_repartition <-mygwas$RSSout
#Partition de variance
png(file="variance_repartition_7steps_k.png")
plot_step_RSS(mygwas)
dev.off()

####Manhattan plots
#Step 1 : pas de cofacteur
png(file="GWAS_MLMM_step1_7_K.png")
plot_fwd_GWAS(x = mygwas, step = 1, snp_info = map, pval_filt = 0.1)
dev.off() 
##Step 2 : 1 cofacteur
png(file="GWAS_MLMM_step2_7_k.png")
plot_fwd_GWAS(x = mygwas, step = 2, snp_info = map, pval_filt = 0.1)
dev.off() 
##Step 3 : 2 cofacteur
png(file="GWAS_MLMM_step3_7_k.png")
plot_fwd_GWAS(x = mygwas, step = 3, snp_info = map, pval_filt = 0.1)
dev.off()

#selection of best model regarding the threshold
png(file="MLMM_best_model_thresh_k.png")
plot_step_table(mygwas, "maxpval")
dev.off()
#qq plot of the best model regarding the threshold
png(file="qqplot_best_model_thresh.png")
qqplot_opt_GWAS(x = mygwas, opt = "thresh")
dev.off()
#Manhatten plot of best model regarding the threshold
png(file="GWAS_MLMM_optimal_thresh.png")
plot_opt_GWAS(x = mygwas,opt = "thresh" , snp_info = map, pval_filt = 0.1)
dev.off()
#QQplots
png(file="qqplot_7steps.png",width = 1000,height = 700,res = 100)
qqplot_fwd_GWAS(x = mygwas, nsteps = 7)
dev.off()

#Illustration of effects
for (i in mygwas$opt_thresh$cof) {
  png(file=paste0("Effect_signif_MLMM_SNP_5_FDR_K",i,".png"))
  boxplot(y ~ genot.ok[, i], varwidth = T, main = i, xlab = "")
  mtext(text = c("n = ", table(genot.ok[, i])), side = 1,
        at = c(0.5, 1, 2), line = 3)
  dev.off()
}

#Coefficients et p-valeurs des cofacteurs
sig_SNP_FDR5 <- mygwas$opt_thresh$coef
write.csv(sig_SNP_FDR5 , "sig_SNP_FDR5_MLMM.csv", row.names=T)
