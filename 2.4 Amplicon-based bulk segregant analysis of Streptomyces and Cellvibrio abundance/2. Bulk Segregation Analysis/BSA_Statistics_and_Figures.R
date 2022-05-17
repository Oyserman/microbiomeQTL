setwd(dir = "~/Desktop/QTL/For_submission/Analysis_Folder/BSA/")
library(readxl)
library("ggpubr")
library("dplyr")
library(gridExtra)
library("ggstar")
library(colorspace)
library(agricolae)

BSA <- read_excel("BSA_RIL1_summary.xlsx", sheet = "Summary 2")


BSA$SNP_2274 <- factor(BSA$SNP_2274, levels = c("M", "P", "A", "B", "Bulk"))
BSA$SNP_464 <- factor(BSA$SNP_464, levels = c("M", "P", "A", "B", "Bulk"))
BSA$SNP_3142 <- factor(BSA$SNP_3142, levels = c("M", "P", "A", "B", "Bulk"))

SNP_2274.aov <- aov(BSA$ASV3	 ~ BSA$SNP_2274)
summary(SNP_2274.aov) # 0,05 is significant
TukeyHSD(SNP_2274.aov)
res_SNP_2274 <- TukeyHSD(SNP_2274.aov, ordered = TRUE)
stats_SNP_2274 <- as.data.frame(res_SNP_2274$`BSA$SNP_2274`)
group_SNP_2274 <- HSD.test(SNP_2274.aov, "BSA$SNP_2274",group = TRUE)


SNP_464.aov <- aov(BSA$ASV9	 ~ BSA$SNP_464)
summary(SNP_464.aov) # 0,05 is significant
TukeyHSD(SNP_464.aov)
res_SNP_464 <- TukeyHSD(SNP_464.aov, ordered = TRUE)
stats_SNP_464 <- as.data.frame(res_SNP_464$`BSA$SNP_464`)
group_SNP_464 <- HSD.test(SNP_464.aov, "BSA$SNP_464",group = TRUE)

SNP_3142.aov <- aov(BSA$ASV9	 ~ BSA$SNP_3142)
summary(SNP_3142.aov) # 0,05 is significant
TukeyHSD(SNP_3142.aov)
res_SNP_3142 <- TukeyHSD(SNP_3142.aov, ordered = TRUE)
stats_SNP_3142 <- as.data.frame(res_SNP_3142$`BSA$SNP_3142`)
group_SNP_3142 <- HSD.test(SNP_3142.aov, "BSA$SNP_3142",group = TRUE, alpha = 0.1)




SNP_2274_metadata <- cbind(table(BSA$SNP_2274),
group_SNP_2274$groups[c(1,3,2,4,5),])

group_SNP_2274$means[c(4,5,1,2,3),7]

SNP_464_metadata <- cbind(table(BSA$SNP_464),
group_SNP_464$groups[c(3,2,4,1,5),])
group_SNP_464$means[c(4,5,1,2,3),7]

SNP_3142_metadata <- cbind(table(BSA$SNP_3142),
group_SNP_3142$groups[c(3,2,4,1,5),])
group_SNP_3142$groups

group_SNP_3142$means[c(4,5,1,2,3),7]


SNP_2274 <- ggplot(data = BSA, aes(SNP_2274, ASV3)) +
  geom_boxplot(data = BSA, aes(SNP_2274, ASV3, col = SNP_2274), outlier.colour="black", outlier.shape=NA, outlier.size=1, notch=FALSE) +
  geom_jitter(width = .2, aes(color = SNP_2274, alpha = 0.25), shape = 1, size = 0.5) +
#  annotate("text", x = c(0.5,0.5,0.5), y=c(1000 * c(1.200,1.150,1.100)), label = c("G =", "N = ", "A = "), cex = 2) +
#  annotate("text", x = c(1,2,3,4,5), y=1000*1.200, label = as.character(floor(SNP_2274_metadata$'BSA$ASV3')), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=1000*1.20, label = as.character(SNP_2274_metadata$Freq), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=1000*1.100, label = as.character(SNP_2274_metadata$groups), cex = 2) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

SNP_464 <- ggplot(data = BSA, aes(SNP_464, ASV9)) +
  geom_boxplot(data = BSA, aes(SNP_464, ASV9, col = SNP_464), outlier.colour="black", outlier.shape=NA, outlier.size=1, notch=FALSE) +
  geom_jitter(width = .2, aes(color = SNP_464, alpha = 0.25), shape = 1, size = 0.5) +
 # annotate("text", x = c(1,2,3,4,5), y=200*1.200, label = as.character(floor(SNP_464_metadata$'BSA$ASV9')), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=200*1.20, label = as.character(SNP_464_metadata$Freq), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=200*1.100, label = as.character(SNP_464_metadata$groups), cex = 2) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

SNP_3142 <- ggplot(data = BSA, aes(SNP_3142, ASV9)) +
  geom_boxplot(data = BSA, aes(SNP_3142, ASV9, col = SNP_3142), outlier.colour="black", outlier.shape=NA, outlier.size=1, notch=FALSE) +
  geom_jitter(width = .2, aes(color = SNP_3142, alpha = 0.25), shape = 1, size = 0.5) +
#  annotate("text", x = c(1,2,3,4,5), y=200*1.200, label = as.character(floor(SNP_3142_metadata$'BSA$ASV9')), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=200*1.20, label = as.character(SNP_3142_metadata$Freq), cex = 2) +
  annotate("text", x = c(1,2,3,4,5), y=200*1.100, label = as.character(SNP_3142_metadata$groups), cex = 2) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") 


### Figure 5
p123 <- ggarrange(SNP_464, SNP_3142, SNP_2274, ncol = 3, nrow = 1) 

