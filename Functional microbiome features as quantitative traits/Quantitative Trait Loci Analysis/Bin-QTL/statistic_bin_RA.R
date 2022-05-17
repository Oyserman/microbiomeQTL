out_folder <- "/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/qtl/metagenome_qtl/bin_level"
big_phenofile <- file.path(out_folder, "phenofile_all_good_bins.csv")

big_RA <- read.csv(big_phenofile)
par(mfrow=c(1,1))
boxplot(big_RA[,c(-1)])

bin_name <- strsplit(colnames(big_RA)[-1], split=".tsv", fixed=TRUE)
b <- unlist(lapply(bin_name,head,1))
colnames(big_RA)[-1] <- b

# test data nomality
normal_distributed <- c()
for (i in colnames(big_RA)[-1]){
  if (shapiro.test(big_RA[,i])$p.value > 0.05){
    normal_distributed <- c(normal_distributed, i)
    cat(i,":", shapiro.test(big_RA[,i])$p.value, "\n")
  }
}

for (i in colnames(big_RA)[-1]){
  if (shapiro.test(log(big_RA[,i]))$p.value < 0.05){
    cat(i,":", shapiro.test(log(big_RA[,i]))$p.value, "\n")
  }
}
# ________________________________________

if(!require(dplyr)){install.packages("dplyr")}
if(!require(FSA)){install.packages("FSA")}
if(!require(car)){install.packages("car")}
if(!require(agricolae)){install.packages("agricolae")}
if(!require(multcomp)){install.packages("multcomp")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(Rmisc)){install.packages("Rmisc")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(pwr)){install.packages("pwr")}

reformated_data <- pivot_longer(big_RA, cols = starts_with("bin."),
             names_to = "bin",
             names_prefix = "",
             values_to = "RA",
             values_drop_na = TRUE)
leveneTest(reformated_data$RA, reformated_data$bin)

library(dplyr)

reformated_data <-
  mutate(reformated_data,
         bin = factor(bin, levels=unique(bin)))

library(FSA)   

reformated_data$log2_RA <- log2(reformated_data$RA)
Summarize(log2_RA~bin,reformated_data)
model <- lm(log2_RA ~ bin, 
           data=reformated_data)
anova(model)
hist(residuals(model),col="darkgray")
plot(fitted(model), 
     residuals(model))

?p.adjust

library(agricolae)
lsd <- LSD.test(model, "bin",   # outer parentheses print result
         alpha = 0.05,       
         p.adj="BH")

TukeyHSD(aov(model), p.adj = "BH" )

x <- seq(1,33)
y <- aggregate(reformated_data$log2_RA, by=list(reformated_data$bin), FUN =max)
y <- y$x
newdf <- data.frame(x,y)
g <- ggplot(reformated_data,aes(x= reorder(bin, log2_RA, FUN = mean), y= log2_RA))
g + geom_boxplot() +
  ylab("log2 relative abundance") +
  xlab("gemomic bin") +
  geom_hline(yintercept = mean(reformated_data$log2_RA), linetype = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  stat_compare_means(method = "anova", label.y = -11, label.x = 3) +
  geom_text(data=newdf, aes(x=x,y=y+0.5,label=lsd$groups$groups))


str(lsd$groups$groups)
