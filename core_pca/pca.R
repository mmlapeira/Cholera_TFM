# Get the working directory
getwd()
# Set the working directory
setwd("./")
getwd()
# Read the eigen vectors table
pca1 <- read.table("./core_pca/288_genomes.eigenvec",sep = "\t", header = T)
# Plot two principal components
plot(data=pca1, PC2~PC1)
# Read the eigen values table
pca_val <- read.table("./core_pca/288_genomes.eigenval", header = F)
# Variance explained by each PC
var_n <- pca_val**2/sum(pca_val**2)
var_n
# Country names in factor string
library(stringr)
country <- str_extract(pca1$IID, "[A-Z][a-z]+")
country[which(is.na(country))] <- "USA"
country
country <- factor(country)
# Year in factor string
year <- str_extract(pca1$IID, "[0-9]+$")
year <- as.numeric(year)
pca1$country <- country
pca1$year <- year
# Continent
continent <- as.character(country)
for (i in 1:277) {
  if ((country[i]) == "Angola") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Bangladesh") {
    continent[i] = "Asia"
  }
  else if ((country[i]) == "Burkina") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "China") {
    continent[i] = "Asia"
  }
  else if ((country[i]) == "Cote") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "India") {
    continent[i] = "Asia"
  }
  else if ((country[i]) == "Japan") {
    continent[i] = "Asia"
  }
  else if ((country[i]) == "Mali") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Mauritania") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Nigeria") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Russia") {
    continent[i] = "Europe"
  }
  else if ((country[i]) == "Sao") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Sri") {
    continent[i] = "Asia"
  }
  else if ((country[i]) == "Uganda") {
    continent[i] = "Africa"
  }
  else if ((country[i]) == "Vietnam") {
    continent[i] = "Asia"
  }
  else {
    continent[i] = "Americas"
  }
}
continent <- factor(continent)
pca1$continent <- factor(continent)
# Plot with variables
library(ggplot2)
ggplot(pca1, aes(x=PC1, y=PC2, color=year, label=country)) + geom_text() + 
  scale_color_gradient(low = "blue", high = "red")
ggsave("raw_pca", device = "png")

ggplot(pca1, aes(x=PC1, y=PC2, color=year, label=country)) + geom_text() + 
  scale_color_gradient(low = "blue", high = "red") +
  xlim(0.005, 0.015) +
  ylim(-0.025,0.01)
ggsave("close1_pca", device = "png")

ggplot(pca1, aes(x=PC1, y=PC2, color=year, label=country)) + geom_text() + 
  scale_color_gradient(low = "blue", high = "red") +
  xlim(0.013, 0.0145) +
  ylim(-0.01,0.005)
ggsave("close2_pca", device = "png")

#Separate populations
pop_1 <- pca1$IID[pca1$PC1>0.005]
length(pop_1)
pop_2 <- pca1$IID[pca1$PC1<0.005]
length(pop_2)
core <- pca1[pca1$PC1<0.005, c(13,14)]
core2 <- pca1[pca1$PC1>0.005, c(13,14)]
# Animated plot
library(gganimate)
library(gifski)
myPlot <- ggplot(pca1, aes(x=PC1, y=PC2, label=country)) + geom_text() +
  labs(title = "Year: {frame_time}") +
  transition_time(year) +
  ease_aes("linear")
animate(myPlot, duration = 10, fps = 20, width = 200, height = 200, renderer = gifski_renderer())
anim_save("./ggplot_animated.gif")

# Plot with continents
ggplot(pca1, aes(x=PC1, y=PC2, color=year, shape=continent)) + geom_point() + 
  scale_color_gradient(low = "blue", high = "red")

ggplot(pca1[pca1$continent == "Asia",], aes(x=PC1, y=PC2, color=year, shape=continent)) + geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +
  xlim(0.005, 0.015) +
  ylim(-0.025,0.01)

ggplot(pca1[pca1$continent == "Americas",], aes(x=PC1, y=PC2, color=year, label=country)) + 
  geom_text() + 
  scale_color_gradient(low = "blue", high = "red") +
  xlim(0.013, 0.0145) +
  ylim(-0.01,0.005)
ggsave("close2_min_america_pca", device = "png")

  # Population metrics
pi.all <- read.table("core_pca/merged_288_genomes_10kb.windowed.pi", header = T)
png("hist_pi.png")
hist(pi.all$PI, br=20)
dev.off()
png("boxplot_pi.png")
boxplot(pi.all$PI,ylab="diversity")
dev.off()
png("pi_position.png")
plot(pi.all$BIN_START,pi.all$PI,xlab="position",ylab="diversity")
abline(h=5e-05, col="red", lty=2)
dev.off()
taj.all <- read.table("core_pca/merged_288_genomes_10kb.Tajima.D", header = T)
png("taj_hist.png")
hist(taj.all$TajimaD,br=20)
dev.off()
png("taj_pos.png")
plot(taj.all$BIN_START,taj.all$TajimaD,xlab="position",ylab="Tajima's D")
dev.off()
fst <- read.table("core_pca/merged_288_genomes_10kb_pop1_vs_pop2.windowed.weir.fst", header = T)
png("fst_hist.png")
hist(fst$WEIGHTED_FST, br=20)
dev.off()
png("fst_boxplot.png")
boxplot(fst$WEIGHTED_FST)
dev.off()
png("fst_pos.png")
plot(fst$BIN_START,fst$WEIGHTED_FST,xlab = "position", ylab = "Fst")
abline(a=0.01, b=0, col="red", lty=2)
abline(a=-0.01, b=0, col="red", lty=2)
abline(a=0.1, b=0, col="blue", lty=2)
abline(a=-0.1, b=0, col="blue", lty=2)
dev.off()
# Population comparison pi
pi_91 <- read.table("core_pca/1991_samples_10kb.windowed.pi", header = T)
pi_95 <- read.table("core_pca/1995_samples_10kb.windowed.pi", header = T)
pi_00 <- read.table("core_pca/2000_samples_10kb.windowed.pi", header = T)
png("pi_comparison.png")
par(mfrow=c(3,1))
plot(pi_91$BIN_START,pi_91$PI,xlab="position",ylab="diversity", main = "Pi in 1991 set")
abline(h=5e-05, col="red", lty=2)
plot(pi_95$BIN_START,pi_95$PI,xlab="position",ylab="diversity", main = "Pi in 1995 set")
abline(h=5e-05, col="red", lty=2)
plot(pi_00$BIN_START,pi_00$PI,xlab="position",ylab="diversity", main = "Pi in 2000 set")
abline(h=5e-05, col="red", lty=2)
dev.off()
# Population comparison Taj
taj_91 <- read.table("core_pca/1991_samples_10kb.Tajima.D", header = T)
taj_95 <- read.table("core_pca/1995_samples_10kb.Tajima.D", header = T)
taj_00 <- read.table("core_pca/2000_samples_10kb.Tajima.D", header = T)
png("taj_comparison.png")
par(mfrow=c(3,1))
plot(taj_91$BIN_START,taj_91$TajimaD,xlab="position",ylab="Tajima's D", main = "1991 set")
plot(taj_95$BIN_START,taj_95$TajimaD,xlab="position",ylab="Tajima's D", main = "1995 set")
plot(taj_00$BIN_START,taj_00$TajimaD,xlab="position",ylab="Tajima's D", main = "2000 set")
dev.off()
# Fst comparison
fst_1 <- read.table("core_pca/1991_vs_1995_10kb.windowed.weir.fst", header = T)
fst_2 <- read.table("core_pca/1995_vs_2000_10kb.windowed.weir.fst", header = T)
png("fst_comparison.png")
par(mfrow=c(2,1))
plot(fst_1$BIN_START,fst_1$WEIGHTED_FST,xlab = "position", ylab = "Fst", main = "1991 vs 1995")
abline(a=0.01, b=0, col="red", lty=2)
abline(a=-0.01, b=0, col="red", lty=2)
abline(a=0.1, b=0, col="blue", lty=2)
abline(a=-0.1, b=0, col="blue", lty=2)
plot(fst_2$BIN_START,fst_2$WEIGHTED_FST,xlab = "position", ylab = "Fst", main = "1995 vs 2000")
abline(a=0.01, b=0, col="red", lty=2)
abline(a=-0.01, b=0, col="red", lty=2)
abline(a=0.1, b=0, col="blue", lty=2)
abline(a=-0.1, b=0, col="blue", lty=2)
dev.off()
