# Get the working directory
getwd()
# Set the working directory
setwd("./acc_pca")
getwd()
# Get the presence-absence matrix
p_a <- read.table("gene_cluster_presence_absence.txt", header = T, sep = "\t")
s <- matrix(data = p_a$value, nrow = length(unique(p_a$layer)))
colnames(s) <- unique(p_a$item)
rownames(s) <- unique(p_a$layer)
s <- as.data.frame(s)
# Get the similarity matrix with Sokal and Michel coeficents
SM.simil <- function(x,y){
   if (length(x) != length(y)){
     stop("los vectores han de tener la misma longitud")} else {
       a <- sum(x==1 & y==1)
       d <- sum(x==0 & y==0)
       (a+d)/length(x)
       }
}
bin.simil <- function(X,simil){
   n <- nrow(X)
   S <- matrix(0,nrow=n,ncol=n)
   for(i in 1:(n-1)) for(j in (i+1):n) S[i,j] <- simil(X[i,],X[j,])
   S <- S + diag(1,n) + t(S)
   rownames(S) <- colnames(S) <- rownames(X)
   as.dist(2*(1-S))
}
d <- bin.simil(s, SM.simil)
# Multidimensional scale of distances
mds <- eigen(d)
# Explained variance
comp_var <- ((mds$values)**2/sum((mds$values)**2))
plot(1:287, comp_var)
which(comp_var > 0.3)
sum(comp_var[which(comp_var > 0.3)])
# PC plot
# Country names in factor string
library(stringr)
country <- str_extract(rownames(s), "[A-Z][a-z]+")
country[which(is.na(country))] <- "USA"
country <- factor(country)
# Year in factor string
year <- str_extract(rownames(s), "[0-9]+$")
year <- as.numeric(year)
mds$country <- country
mds$year <- year
# Continent
continent <- as.character(country)
for (i in 1:287) {
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
mds$continent <- factor(continent)
# Plot with variables
pca <- data.frame(mds$vectors[,which(comp_var > 0.3)], country, year, continent)
library(ggplot2)
ggplot(pca, aes(x=X1, y=X2, color=year, label=country)) + geom_text() + 
   scale_color_gradient(low = "blue", high = "red")
ggsave("raw_pca", device = "png")

ggplot(pca, aes(x=X1, y=X2, color=year, label=country)) + geom_text() + 
   scale_color_gradient(low = "blue", high = "red") +
   xlim(-0.123,-0.095) +
   ylim(0.09,0.14)
ggsave("close_left_up", device = "png")
acc <- pca[pca$X1 < -0.095 & pca$X2 > 0.09, c(3,4)]
coin <- matrix(0, ncol = ncol(acc))
colnames(coin) = colnames(acc)
for (i in 1:nrow(acc)) {
   for (j in 1:nrow(core)) {
      if (acc$country[i] == core$country[j] & acc$year[i] == core$year[j]) {
         coin <- rbind(coin, acc[i,])
         num_coin <- append(num_coin, i)}
      # else {print(acc[i,])}
   }
}
unique(coin)
ggplot(pca, aes(x=X1, y=X2, color=year, label=country)) + geom_text() + 
   scale_color_gradient(low = "blue", high = "red") +
   xlim(-0.07, -0.04) +
   ylim(0,-0.05)
ggsave("close_right_down", device = "png")
ggplot(pca[pca$continent != "Americas",], aes(x=X1, y=X2, color=year, label=country)) + geom_text() + 
   scale_color_gradient(low = "blue", high = "red") +
   xlim(-0.07, -0.04) +
   ylim(0,-0.05)
acc2 <- pca[pca$X1 > -0.095 & pca$X2 < 0.09, c(3,4)]
coin2 <- matrix(0, ncol = ncol(acc))
colnames(coin2) = colnames(acc2)
for (i in 1:nrow(acc2)) {
   for (j in 1:nrow(core2)) {
      if (acc2$country[i] == core2$country[j] & acc2$year[i] == core2$year[j]) {
         coin2 <- rbind(coin2, acc2[i,])
         num_coin2 <- append(num_coin2, i)}
      # else {print(acc[i,])}
   }
}
unique(coin2)
   ggplot(pca, aes(x=X1, y=X2, color=year, label=country)) + geom_text() + 
   scale_color_gradient(low = "blue", high = "red") +
   xlim(-0.07, -0.04) +
   ylim(0,-0.05)
