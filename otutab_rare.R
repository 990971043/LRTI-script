
print(paste("The input feature table is ", opts$input,  sep = ""))
# print(paste("Normalized filename: ", opts$normalize,  sep = ""))
# print(paste("Output alpha diversity: ", opts$output, sep = ""))

suppressWarnings(dir.create("alpha/"))

package_list <- c("vegan")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


species = read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="") 
# colSums(taxonomy)

print(paste0("Samples size are:"))
colSums(species)
min = min(colSums(species))
if (opts$depth==0){
  opts$depth=min}
print(paste("Rarefaction depth ", opts$depth, ". If depth set 0 will using sample minimum size ", min, sep = ""))

print(paste("Random sample number: ", opts$seed,  sep = ""))
set.seed(opts$seed)
otu = vegan::rrarefy(t(species), opts$depth)
# print(paste0("All sample rarefaction as following"))
# rowSums(otu)

estimateR = t(estimateR(otu))[,c(1,2,4)]
colnames(estimateR) = c("richness", "chao1", "ACE")

shannon = diversity(otu, index = "shannon")
simpson = diversity(otu, index = "simpson")
invsimpson = diversity(otu, index = "invsimpson")

alpha_div = cbind(estimateR, shannon, simpson, invsimpson)
print(paste0("Calculate six alpha diversities by estimateR and diversity"))
head(alpha_div, n=1)

write.table("#OTUID\t", file=paste(opts$normalize,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)

suppressWarnings(write.table(t(otu), file=paste(opts$normalize,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

write.table("Alpha_diversity\t", file=paste(opts$output,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)

suppressWarnings(write.table(alpha_div, file=paste(opts$output,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

print(paste("Name of rarefaction file ", opts$normalize,  sep = ""))
print(paste("Output alpha diversity filename ", opts$output, sep = ""))
