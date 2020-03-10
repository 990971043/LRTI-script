
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# head(otutab)

design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F) 

design$group=design[,opts$group]

norm = t(otutab)/colSums(otutab,na=T)*100
# rowSums(norm)
idx = colMeans(norm) > opts$thre
HA = norm[,idx]
# dim(HA)
# rowSums(HA)

merge=cbind(HA, design[,c("group"),drop=F])
HA_group_mean = merge %>% group_by(group) %>% summarise_all(mean)
HA_t = as.data.frame(t(HA_group_mean))
HA_t = cbind(HA_t,  c("All", colMeans(norm)))
# HA_t$V4 = c("All", colMeans(norm))
rownames(HA_t)[1] = "OTUID"
write.table(HA_t, file=paste(opts$output, "", sep = ""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=F)

