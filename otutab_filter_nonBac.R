
otutab_dir = "result/raw/otutab.txt"
sintax_dir = "result/raw/otus.sintax"

output_dir = "result/otutab.txt"
stat_dir = "result/raw/otutab_nonBac.txt"

otutab = read.table(otutab_dir, header=T, row.names=1, sep="\t", comment.char="")
sintax = read.table(sintax_dir, header=F, row.names=1, sep="\t", comment.char="")

print(paste0("Input feature table is ", otutab_dir))
print(paste0("Input sintax taxonomy table is ", sintax_dir))

total_reads = colSums(otutab)

idx = grepl("Bacteria|Archaea", sintax$V4, perl = T)
nonspecific = sintax[!idx,]
nonspecific_reads = colSums(otutab[rownames(nonspecific),])
sintax = sintax[idx,]

idx = grepl("Chloroplast", sintax$V2, perl = T)
chloroplast = sintax[idx,]
chloroplast_reads = colSums(otutab[rownames(chloroplast),])
sintax = sintax[!idx,]

idx = grepl("Mitochondria", sintax$V2, perl = T)
mitochondria = sintax[idx,]
mitochondria_reads = colSums(otutab[rownames(mitochondria),])
sintax = sintax[!idx,]

# otutab = otutab[rownames(sintax),]
idx = rownames(otutab) %in% rownames(sintax)
otutab = otutab[idx,]
idx = order(rowSums(otutab), decreasing = T)
otutab = otutab[idx,]
filtered_reads = colSums(otutab)

write.table(rbind(nonspecific, chloroplast, mitochondria), file=stat_dir, append = F, sep="\t", quote=F, row.names=T, col.names=F, eol = "\n")
df = as.data.frame(cbind(total_reads, nonspecific_reads, chloroplast_reads, mitochondria_reads, filtered_reads))
suppressWarnings(write.table(paste("\nSampleID\t",  sep=""), file=stat_dir, append = T, sep="\t", quote=F, row.names=F, col.names=F, eol = ""))
suppressWarnings(write.table(df, file=paste(stat_dir, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

write.table(paste("#OTUID\t",  sep=""), file=paste(output_dir, sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
suppressWarnings(write.table(otutab, file=paste(output_dir, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

print(paste0("Summary of samples size in final feature table: "))
summary(filtered_reads)

print(paste0("Onput feature table is ", output_dir))
print(paste0("Detail and statistics in ", stat_dir))
