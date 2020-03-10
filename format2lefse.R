
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F) 

metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="")

suppressWarnings(suppressMessages(library(amplicon)))
format2lefse(otutab, tax, metadata, opts$threshold, opts$group, opts$output)

