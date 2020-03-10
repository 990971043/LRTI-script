
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)


suppressWarnings(suppressMessages(library(amplicon)))
format2stamp(otutab, tax, opts$threshold, opts$output)
