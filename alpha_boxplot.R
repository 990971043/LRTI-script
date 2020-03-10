#!/usr/bin/env Rscript

# Functions: Alpha boxplot
options(warn = -1) # Turn off warning

suppressWarnings(suppressMessages(library(amplicon)))


# 读取
alpha_div = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
p = alpha_boxplot(alpha_div, index = opts$alpha_index, metadata, groupID = opts$group)

# Saving figure
ggsave(paste0(opts$output,"/alpha_boxplot_",opts$alpha_index,".pdf"), p, width = opts$width, height = opts$height, units = "mm")
