#!/usr/bin/env Rscript
library(plyr)
library(stringr)
library(optparse)

loadSamples = function(fns) {
    
    sample_data = ldply(fns, function(fn) {
        df = read.delim(fn, header=FALSE, stringsAsFactors=FALSE)
        colnames(df) = c("gene_id", "gene_length", "count", "unique", "contiguous", "sep_length", "insert_size",
                         "sep1", "sep2")
        df$contiguous = as.logical(toupper(df$contiguous))
        df$sample = str_split(basename(fn), "\\.")[[1]][1]
        return(df)
    })
    return(sample_data)
}


fpkm = function(samples) {
    total = sum(samples$count)
    
    df = ddply(samples, .(gene_id, sample), function(x) {
        gene_length = x$gene_length[1]
        data.frame(
            fpkm = (10**9)*(sum(x$count)/(total*gene_length)),
            gene_length = gene_length
            )
    })
    
    dcast(df, gene_id + gene_length ~ sample, value.var="fpkm")
}

option_list = list(
    
    make_option(c("--level"), type="character", default="gene",
                help="output file name [default: %default]",
                metavar="level"),
    
    )

parser = OptionParser("%prog [options] sep_files",

    description=c(
        "This script takes one or more 'sep_files' and estimates the gene-level FPKMs."),

    epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2014 "),

    option_list=option_list
    )

opt =  parse_args(parser, positional_arguments=TRUE)

if (length(opt$args) < 1) {
    write("sepath-fpkm.r requires at least one 'sep' file.", stderr())
    write("Incorrect number of required positional arguments.\n", stderr())
    print_help(parser)
    quit("no", 1)
}

samples = loadSamples(opt$args)


if (opt$options$out == "stdout") {
##     write.table(df, stdout(), sep="\t", quote=FALSE, col.names=FALSE)
## } else {
##     write.table(df, opt$options$out, sep="\t", quote=FALSE, col.names=FALSE)
## }
