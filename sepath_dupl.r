#!/usr/bin/env Rscript
library(plyr)
library(stringr)
library(optparse)

loadSamples = function(fns) {
    sample_data = ldply(fns, function(fn) {
        df = read.delim(fn, header=FALSE, stringsAsFactors=FALSE)
        colnames(df) = c("gene_id", "sep1", "sep2", "count", "unique")
        df$sample = str_split(basename(fn), "\\.")[[1]][1]
        return(df)
    })
    return(sample_data)
}

option_list = list(
    
    make_option(c("--level"), type="character", default="gene",
                help="output file name [default: %default]",
                metavar="out"),
    
    )

parser = OptionParser("%prog [options] sep_files",

    description=c(
        "This script takes one or more 'sep_files' and estimates the percentage of duplicate reads"),

    epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2014 "),

    option_list=option_list
    )

opt =  parse_args(parser, positional_arguments=TRUE)

if (length(opt$args) < 1) {
    write("sepath-dup.r requires at least one 'sep' file.", stderr())
    write("Incorrect number of required positional arguments.\n", stderr())
    print_help(parser)
    quit("no", 1)
}


samples = loadSamples(opt$args, groups)
sepdiv = sepDivergences(samples)

df = data.frame(SEPDIV=as.matrix(unlist(sepdiv)))

if (opt$options$out == "stdout") {
    write.table(df, stdout(), sep="\t", quote=FALSE, col.names=FALSE)
} else {
    write.table(df, opt$options$out, sep="\t", quote=FALSE, col.names=FALSE)
}
