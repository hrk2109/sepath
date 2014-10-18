#!/usr/bin/env Rscript
library(optparse)

parser = OptionParser("%prog [options] gff_file dst_file",

    description=c(
        "This script takes an 'annotation file' in GFF format and a 'distance file'",
        "in DST format and calculates the 5'/3' bias."
        ),
    
    epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2014 "
        ),

    option_list = list(
        make_option("--out", type="character", default="stdout",
                help="output file name [default: %default]", metavar="out")
        ))

opt = parse_args(parser, positional_arguments=TRUE)
if (length(opt$args) != 2) {
    
    write("error: sepath_bias.r requires an annotation 'gff_file' and a distance 'dst_file'", stderr())
    print_help(parser)
    quit("no", 1)
    
}
