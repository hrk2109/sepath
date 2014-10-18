#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(rtracklayer)

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

write("info: loading DST file", stderr())
dst = fread("test/bias/test.dst")
setnames(dst, c("gene_id", "insert", "start", "end"))
write("info: loading GFF file", stderr())
gff = import("test/bias/gencode.v19.annotation.gtf")

e = gff[gff$type == "exon"]
e_g = split(e, e$gene_id)
gl = data.table(gene_id = names(e_g),
                length = sapply(width(reduce(e_g)), sum))
g = gff[gff$type == "gene"]
ga = data.table(gene_id = as.character(g$gene_id),
    gene_name = as.character(g$gene_name),
    gene_type = as.character(g$gene_type),
    strand = as.character(strand(g)))
ga = merge(ga, gl, by="gene_id")
dst = merge(dst, ga, by="gene_id")

dst$rel_start = ifelse(dst$strand == "+", dst$start, dst$end)
dst$rel_end = ifelse(dst$strand == "+", dst$end, dst$start)


sel = dst[gene_type == "protein_coding"]
sel[,size:=cut(length, c(0,1000,2000,5000,10000,1e6)),]
sel[,bias:=rel_start/length,]

plt = ggplot(sel) + aes(x=bias, group=size, color=size) + geom_density(adjust=10) + theme_bw()
ggsave("t.png", plt)
























a = dst[gene_name == "JUNB"]
tbl = a[insert+start+end == length]

library(ggplot2)
plt = ggplot(tbl) + theme_bw() +
    aes(x=start, y=insert) +
    geom_point()

