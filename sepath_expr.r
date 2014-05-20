#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(optparse)

strReverse = function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

loadSepTbl = function(fns) {
    tbls = lapply(fns, function(fn) {
        df = fread(fn)
        df$sample = str_split(basename(fn), "\\.")[[1]][1]
        return(df)
    })
    tbl = do.call("rbind", tbls)
    return(tbl)
}

loadBed = function(bedfn) {
    bed = read.delim(bedfn, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(bed) = c("chr", "start", "end", "name", "score", "strand")
    return(bed)
}

seTbl = function(sep, bed) {
    splits = str_split(strReverse(bed$name), "_", 2)    
    gse = data.table(gene_id=strReverse(sapply(splits, "[" ,2)),
                         sep=strReverse(sapply(splits, "[" ,1)))
    gse$length = bed$end - bed$start
    setkey(gse, "gene_id", "sep")
    seps = lapply(str_split(paste(sep$sep1, sep$sep2, sep="-"), "-"), function(x) unique(x[x!=""]))
    unrolled_seps = unlist(seps)
    unrolled_idx = rep(1:nrow(sep), sapply(seps, length))
    unrolled_sep = sep[unrolled_idx]
    unrolled_sep$sep = unrolled_seps
    counts = unrolled_sep[, list(total=sum(total),unique=sum(unique)), by=list(gene_id, sep, sample)]
    setkey(counts, "gene_id", "sep", "sample")
    se = merge(counts, gse, all.y=TRUE)
    se$sample = unique(counts$sample)
    se[is.na(se)] = 0
    return(se)
}

seLength = function(se, quant=0.5, rl=100) {
    sel = copy(se) # copy
    sel[,nfpb:=(total+1)/(length+(2*rl))]
    sel[,se_weight:=pmin(1.0, nfpb/quantile(nfpb, quant)), by=list(gene_id, sample)]
    sel[,elength:=length*se_weight,by=list(gene_id, sample)]
    setkey(sel, gene_id, sample)
    return(sel)
}

fpkm = function(sep, sel) {
    lent = sel[,list(length=as.numeric(sum(length)), elength=sum(elength)),by=list(gene_id, sample)]
    gent = sep[, list(gene_total=sum(total)), by=list(sample, gene_id)]
    setkey(lent, gene_id, sample)
    setkey(gent, gene_id, sample)
    gent[, sample_total:=sum(gene_total), by=list(sample)]
    tmp = merge(gent, lent)
    tmp[, fpkm:=10**9 * gene_total / (sample_total *  length),]
    tmp[,efpkm:=10**9 * gene_total / (sample_total * elength),]
    return(tmp)
}

parser = OptionParser("%prog [options] bed_file cnt_file",

    description = c(
        "This script takes a sepath_count output BED file and a CNT file and estimates expression (e)FPKM"),
    
    epilogue = c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2014 "),

    option_list = list(
        make_option(c("--out"), type="character", default="stdout",
                help="output file name [default: %default]",
                metavar="out"),
        make_option(c("--quant"), type="numeric", default=0.5,
                help="Exon coverage penalty quantile [default: %default]",
                    metavar="quant"),
        make_option(c("--rl"), type="numeric", default=100,
                help="Averaged read length [default: %default]",
                metavar="quant"),
        make_option(c("--level"), type="character", default="gene",
                help="output file name [default: %default]",
                metavar="level")
        )
    )

opt = parse_args(parser, positional_arguments=TRUE)

if (length(opt$args) != 2) {
    write("sepath_expr.r requires a sepath_count 'BED' file and a 'CNT' file.", stderr())
    print_help(parser)
    quit("no", 1)
}

write("info: loading BED file", stderr())
bed = loadBed(opt$args[1])
write("info: loading CNT file", stderr())
sep = loadSepTbl(opt$args[2])
write("info: computing exon counts", stderr())
se = seTbl(sep, bed)
write("info: computing corrected exon lengths", stderr())
sel = seLength(se, quant=opt$options$quant, rl=opt$options$rl)
write("info: computing FPKM", stderr())
res = fpkm(sep, sel)

df = res[,c("gene_id", "fpkm", "efpkm"),with=FALSE]
if (opt$options$out == "stdout") {
    write.table(df, stdout(), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
} else {
    write.table(df, opt$options$out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}



