#!/usr/bin/env Rscript
library(plyr)
library(stringr)
library(optparse)
library(data.table)
# library(parallel)
# library(doMC)
# registerDoMC(detectCores())
library(MultinomialCI)

merge.list =function (x, y, ...) {
    if (length(x) == 0) 
        return(y)
    if (length(y) == 0) 
        return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) 
        x[names(y)[which(i)]] = y[which(i)]
    x
}

shannonEntropy = function(freq, base=2) {
    tmp = log(freq, base)
    tmp[is.infinite(tmp)] = NA
    Hs = -rowSums(tmp*freq, na.rm=TRUE)
    return(Hs)
}

JSD = function(freq, weights=rep(1/nrow(freq), nrow(freq)), base=2) {
    wmean = matrix(colSums(sweep(freq, 1, weights, "*")), nrow=1)
    H_avg_freq = shannonEntropy(wmean, base)    
    H_freq = shannonEntropy(freq, base)
    avg_H_freq = weighted.mean(H_freq, weights)
    H_avg_freq - avg_H_freq
}

JSDg = function(freq, groups, weights=rep(1/nrow(freq), nrow(freq)), base=2) {
    ## this the first term of the JSD
    wmean = matrix(colSums(sweep(freq, 1, weights, "*")), nrow=1)
    H_avg_freq = shannonEntropy(wmean, base)
    ## now rather than calculateing the avg shannon entroypy
    ## of each of the freqs we calculate mean frequency for
    ## each group and for it the shannon entropy the group Hs
    ## are weighted by the group size to give the second term
    ## of the JSD
    avg_H_freq = 0
    for (group in unique(groups)) {
        group_sel = (groups == group)
        freq_g = freq[group_sel, , drop=FALSE]
        weights_g = weights[group_sel]
        weights_freq_g = weights_g / sum(weights_g)
        wmean_g = matrix(colSums(sweep(freq_g, 1, weights_freq_g, "*")), nrow=1)
        H_avg_freq_g = shannonEntropy(wmean_g, base)
        avg_H_freq = avg_H_freq + sum(weights_g) * H_avg_freq_g
    }
    H_avg_freq - avg_H_freq
}

sepDivergence = function(sep_count, groups, alpha, sepcut, rowcut, colcut) {
    ## For each sample determine which seps are present above 'sepcut' in at least 'colcut'
    ## samples. And which samples contain at least 'rowcut' 
    ## at a 'alpha' false positive rate
    sep_count = as.matrix(sep_count)
    tmp = alply(sep_count, 1, function(row) {
        if (sum(row) >= length(row)) {
            multinomialCI(row, alpha)[,1]
        } else {
            rep(0, length(row))
        }
    }, .parallel=TRUE)
    low = splat(rbind)(tmp)
    ## filter
    lowcut = (low > sepcut)
    col_pass = colSums(lowcut) > (colcut * nrow(low))
    row_pass = rowSums(lowcut) > (rowcut * ncol(low))
    sep_count_pass = sep_count[row_pass,col_pass,drop=FALSE]
    groups_pass = groups[row_pass]
    ## transform rows of counts into rows of frequencies
    gene_sample_coverage = rowSums(sep_count_pass)
    sep_freq = sweep(sep_count_pass, 1, gene_sample_coverage, "/")
    ## calculate groupwise JSD
    JSDg(sep_freq, groups_pass)
}

loadSamples = function(fns, groups) {
    sample_data = ldply(fns, function(fn) {
        df = read.delim(fn, header=TRUE, stringsAsFactors=FALSE)
        df$sample = str_split(basename(fn), "\\.")[[1]][1]
        return(df)
    })
    samples = sapply(str_split(basename(fns), "\\."), "[", 1)
    sample_groups = data.frame(group=groups, stringsAsFactors=FALSE)
    rownames(sample_groups) = samples
    return(list(sample_data, sample_groups))
}

sepCounts = function(df) {
    samples = unique(df$sample)
    sep_counts = dlply(df, c("gene_id"), function(gdf) {
        ggdf = dlply(gdf, c("sample"), function(sdf) {
            gd = data.frame(matrix(sdf$total, ncol=nrow(sdf)))
            colnames(gd) = paste(sdf$sep1, sdf$sep2, sep="--")
            rownames(gd) = sdf$sample[[1]]
            return(gd)
        })
        res = splat(rbind.fill)(ggdf)
        res[is.na(res)] = 0
        rownames(res) = unlist(lapply(ggdf, row.names))
        res[setdiff(samples, rownames(res)),] = 0
        res = res[samples,]
        return(res)
    }, .parallel=TRUE)
    return(sep_counts)
}           

sepDivergences = function(samples, alpha=0.05, sepcut=0, colcut=1/3, rowcut=1/10) {
    sep_counts = sepCounts(samples[[1]])
    groups = samples[[2]]
    llply(sep_counts, function(sep_count) {
        gene_groups = groups[rownames(sep_count),]
        sepDivergence(sep_count, gene_groups, alpha=alpha, sepcut=sepcut, rowcut=rowcut, colcut=colcut)
    })
}
    
option_list = list(
    
    make_option(c("--groups"), type="character",
                help="comma-delimited groups associated with the samples [default: A,B,C...]",
                metavar="groups"),

    make_option(c("--out"), type="character", default="stdout",
                help="output file name [default: %default]",
                metavar="out"),

    make_option(c("--alpha"), type="numeric", default=0.05,
                help="confidence level [default: %default]",
                metavar="alpha"),

    make_option(c("--sepcut"), type="numeric", default=0,
                help="sep cut-off [default: %default]",
                metavar="cutoff"),

    make_option(c("--colcut"), type="numeric", default=0,
                help="col cut-off [default: %default]",
                metavar="cutoff"),

    make_option(c("--rowcut"), type="numeric", default=0,
                help="row cut-off [default: %default]",
                metavar="cutoff"),

    
    make_option(c("--fract"), type="numeric", default=0.33,
                help="fract [default: %default]",
                metavar="fract")
    
    )

parser = OptionParser("%prog [options] sep_files",

    description=c(
        "This script takes two or more 'sep_files' belonging to two or more 'groups' and estimates",
        "for each gene differences in splice-exon path distributions, which is a proxy for alternative",
        "splicing. Specifically, it calculates the excess Jensen-Shannon divergence"),

    epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2014 "),

    option_list=option_list
    )

opt =  parse_args(parser, positional_arguments=TRUE)

if (length(opt$args) < 2) {
    write("sepath-diff.r requires at least two 'sep' files.", stderr())
    write("Incorrect number of required positional arguments.\n", stderr())
    print_help(parser)
    quit("no", 1)
}

if (is.null(opt$options$groups)) {
    groups = LETTERS[1:length(opt$args)]
    write(paste("warning: groups not specified assigning each sample to separate groups:",
                paste(groups, collapse=","), "\n"), stderr())
} else {
    groups = str_split(opt$options$groups, ",")[[1]]
}

if (length(opt$args) != length(groups)) {
    write("Number of 'groups' does not match the number of 'sep' files.\n", stderr())
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
