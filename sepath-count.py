#!/usr/bin/env python2
from collections import defaultdict
import re
import optparse
import sys
import json
import pysam
import HTSeq

# "M":0, # aligned, add
# "I":1, # insert, skip
# "D":2, # deleted, add
# "N":3, # skipped, add
# "S":4, # soft-clipped, skip
# "H":5, # hard-clipped, skip
# "P":6, # padded, skip
# "=":7, # base match, add
# "X":8, # base mis-match, add
COUNT_ON_REF = set((0, 2, 3, 7, 8))

def parse_gtf(gtf, stranded=False):

    gtf_file = HTSeq.GFF_Reader(gtf, end_included=True)
    ges_ga = HTSeq.GenomicArrayOfSets("auto", stranded=stranded) # gene-exon genomic array
    g_ivs = {}

    for feature in gtf_file:
        if feature.type in ("gene", "exon"):
            gene_id = feature.attr["gene_id"]
            if feature.type == "exon":
                ges_ga[feature.iv] += (gene_id, feature.attr["exon_id"])
            elif feature.type == "gene":
                g_ivs[gene_id] = feature.iv
    
    se_ga = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    for gene_id, g_iv in g_ivs.iteritems():
        i, j = 0, 0
        for ge_iv, ge_ids in ges_ga[g_iv].steps():
            g_ids = [g_id for g_id, e_id in ge_ids]
            if gene_id in g_ids:
                se_ga[ge_iv] += (gene_id, (i, j))
                j += 1
            elif j:
                i += 1
                j = 0

    return se_ga

def write_bed(se_ga, fn):
    with open(fn, "wb") as fh:
        for iv, ses in se_ga.steps():
            for gene_id, e_rank in ses:
                seg_id = "%s_%s:%s" % (gene_id, e_rank[0], e_rank[1])
                fh.write("\t".join(map(str, [iv.chrom, iv.start, iv.end, seg_id, 0, 
                                             iv.strand, "\n"])))

def write_sep(sep_counts, fh):
    for gene_id, seps_counts in sep_counts["assigned"].iteritems():
        for seps, count in seps_counts.iteritems():
            ses1 = "-".join(["%s:%s" % se for se in seps[0]])
            ses2 = "-".join(["%s:%s" % se for se in seps[1]])
            line = "\t".join([gene_id, ses1, ses2, str(count)]) + "\n"
            fh.write(line)

def mapped_ranges(read):
    ranges = []
    offset = read.pos
    for (op, length) in read.cigar:
        if op in COUNT_ON_REF: 
            if op == 0:
                ranges.append((offset, offset+length))
            offset += length
    return ranges

def subexonbag(read, read_chr, read_strand, mate, mate_chr, mate_strand, ga):

    def bag(read, read_chr, read_strand, ga):
        subexonbag = set() # rngs mapping to the same sub-exon will be merged
        rngs = mapped_ranges(read)
        for rng in rngs:
            rng_iv = HTSeq.GenomicInterval(read_chr, rng[0], rng[1], read_strand)
            for iv, subexons in ga[rng_iv].steps():
                subexonbag.update(subexons)
        return tuple(sorted(subexonbag))

    return tuple(sorted((bag(read, read_chr, read_strand, ga), bag(mate, mate_chr, mate_strand, ga))))

def count_subexonbags(sf, ga, progress, qc):
    #
    reads = {}
    counts = defaultdict(int)
    i_pair = 0
    i_pass = 0
    for i_any, read in enumerate(sf):
        # LOG
        if (i_any + 1) % progress == 0:
            sys.stderr.write("reads: %d (processed) %d (pass QC) %d (paired) %d (unpaired)\n" %
                             (i_any + 1, i_pass, i_pair*2, len(reads)))
        # QC
        if qc == "strict":
            if not read.is_proper_pair:
                continue
        elif qc == "loose":
            if read.is_unmapped or read.mate_is_unmapped:
                continue
        if read.is_secondary:
            continue
        i_pass += 1

        read_strand = "-" if read.is_reverse else "+"
        read_chr = sf.getrname(read.tid)
        read_id = (read_chr, read.pos, read_strand, read.qname, read.is_read2)
        mate_strand = "-" if read.mate_is_reverse else "+"
        mate_chr = sf.getrname(read.mrnm)
        mate_id = (mate_chr, read.mpos, mate_strand, read.qname, read.is_read1)
        mate = reads.pop(mate_id, None)

        if mate:
            # read was second in alignment, mate was first in alignment (irrelevant now)
            # for stranded libraries:
            #  - always process as (read1, read2)
            #  - read2 is mapped to the + strand
            if read.is_read1:
                counts[subexonbag(read, read_chr, mate_strand, mate, mate_chr, mate_strand, ga)] +=1
            else:
                counts[subexonbag(mate, mate_chr, read_strand, read, read_chr, read_strand, ga)] +=1
            # log
            i_pair+=1
        else:
            # first read in alignment
            reads[read_id] = read
    counts = dict(counts)
    return counts

def unique_subexons(ga):
    se_uniq = set()
    for iv, subexons in ga.steps():
        if len(subexons) == 1:
            se_uniq.update(subexons)
    return se_uniq
 
def counts_subexonpaths(seb_counts, se_uniq):

    NULL_SEP = ((),())

    uniq_counts = defaultdict(lambda: defaultdict(int)) # uniquely assigned
    isep_counts = defaultdict(lambda: defaultdict(int)) # inseparable
    uacc_counts = defaultdict(lambda: defaultdict(int)) # unaccountable
    uass_counts = {NULL_SEP:0}                          # unassigned

    for seb, seb_count in seb_counts.iteritems():
        ## split subexons by genes count unique ones
        gene_subexonbags = defaultdict(lambda:[[],[]]) # keyed by gene
        gene_subexonuniq = defaultdict(int) # keyed by gene
        for i, read_exons in enumerate(seb):
            for subexon_id in read_exons:
                gene_id = subexon_id[0]
                gene_subexonbags[gene_id][i].append(subexon_id[1])
                gene_subexonuniq[gene_id] += (subexon_id in se_uniq)

        uniqgene_ids = [g_id for g_id in gene_subexonuniq if gene_subexonuniq[g_id]]
        if len(uniqgene_ids) == 1:
            # 
            gene_id = uniqgene_ids.pop()
            seb = gene_subexonbags[gene_id]
            sep = tuple(sorted([tuple(sorted(seb[0])), 
                                tuple(sorted(seb[1]))]))
            uniq_counts[gene_id][sep] += seb_count
        elif len(uniqgene_ids) == 0 and len(gene_subexonbags):
            # inseperable - mapped only to non-unique positions
            gene_seps = []
            for gene_id in gene_subexonbags:
                seb = gene_subexonbags[gene_id]
                sep = tuple(sorted([tuple(sorted(seb[0])), 
                                    tuple(sorted(seb[1]))]))
                gene_seps.append(sep)
            isep_counts[tuple(sorted(gene_subexonbags))][tuple(sorted(gene_seps))] += seb_count
        elif len(uniqgene_ids):
            # unaccountable - mapped to multiple unique positions
            gene_seps = []
            for gene_id in uniqgene_ids:
                seb = gene_subexonbags[gene_id]
                sep = tuple(sorted([tuple(sorted(seb[0])), 
                                    tuple(sorted(seb[1]))]))
                gene_seps.append(sep)
            uacc_counts[tuple(sorted(uniqgene_ids))][tuple(sorted(gene_seps))] += seb_count
        else:
            # unassigned - mapped outside of counting regions
            uass_counts[NULL_SEP] += seb_count
    
    sep_counts = {"assigned":uniq_counts, "unassigned":uass_counts, 
                  "inseparable":isep_counts, "unaccountable":uacc_counts}   
    return sep_counts


if __name__ == "__main__":

    optParser = optparse.OptionParser( 
        usage = "%prog [options] alignment_file annotation_file",
        
        description = \
        "This script takes a paired-end 'alignment_file' in BAM/SAM format and a" +
        "'annotation_file' in GTF/GFF format and counts how many times a fragment was" +
        "mapped to a specific order of sub-exons a.k.a its sub-exon path",
        
        epilog = \
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) " +
        "Michigan Center for Translational Pathology (c) 2014 " +
        "Built using 'HTSeq' (%s)." % HTSeq.__version__
    )

    optParser.add_option("--stranded", action="store_true", dest="stranded",
                         default=False, help="turn on strand-specific analysis (fr-firststrand)")

    optParser.add_option("--qc", type="string", dest="qc",
                         default="strict", help="read QC filtering 'strict' or 'loose'")

    optParser.add_option("--out", type="string", dest="out", 
                         help="sub-exon path output file (tsv)"),

    optParser.add_option("--se_bed", type="string", dest="se_bed",
                         help="derived sub-exon annotation (bed"),

    optParser.add_option("--sep_json", type="string", dest="sep_json",
                         help="full sub-exon path output file (json)"),

    optParser.add_option("--seb_json", type="string", dest="seb_json",
                         help="full sub-exon bag output file (json)"),

    optParser.add_option("--eattr", type="string", dest="eattr",
                         default="exon_id", help="GFF attribute to be used as exon id (default, " +
                         "suitable for Ensembl GTF files: exon_id)"),
         
    optParser.add_option("--gattr", type="string", dest="gattr",
                         default="gene_id", help="GFF attribute to be used as gene id (default, " +
                         "suitable for Ensembl GTF files: gene_id)"),

    optParser.add_option("--verbose", action="store_true", dest="verbose",
                          help="run-time messages printed to stderr")

    optParser.add_option("--progress", type="int", dest="progress", default=100000,
                          help="progress on BAM processing printed every n lines")

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()
    
    if len(args) != 2:
        optParser.print_help()
        sys.exit(1)

    with pysam.Samfile(args[0]) as sf:
        
        try:
            order = re.search("SO:(.*)", sf.text).groups()[0]
        except Exception, e:
            order = None
        if not order in ("queryname", "coordinate"):
            sys.stderr.write("error: alignment_file should be sorted by queryname (better) or coordinate.\n")
            sys.exit(1)
    
        sys.stderr.write("parsing GTF file\n")
        se_ga = parse_gtf(args[1], stranded=opts.stranded)

        if opts.se_bed:
            sys.stderr.write("writing sub-exon BED file\n")
            write_bed(se_ga, opts.se_bed)

        sys.stderr.write("counting sub-exon bags\n")
        seb_counts = count_subexonbags(sf, se_ga, opts.progress, opts.qc)

        sys.stderr.write("finding unique sub-exons\n")
        se_unique = unique_subexons(se_ga)

        sys.stderr.write("counting sub-exon paths\n")
        sep_counts = counts_subexonpaths(seb_counts, se_unique)

        out = open(opts.out, "wb") if opts.out else sys.stdout
        write_sep(sep_counts, out)

