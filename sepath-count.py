#!/usr/bin/env python2
from collections import defaultdict
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
                fh.write("\t".join(map(str, [iv.chrom, iv.start, iv.end, seg_id, 0, iv.strand, "\n"])))

def write_json(x, fn):
    with open(fn, "wb") as fh:
        json.dump(x, fh)

def write_sep(sep_counts, fn):
    with open(fn, "wb") as fh:
        for gene_id, seps_counts in sep_counts["assigned"].iteritems():
            for seps, count in seps_counts.iteritems():
                ses1 = "-".join(["%s:%s" % se for _, se in seps[0]])
                ses2 = "-".join(["%s:%s" % se for _, se in seps[1]])
                line = "\t".join([gene_id, ses1, ses2, str(count), "\n"])
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
        subexonbag = []
        rngs = mapped_ranges(read)
        for rng in rngs:
            rng_iv = HTSeq.GenomicInterval(read_chr, rng[0], rng[1], read_strand)
            for iv, subexons in ga[rng_iv].steps():
                subexonbag.extend(subexons)
        return tuple(sorted(subexonbag))

    return tuple(sorted((bag(read, read_chr, read_strand, ga), bag(mate, mate_chr, mate_strand, ga))))




def count_subexonbags(sam, ga, verbose=9223372036854775783):
    sf = pysam.Samfile(sam)
    reads = {}
    counts = defaultdict(int)
    i = 0
    for read in sf:
        read_strand = "-" if read.is_reverse else "+"
        read_chr = sf.getrname(read.tid)
        read_id = (read_chr, read.pos, read_strand, read.qname, read.is_read1)
        mate_strand = "-" if read.mate_is_reverse else "+"
        mate_chr = sf.getrname(read.mrnm)
        mate_id = (mate_chr, read.mpos, mate_strand, read.qname, read.is_read2)
        mate = reads.pop(mate_id, None)
        if mate:
            # this read was second in alignment, mate was first in alignment
            # for stranded libraries we fix the order of the reads and invert 
            # the strand of the second read.
            if read.is_read1:
                counts[subexonbag(read, read_chr, read_strand, mate, read_chr, mate_strand, ga)] +=1
            else:
                counts[subexonbag(mate, mate_chr, mate_strand, read, mate_chr, read_strand, ga)] +=1

            # log
            i+=1
            if i % verbose == 0:
                sys.stderr.write("%d fragments processed.\n" % i)
                sys.stderr.write("%d reads in queue.\n" % len(reads))

        else:
            # first read in alignment
            if read.is_proper_pair and not read.is_secondary:
                reads[read_id] = read
            else:
                print read

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
                gene_subexonbags[gene_id][i].append(subexon_id)
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


def main():

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

    optParser.add_option("-o", "--out", type="string", dest="out",
                         default="", help="sub-exon path output file (tsv)"),


    optParser.add_option("-p", "--sep_json", type="string", dest="sep_json",
                         help="full sub-exon path output file (json)"),

    optParser.add_option("-b", "--seb_json", type="string", dest="seb_json",
                         help="full sub-exon bag output file (json)"),

    optParser.add_option("-e", "--eattr", type="string", dest="eattr",
                         default="exon_id", help="GFF attribute to be used as exon id (default, " +
                         "suitable for Ensembl GTF files: exon_id)"),
         
    optParser.add_option("-g", "--gattr", type="string", dest="gattr",
                         default="gene_id", help="GFF attribute to be used as gene id (default, " +
                         "suitable for Ensembl GTF files: gene_id)"),

    optParser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                          help="run-time messages printed to stderr")

    optParser.add_option("-r", "--progress", type="int", dest="progress", default=100000,
                          help="progress on BAM processing printed every n lines")

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()
    
    if len(args) != 2:
        optParser.print_help()
        #sys.exit(1)

    print opts
    print args
      
if __name__ == "__main__":
    main()




# c_bam = "test/bam/coordinate.bam"
# n_bam = "test/bam/name.bam"

# gtf = "test/gencode.v19.annotation.gtf.gz"

# se_ga = parse_gtf(gtf, stranded=False)
# write_bed(se_ga, "se_ga_false.bed")

# se_unique = unique_subexons(se_ga)
# c_seb_counts = count_subexonbags(c_bam, se_ga, 1000) 
# n_seb_counts = count_subexonbags(n_bam, se_ga, 1000) 
# print c_seb_counts == n_seb_counts

# b1 = count_subexonbags("test/sam/1.bam", se_ga) 
# p1 = counts_subexonpaths(b1, se_unique)
# b5 = count_subexonbags("test/sam/5.bam", se_ga) 
# p5 = counts_subexonpaths(b5, se_unique)

# b9 = count_subexonbags("test/sam/9.bam", se_ga) 
# p9 = counts_subexonpaths(b9, se_unique)
# b13 = count_subexonbags("test/sam/13.bam", se_ga) 
# p13 = counts_subexonpaths(b9, se_unique)
# b233_1 = count_subexonbags("test/sam/233_1.bam", se_ga) 
# p233_1 = counts_subexonpaths(b233_1, se_unique)
# b233_2 = count_subexonbags("test/sam/233_2.bam", se_ga) 
# p233_2 = counts_subexonpaths(b233_2, se_unique)


# sep_counts = counts_subexonpaths(c_seb_counts, se_unique)
# se_ga_stranded = parse_gtf(gtf, stranded=True)

