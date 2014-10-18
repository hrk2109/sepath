#!/usr/bin/env python2
from sepath_core import *
from itertools import chain, islice, izip
import re
import optparse
import sys
import numpy
import pysam

if __name__ == "__main__":

    optParser = optparse.OptionParser( 
        usage = "%prog [options] alignment_file annotation_file",
        
        description = \
        "This script takes a paired-end 'alignment_file' in BAM/SAM format and an" +
        "'annotation_file' in GTF/GFF format and for each fragment calculates the distances:" +
        "(1) start of left read to end of right read i.e. insert size," + 
        "(2) start of gene to start of left read," + 
        "(3) end of right read to end of gene",

        epilog = \
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) " +
        "Michigan Center for Translational Pathology (c) 2014 " +
        "Built using 'HTSeq' (%s)." % HTSeq.__version__
    )

    optParser.add_option("--stranded", action="store_true", dest="stranded",
                         default=False, help="turn on strand-specific analysis (fr-firststrand)")

    optParser.add_option("--qc", type="string", dest="qc",
                         default="strict", help="read QC filtering 'strict' or 'loose'")

    optParser.add_option("--n", type="int", dest="n", default=sys.maxint,
                          help="minimum number of fragments to process")

    optParser.add_option("--out", type="string", dest="out", 
                         help="distance output file (dst)"),

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
            sys.stderr.write("warning: missing SO SAM header flag. " +\
                             "Alignment_file should be sorted by queryname (better) or coordinate.\n")

        sys.stderr.write("info: parsing GTF file\n")
        se_ga, se_gm, se_gl, se_gs = parse_gtf(args[1], stranded=opts.stranded)

        sys.stderr.write("info: finding unique sub-exons\n")
        se_unique = unique_subexons(se_ga)

        def count(cargo, seb, split):
            sep = subexonpath(seb, se_unique)
            if sep:
                gene_id = sep[0]
                ses = sep[1]
                try:
                    first_left = ses[0][0]
                    last_left = ses[0][-1]
                    first_right = ses[1][0]
                    last_right = ses[1][-1]
                except IndexError:
                    return
                    
                l_se = (gene_id,) + first_left
                r_se = (gene_id,) + last_right
                first_sep = se_gm[l_se]
                last_sep = se_gm[r_se]
                start = split[1]
                end = split[2]

                enclosed = \
                    (first_sep.start < start < first_sep.end) and \
                    (last_sep.start  <   end < last_sep.end)

                if enclosed:
                    gsg = se_gs[gene_id]
                    l_tail = sum([se_gm[n].length for n in gsg[:gsg.index(l_se)]])
                    r_tail = sum([se_gm[n].length for n in gsg[gsg.index(r_se)+1:]])
                    l_dist = start - se_gm[l_se].start
                    r_dist = se_gm[r_se].end - end
                    l_gap = l_tail + l_dist
                    r_gap = r_tail + r_dist

                    contiguous = (first_right[-1] - last_left[-1] <= 1)
                    if contiguous:
                        ses_lengths = [[se_gm[(gene_id,) + se].length for se in sei] for sei in ses]
                        sep_length = sum(set(chain.from_iterable(ses_lengths)))
                        insert = sep_length - l_dist - r_dist
                    else:
                        insert = -1
                    cargo["gene"].append(gene_id)
                    cargo["data"].append((insert, l_gap, r_gap))

        sys.stderr.write("info: processing BAM file\n")
        cargo = scanBAM(sf, se_ga, count, {"gene":[], "data": []}, opts.progress, opts.qc, "qname", 1, opts.n)
        sys.stderr.write("info: writing dst file '%s'\n" % opts.out)
        out = open(opts.out, "wb") if opts.out else sys.stdout
        for gene, data in izip(cargo["gene"], cargo["data"]):
            out.write("%s\t%d\t%d\t%d\n" % ((gene,) + data))
        out.close()

