#!/usr/bin/env python2
from sepath_core import *
from itertools import chain, islice
import optparse
import sys
import numpy
import pysam

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

    optParser.add_option("--n", type="int", dest="n", default=100000,
                          help="minimum number of fragments to process")

    optParser.add_option("--out", type="string", dest="out", 
                         help="insert size output file (tsv)"),

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
        se_ga, se_gm, se_gl = parse_gtf(args[1], stranded=opts.stranded)

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
                    
                start = split[1]
                end = split[2]

                contiguous = (first_right[-1] - last_left[-1] <= 1)

                first_sep = se_gm[(gene_id,) + first_left]
                last_sep = se_gm[(gene_id,) + last_right]

                enclosed = \
                    (first_sep.start < start < first_sep.end) and \
                    (last_sep.start  <   end < last_sep.end)

                if contiguous and enclosed:
                    ses_lengths = [[(se_gm[(gene_id,) + se]).length for se in sei] for sei in ses]
                    sep_length = sum(set(chain.from_iterable(ses_lengths)))
                    insert = sep_length - (start - se_gm[(gene_id,) + first_left].start) - \
                                          (se_gm[(gene_id,) + last_right].end - end)
                    cargo.append(insert)

        sys.stderr.write("info: processing BAM file\n")
        cargo = scanBAM(sf, se_ga, count, [], opts.progress, opts.qc, "qname", 1, opts.n)

        sys.stderr.write("info: writing tsv file '%s'\n" % opts.out)
        out = open(opts.out, "wb") if opts.out else sys.stdout
        
        std = numpy.std(cargo)
        avg = numpy.mean(cargo)
        med = numpy.median(cargo)
        n = len(cargo)
        out.write("n\t%s\nmean\t%.2f\nstd\t%.2f\nmed\t%.2f\n" % (n, avg, std, med))
        out.close()

        #  for \
        #            (chr, start, end, strand), count in cargo.iteritems()]
        # isize = "%d:%d:%s:%s" % (numpy.mean(inserts), numpy.std(inserts), 
        #                          min(inserts), max(inserts))


        # sys.stderr.write("info: calculating sub-exon paths\n")
        # sep_cargos = seb2sep(seb_cargos, se_unique)

        # sys.stderr.write("info: counting sub-exon paths\n")
        # sep_counts = count_subexonpaths(sep_cargos)
        # sep_lengths = measure_subexonpaths(sep_cargos, se_gm)


