import sys
import HTSeq
from collections import defaultdict

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
    se_gm = {}
    se_gl = defaultdict(int)
    
    for gene_id, g_iv in g_ivs.iteritems():
        # exon, subexon, n
        i, j, n = 0, 0, 0
        for ge_iv, ge_ids in ges_ga[g_iv].steps():            
            g_ids = [g_id for g_id, e_id in ge_ids]
            if gene_id in g_ids:
                se_gl[gene_id] += ge_iv.length
                se = (gene_id, i, j, n)
                se_ga[ge_iv] += se
                se_gm[se] = ge_iv
                j += 1
                n += 1
            elif j:
                i += 1
                n += 1
                j = 0

    return se_ga, se_gm, se_gl

def unique_subexons(ga):
    se_uniq = set()
    for iv, subexons in ga.steps():
        if len(subexons) == 1:
            se_uniq.update(subexons)
    return se_uniq


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

def subexonpath(seb, se_uniq):
    
    gene_subexonbags = defaultdict(lambda:[[],[]]) # keyed by gene
    gene_subexonuniq = defaultdict(int) # keyed by gene

    for i, read_exons in enumerate(seb):
        for subexon_id in read_exons:
            gene_id = subexon_id[0]
            gene_subexonbags[gene_id][i].append(subexon_id[1:])
            gene_subexonuniq[gene_id] += (subexon_id in se_uniq)

    uniqgene_ids = [g_id for g_id in gene_subexonuniq if gene_subexonuniq[g_id]]
    if len(uniqgene_ids) == 1:
        # uniquely assigned
        gene_id = uniqgene_ids.pop()
        seb = gene_subexonbags[gene_id]
        path = tuple(sorted([tuple(sorted(seb[0])), 
                            tuple(sorted(seb[1]))]))
        sep = (gene_id, path)
    else:
        sep = None
    return sep

def scanBAM(sf, ga, func, progress, qc, split_mode):
    split = (None, None, None, None)
    reads = {}
    cargo = defaultdict(lambda: defaultdict(int))
    i_pair = 0
    i_pass = 0
    for i_any, i_read in enumerate(sf):
        # LOG
        if (i_any + 1) % progress == 0:
            sys.stderr.write("reads: %d (processed) %d (pass QC) %d (paired) %d (unpaired)\n" %
                             (i_any + 1, i_pass*2, i_pair*2, len(reads)))
        # QC
        if qc == "strict":
            if not i_read.is_proper_pair:
                continue
        elif qc == "loose":
            if i_read.is_unmapped or i_read.mate_is_unmapped:
                continue

        if i_read.is_secondary:
            continue

        i_read_chr = sf.getrname(i_read.tid)
        j_read_chr = sf.getrname(i_read.mrnm)

        if i_read_chr != j_read_chr:
            continue

        i_read_strand = "-" if i_read.is_reverse else "+"
        j_read_strand = "-" if i_read.mate_is_reverse else "+"

        if i_read_strand == j_read_strand:
            continue

        i_read_id = (i_read_chr, i_read.pos,  i_read_strand, i_read.qname, i_read.is_read2)
        j_read_id = (j_read_chr, i_read.mpos, j_read_strand, i_read.qname, i_read.is_read1)

        j_read = reads.pop(j_read_id, None)
        if j_read:

            i_pair+=1

            if i_read_strand == "+":
                if (i_read.positions[0]  > j_read.positions[0]) or \
                   (i_read.positions[-1] > j_read.positions[-1]):
                    continue
                frag_start = i_read.positions[0]
                frag_end = j_read.positions[-1]
            else:
                if (i_read.positions[0]  < j_read.positions[0]) or \
                   (i_read.positions[-1] < j_read.positions[-1]):
                    continue
                frag_start = j_read.positions[0]
                frag_end = i_read.positions[-1]
                
            i_pass+=1

            # For stranded libraries:
            #  - read2 is from the "+" strand
            frag_chr, frag_strand = (i_read_chr, i_read_strand) if i_read.is_read2 else \
                                    (j_read_chr, j_read_strand)
            seb = subexonbag(i_read, frag_chr, frag_strand, j_read, frag_chr, frag_strand, ga)
            # i_read was second in alignment, j_read was first in alignment
            # there is no information associated with this order.
            if split_mode == "fragment":
                frag_start = min(i_read.positions[0], j_read.positions[0])
                frag_end = max(i_read.positions[-1], j_read.positions[-1])
                split = (frag_chr, frag_start, frag_end, frag_strand)
            elif split_mode == "qname":
                split = i_read.qname
            
            ##
            func(cargo, seb, split)
        else:
            # first read in alignment
            reads[i_read_id] = i_read
    return dict(cargo)

def seb2sep(seb_cargos, se_uniq):

    sep_cargos = defaultdict(dict) # uniquely assigned

    for seb, cargo in seb_cargos.iteritems():
        sep = subexonpath(seb, se_uniq)
        if sep:
            ## split subexons by genes count unique ones
            sep_cargos[sep[0]][sep[1]] = cargo

    return dict(sep_cargos)
