#!/usr/bin/env python2
from collections import defaultdict
from itertools import chain
import numpy as np
import optparse
import re
import os
import sys
import json

    

# multinomialCI <- function(x,alpha,verbose=FALSE){
  
#   n = sum(x, na.rm=TRUE);
#   k = length(x);
#   p = x/n;
#   c = 0;
#   pold=0;
#   for(cc in 1:n){
#    p = .truncpoi(cc,x,n,k);
#    if(p > 1-alpha && pold < 1-alpha) { c = cc; break; };
#    pold=p;
#   }

#   salida = matrix(0,k,2);
#   delta=(1-alpha-pold)/(p-pold);
#   out=matrix(0,k,5);
#   num=matrix(0,k,1);
#   c=c-1;  
#   vol1=1;
#   vol2=1;
#   for(i in 1:k){
#    num[i,1]=i;
#    obsp=x[i]/n;
#    out[i,1]=obsp;
#    out[i,2]=obsp-c/n;
#    out[i,3]=obsp+c/n+2*delta/n;
#    if(out[i,2]<0){ out[i,2]=0; }
#    if(out[i,3]>1){ out[i,3]=1; }
#    out[i,4]=obsp-c/n-1/n;
#    out[i,5]=obsp+c/n+1/n;
#    if(out[i,2]<0){ out[i,2]=0; }
#    if(out[i,3]>1){ out[i,3]=1; }
#    vol1=vol1*(out[i,3]-out[i,2]);
#    vol2=vol2*(out[i,5]-out[i,4]);
   
#    salida[i,1] = out[i,2];
#    salida[i,2] = out[i,3];
#   }
#   c1=c('PROPORTION', 'LOWER(SG)','UPPER(SG)','LOWER(C+1)','UPPER(C+1)');
#   cov=100*(1-alpha);
#   sg=(x+delta)/n;
#   c2=c('SG-midpoint');
#   if(verbose==TRUE){
#     print('-------------------------------------------------------------');
#     print(paste('    ',cov,'% SIMULTANEOUS CONFIDENCE INTERVALS'));
#     print('       BASED ON THE METHODS OF SISON AND GLAZ');
#     print('-------------------------------------------------------------');
#     print(paste('C = ',c));
#     print(paste('P(c+1) = ',p));
#     print(paste('P(c)   = ',pold));
#     print(paste('delta =  ',delta));
#     print(paste('Volume(SG) = ',vol1));
#     print(paste('Volume(C+1)= ',vol2));
#     print(paste(c1));
#     print(out);
#     #print(paste(c(num,out)));
#     print(paste(c2));
#     #print(paste(c(num,sg)));
#     print(sg);
#   }
#   return(salida);
#}


def get_base(unit='bit'):
    if unit == 'bit':
        log = np.log2
    elif unit == 'nat':
        log = np.log 
    elif unit in ('digit', 'dit'):
        log = np.log10  
    else:
        raise ValueError('The "unit" "%s" not understood' % unit)
    return log

def shannon_entropy(freq, unit='bit'):
    """Calculates the Shannon Entropy (H) of a frequency.
    
    Arguments:
    
        - freq (``numpy.ndarray``) A ``Freq`` instance or ``numpy.ndarray`` with 
          frequency vectors along the last axis.
        - unit (``str``) The unit of the returned entropy one of 'bit', 'digit' 
          or 'nat'.
    """
    log = get_base(unit)
    shape = freq.shape # keep shape to return in right shape
    Hs = np.ndarray(freq.size / shape[-1]) # place to keep entropies
    # this returns an array of vectors or just a vector of frequencies
    freq = freq.reshape((-1, shape[-1])) 
    # this makes sure we have an array of vectors of frequencies
    freq = np.atleast_2d(freq)
    # get fancy indexing
    positives = freq != 0.
    for i, (freq, idx) in enumerate(zip(freq, positives)):
        freq = freq[idx] # keep only non-zero
        logs = log(freq) # logarithms of non-zero frequencies
        Hs[i] = -np.sum(freq * logs)
    Hs.reshape(shape[:-1])
    return Hs
 

def jensen_shannon_divergence(freq, weights =None, unit='bit'):
    """
    Calculates the Jensen-Shannon Divergence (Djs) of two or more frequencies.
    The weights are for the relative contribution of each frequency vector. 
    
    Arguments:
    
        - freq (``numpy.ndarray``) rank-2 array of frequencies along the last dimension
        - weights (``numpy.ndarray``) A rank-1 array with a weight for each  frequency vector.
        - unit (``str``) see: the function ``shannon_entropy``.
    """
    if weights is not None:
        if len(freq) != len(weights):
            raise ValueError('The number of frequencies and weights do not match.')
    if freq.ndim != 2:
        raise ValueError('At least two frequencies in a rank-2 array expected.')
    weighted_average = np.average(freq, axis=0, weights=weights)
    H_avg_freq = shannon_entropy(weighted_average, unit)
    H_freq = shannon_entropy(freq, unit)
    avg_H_freq = np.average(H_freq, weights=weights)
    JSD = H_avg_freq - avg_H_freq
    return JSD

def parse_sep(fh):
    seps_counts = defaultdict(lambda: defaultdict(dict))
    for line in fh:
        gene_id, seps0, seps1, count = line.rstrip().split("\t")
        sep = [[], []]
        for i, sep_str in enumerate((seps0, seps1)):
            if sep_str:
                for se_str in sep_str.split("-"):
                    e, part = map(int, se_str.split(":"))
                    sep[i].append((e, part))
            sep[i] = tuple(sep[i])
        seps_counts[gene_id][tuple(sep)] = count
    return {"assigned":seps_counts}

def load_samples(fns, group_ids, format):
    sample_ids = [fn.rsplit("/", 1)[-1].rsplit(".", 1)[0] for fn in fns]
    sep_ids = zip(group_ids, sample_ids)
    sample_seps = {}
    for sep_id, fn in zip(sep_ids, fns):
        fn = os.path.abspath(os.path.expanduser(fn))
        if format == "json":
            sep = json.load(open(fn))
        else:
            sep = parse_sep(open(fn))
        sample_seps[sep_id] = sep["assigned"]
    return sample_seps

def seps_to_arrs(sample_seps):
    sample_ids = tuple(sample_seps.keys())
    # union of gene_ids across samples
    gene_ids = tuple(sorted(set(chain.from_iterable(sample_sep.keys() for sample_sep in sample_seps.values()))))
    sep_arrs = {}
    for gene_id in gene_ids:
        sep_ids = tuple(sorted(set(chain.from_iterable(sample_sep[gene_id] for sample_sep in sample_seps.values()))))
        sep_arr = np.zeros((len(sample_ids), len(sep_ids)), dtype=long)
        for row_id, sample_id in enumerate(sample_ids):
            for col_id, sep_id in enumerate(sep_ids):
                # make sure not to modify the defaultdicts
                if (gene_id in sample_seps[sample_id]) and (sep_id in sample_seps[sample_id][gene_id]):
                    sep_arr[row_id, col_id] = sample_seps[sample_id][gene_id][sep_id]
        sep_arrs[gene_id] = (sep_ids, sep_arr)
    return (sample_ids, sep_arrs)

def arr_valid(min_count):
    pass_count = arr >= min_count

    #

# freq = (arr.transpose() / np.array(arr.sum(axis=1), dtype=float)).transpose()

# pass_freq = (arr.transpose() / np.array(arr.sum(axis=1), dtype=float) > 0).transpose()
#     #
#     arr_pass = pass_count & pass_freq
#     row_pass = arr_pass.any(axis=1)






# sep_counts = arr.sum(axis=0)
# sample_counts = arr.sum(axis=1)



    



# # test
# fns = ["test/diff/test1.sep", "test/diff/test2.sep", "test/diff/test3.sep"]
# group_ids = ["A", "B", "C"]
# format = "sep"

# sample_seps = load_samples(fns, group_ids, format)
# sample_arrs = seps_to_arrs(sample_seps)


# # a = seps_to_arrs(seps)
    

# #     print gene_id, len(sep_ids)
                 


# # [0].keys()




# # a["assigned"].keys()


# #     for (row_idx, sep_id) in enumerate(sep_ids):
# #         for (col_idx, col_id) in enumerate(
# #         sep_arr[]




# # seps1 = a["assigned"]["ENSG00000246228.2"]
# # seps2 = a["assigned"]["ENSG00000246228.2"]
# # seps3 = a["assigned"]["ENSG00000246228.2"]
# # gene_seps = [seps1, seps2, seps3]



    


# a = parse_sep(open("t"))



if __name__ == "__main__":
    optParser = optparse.OptionParser(

        usage = "%prog [options] sep_file(s)",
        
        description = \
        "This script takes sepath-count output and identifies " +
        "differentially spliced genes.",
        
        epilog = \
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) " +
        "Michigan Center for Translational Pathology (c) 2014 "
    )

    optParser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                         help="run-time messages printed to stderr")

    optParser.add_option("-g", "--groups", type="string", dest="groups", 
                         help="")

    optParser.add_option("-f", "--format", type="string", dest="format", 
                         default="sep", help="input file format 'sep' or 'json'")

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()
    
    if len(args) < 2:
        optParser.print_help()
        sys.exit(1)

    groups = opts.groups.split(",") if opts.groups else map(str, range(1, 1+len(args))) 

    if len(groups) != len(args):
        sys.stderr.write("error: The length of 'groups' does not match the " +
                         "number of 'sep' files.\n")
        sys.exit(1)
    
    seps = load_seps(args, groups, opts.format)


    
