#!/usr/bin/env python

import pandas as pd
from pybedtools import BedTool, helpers
import sys
import os
from multiprocessing import Pool

""" Adapted from Gabrielle's script """
""" Compute p-values for overlapping target interactions for 2 entities: snoRNA-snoRNA, RBP-RBP, snoRNA-RBP """

def shuffle_val(args):
    bed_int, genome_file, pc_bed, i, outfile_prefix = args
    # shuffle interaction bed
    shuffled = bed_int.shuffle(incl=pc_bed.fn, seed=i, g=genome_file)
    shuffled.moveto(outfile_prefix + str(i) + '.bed')


def calc_overlap(args):
    shuffled_prefix, real_nb_overlaps, bed_target, i = args
    shuffled = BedTool(shuffled_prefix + str(i) + '.bed')
    # intersect shuffled int_bed target
    bed_out = bed_target.intersect(shuffled, s=True, u=True)
    # count nb of shuffled int
    shuffle_nb_overlaps = bed_out.count()
    if shuffle_nb_overlaps >= real_nb_overlaps:
        n = 1
    else:
        n = 0
    out_fn = bed_out.fn
    os.remove(out_fn)
    return n

def sno_ovlp(int_file): # format snoGloBe predictions for snoRNA overlap calculations
    df = pd.read_csv(int_file,
                     sep='\t',
                     names=['seqname', 'target_window_start', 'target_window_end', 'snoid', 'count', 'strand']) 
    return df

def rbp_rbp_ovlp(int_file): # format ENCODE RBP binding interactions for RBP overlap calculations
    df = pd.read_csv(int_file,
                     sep='\t',
                     names=['seqname', 'target_window_start', 'target_window_end',
                            'rep', 'score', 'strand']) 
    df.rename(columns={"rep": "rbp_name"},inplace=True)
    rbp = os.path.basename(int_file).split("_")
    df['rbp_name'] = rbp[0]
    return df

def get_name(f_name): # extract snoRNA or RBP name from file name
    name = os.path.basename(f_name).split("_")
    if name[0] == "NR": # snoRNAs with id staring with 'NR'
        name = '_'.join([name[0],name[1]])
    else: # snoRNAs with id starting with 'ENSG00000' or RBP
        name = name[0]
    return name

def main():
    int_file = sys.argv[1] # source file (snoGloBe predicted interactions for sno-RBP overlap calculation)
    target_path = sys.argv[2] # target file (path to all rbp merged bed for sno-RBP overlap calculation)
    gtf_file = sys.argv[3] # full human gtf file in tsv format
    genome_file = sys.argv[4] # tsv file with name and length of every chromosome
    outfile = sys.argv[5] # output file name with complete path (one file per source)
    nb_threads = int(sys.argv[6]) # nb of threads to use in parallel
    temp_dir = sys.argv[7] # $SLURM_TMPDIR
    ovlp_type = sys.argv[8] # choose between: sno, rbp, sno_rbp

    helpers.set_tempdir(temp_dir)

    # create pc_region +- 12 nts bed
    df_gtf = pd.read_csv(gtf_file, sep='\t', usecols=['seqname', 'start', 'end', 'strand', 'feature', 'gene_biotype'], # gtf
                         dtype={'seqname': str})
    df_gtf = df_gtf[(df_gtf.feature == 'gene') & (df_gtf.gene_biotype == 'protein_coding')]
    df_gtf['start'] -= 12
    df_gtf['end'] += 12
    df_gtf['name'] = 'pc'
    df_gtf['score'] = 0
    df_gtf[['start', 'end']] = df_gtf[['start', 'end']].astype(int)
    df_gtf = df_gtf.sort_values(['seqname', 'start', 'end'])
    pc_bed = BedTool.from_dataframe(df_gtf[['seqname', 'start', 'end', 'name', 'score', 'strand']])
    del df_gtf

    target_list = os.listdir(target_path)

    # intersect source & target interactions depending on overlap type
    df = None
    if ovlp_type == "sno_rbp" or ovlp_type == "sno":
        df = sno_ovlp(int_file)
    elif ovlp_type == "rbp":
        df = rbp_rbp_ovlp(int_file)
    else:
        sys.exit('Option not supported by this script. Please choose between: sno_rbp, sno, rbp.')

    bed_int = BedTool.from_dataframe(df) # source interactions in bed format

    max_iteration = 100000

    shuffle_prefix = os.path.join(temp_dir, 'shuffled_%s' % (os.path.basename(int_file)))

    # create all the shuffled interaction bedfiles
    with Pool(nb_threads) as p:
        res_shuffle = p.map(shuffle_val, [[bed_int, genome_file, pc_bed, i, shuffle_prefix]
                                  for i in range(max_iteration)])

    # avoid re-running the source-target pair if was already computed (example if job timed out)
    already_done = ''
    try:
        with open(outfile) as f:
            for line in f:
                already_done += (line.strip() + ',')
    except FileNotFoundError:
        pass

    for target_file in target_list:
        if os.path.basename(int_file) == os.path.basename(target_file): # don't run self pairs for sno-sno and rbp-rbp overlap calculations
            continue
        if os.path.basename(target_file) in already_done:
            continue
        bed_target = None
        if ovlp_type == "sno":
            bed_target = BedTool.from_dataframe(sno_ovlp(os.path.join(target_path, target_file)))
        else: # rbp file as target
            bed_target = BedTool(os.path.join(target_path, target_file))
        # intersect bed int target
        bed_out = bed_target.intersect(bed_int, s=True, u=True)
        # count nb of intersections
        real_nb_overlaps = bed_out.count()

        prev_nb_shuffle = 0
        nb_shuffle = 10
        n = 0


        # to limit the computation time, do as few tests as possible. Start by 10, if suffled >= real (n) at most once,
        # then continue to 100 tests, and so on until max_iteration is reached
        while n <= 1 and nb_shuffle <= max_iteration :
            with Pool(nb_threads) as p:
                res = p.map(calc_overlap, [[shuffle_prefix, real_nb_overlaps, bed_target, i]
                                          for i in range(prev_nb_shuffle, nb_shuffle)])
            n += sum(res)
            prev_nb_shuffle = nb_shuffle
            nb_shuffle *= 10

        # compute pval
        # pval = n / nb_shuffle
        curr_source = get_name(int_file)
        curr_target = get_name(target_file)
        with open(outfile, 'a+') as w:
            w.write('\t'.join([curr_source, curr_target, str(n), str(prev_nb_shuffle)]) + '\n')
        # the resulting file will have the followin columns:
        # source filename, target filename, nb of times shuffled >= real, number of tests done
        # the p-val = nb of times shuffled >= real / number of tests done

if __name__ == '__main__':
    main()
