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
    shuffled_prefix, real_nb_overlaps, bed_rbp, i = args
    shuffled = BedTool(shuffled_prefix + str(i) + '.bed')
    # intersect shuffled int_bed rbp
    bed_out = bed_rbp.intersect(shuffled, s=True, u=True)
    # count nb of shuffled int
    shuffle_nb_overlaps = bed_out.count()
    if shuffle_nb_overlaps >= real_nb_overlaps:
        n = 1
    else:
        n = 0
    out_fn = bed_out.fn
    os.remove(out_fn)
    return n


def main():
    int_file = sys.argv[1] # file with predicted interactions
    rbp_path = sys.argv[2] # path to all rbp merged bed
    gtf_file = sys.argv[3] # full human gtf file in csv format
    genome_file = sys.argv[4] # tsv file with name and length of every chromosome
    outfile = sys.argv[5] # output file name with complete path
    nb_threads = int(sys.argv[6]) # nb of threads to use in parallel
    temp_dir = sys.argv[7] # $SLURM_TMPDIR

    helpers.set_tempdir(temp_dir)

    # create pc_region +- 12 nts bed
    df_gtf = pd.read_csv(gtf_file, sep='\t', usecols=['seqname', 'start', 'end', 'strand', 'feature', 'gene_biotype'],
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

    rbp_list = os.listdir(rbp_path)

    # intersect bed int rbp
    df = pd.read_csv(int_file,
                     sep='\t',
                     names=['seqname', 'target_window_start', 'target_window_end',
                            'snoid', 'count', 'strand', 'sno_window_start',
                            'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
    df = df[df.snoid != 'snoid']  # drop header if present
    df = df.drop(['target'], axis=1)
    df = df.drop_duplicates()
    df = df[['seqname', 'target_window_start', 'target_window_end', 'snoid', 'count', 'strand']]
    bed_int = BedTool.from_dataframe(df)

    max_iteration = 100000

    shuffle_prefix = os.path.join(temp_dir, 'shuffled_%s' % (os.path.basename(int_file)))

    # create all the shuffled interaction bedfiles
    with Pool(nb_threads) as p:
        res_shuffle = p.map(shuffle_val, [[bed_int, genome_file, pc_bed, i, shuffle_prefix]
                                  for i in range(max_iteration)])

    # avoid re-running the sno-rbp pair if was already computed (example if job timed out)
    already_done = ''
    try:
        with open(outfile) as f:
            for line in f:
                already_done += (line.strip() + ',')
    except FileNotFoundError:
        pass

    for rbp_file in rbp_list:
        if os.path.basename(rbp_file) in already_done:
            continue
        bed_rbp = BedTool(os.path.join(rbp_path, rbp_file))
        # intersect bed int rbp
        bed_out = bed_rbp.intersect(bed_int, s=True, u=True)
        # count nb of intersections
        real_nb_overlaps = bed_out.count()

        prev_nb_shuffle = 0
        nb_shuffle = 10
        n = 0


        # to limit the computation time, do as few tests as possible. Start by 10, if suffled >= real (n) at most once,
        # then continue to 100 tests, and so on until max_iteration is reached
        while n <= 1 and nb_shuffle <= max_iteration :
            with Pool(nb_threads) as p:
                res = p.map(calc_overlap, [[shuffle_prefix, real_nb_overlaps, bed_rbp, i]
                                          for i in range(prev_nb_shuffle, nb_shuffle)])
            n += sum(res)
            prev_nb_shuffle = nb_shuffle
            nb_shuffle *= 10

        # compute pval
        #pval = n / nb_shuffle
        with open(outfile, 'a+') as w:
            w.write('\t'.join([os.path.basename(int_file), os.path.basename(rbp_file), str(n), str(prev_nb_shuffle)]) + '\n')
        # the resulting file will have the followin columns:
        # sno filename, rbp filename, nb of times shuffled >= real, number of tests done
        # the p-val = nb of times shuffled >= real / number of tests done

if __name__ == '__main__':
    main()
