import argparse
import functools
import gzip
import itertools
import os
import sys
import tempfile
import tomllib
from collections import Counter, defaultdict
from glob import glob
from typing import Callable, Union, Any, Optional

import scipy
from multiprocess import Pool  # pip install multiprocess
import rapidfuzz  # pip install rapidfuzz
import numpy as np  # pip install numpy
import pandas as pd  # pip install pandas
from Bio.Seq import reverse_complement  # pip install biopython
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from scipy.stats import combine_pvalues, hmean, entropy
from statsmodels.stats.multitest import multipletests  # pip install statsmodels
from tqdm.auto import tqdm  # pip install tqdm
import matplotlib.pyplot as plt  # pip install matplotlib
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import seaborn as sns  # pip install seaborn


def infer_position_and_variant_locations(filename, capture_seq, rev_capture_seq, barcode_offset):
    with gzip.open(filename, mode='rt') as f:
        position_locations = []
        variant_locations = []
        count = 0
        missing_count = 0
        max_count = 1_000
        for (title, sequence, quality) in FastqGeneralIterator(f):
            if count + 1 > max_count:
                break
            sequence = sequence.upper()
            if capture_seq in sequence:  # Need to get the reverse compliment
                capture_location = sequence.find(capture_seq)
                position_location = (capture_location - barcode_offset, capture_location, True)
                variant_location = (capture_location - (barcode_offset + barcode_offset + 4),
                                    capture_location - (barcode_offset + 4), True)
            elif rev_capture_seq in sequence:
                sequence = reverse_complement(sequence)
                capture_location = sequence.find(capture_seq)
                position_location = (capture_location - barcode_offset, capture_location, False)
                variant_location = (capture_location - (barcode_offset + barcode_offset + 4),
                                    capture_location - (barcode_offset + 4), False)
            else:  # The capture seq still isnt in the sequence
                # print("WARNING: MISSING CAPTURE SEQ IN READ " + title + " Skipping")
                # print(sequence)
                missing_count += 1
                continue
            position_locations.append(position_location)
            variant_locations.append(variant_location)
            count += 1

        position_location, position_count = Counter(position_locations).most_common(1)[0]
        variant_location, variant_count = Counter(variant_locations).most_common(1)[0]

        print(f"Most common position location: {position_location} ({position_count}/{max_count} valid reads)")
        print(f"Most common variant location: {variant_location} ({variant_count}/{max_count} valid reads)")
        print(f"Reads missing a capture sequence: {missing_count}")

    return position_location, variant_location


def map_variants(barcodes_fasta_path: str, twist_barcodes_fasta_path: str, twist_position_only_barcodes_fasta_path: str,
                 fastqs: list[str], capture_seq, barcode_offset, n_processes):
    # Generate a temporary FASTQ file with just the barcode region
    # Note that we are concatenating the files here
    rev_capture_seq = reverse_complement(capture_seq)

    def read_fastq():
        for f in fastqs:
            # position_location, variant_location = infer_position_and_variant_locations(f, capture_seq,
            #                                                                            rev_capture_seq,
            #                                                                            barcode_offset)
            #
            # position_start, position_end, is_forward = position_location
            # variant_start, variant_end, is_forward = variant_location

            for (title, sequence, quality) in FastqGeneralIterator(gzip.open(f, mode='rt')):
                capture_seq_loc = sequence.find(capture_seq, (
                            barcode_offset * 2) + 4)  # Skip the first 36 (2 barcodes + spacer) bases since we need at least 36 bases to extract a barcode...
                if capture_seq_loc == -1:
                    capture_seq_loc = sequence.find(rev_capture_seq, (barcode_offset * 2) + 4)
                    if capture_seq_loc == -1:
                        continue
                    sequence = reverse_complement(sequence)
                position_end = capture_seq_loc
                # position_start = position_end - barcode_offset
                variant_start = position_end - (barcode_offset + barcode_offset + 4)
                # variant_end = position_end - barcode_offset
                yield title, sequence[variant_start:position_end]
        # yield from itertools.chain(*[FastqGeneralIterator(gzip.open(f, mode='rt')) for f in fastqs])

    def map_barcode_seqs(barcodes_file: str, position_only: bool, to_ignore=None):
        counter = defaultdict(int)
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_fasta = os.path.join(temp_dir, 'temp.fasta')
            for (title, bc) in tqdm(read_fastq(),
                                    desc=f'Extracting {"Position" if position_only else "Clonal"} barcodes',
                                    unit='reads', leave=True):
                if to_ignore is not None and bc in to_ignore:
                    continue
                counter[bc] += 1
            with open(temp_fasta, 'w') as temp_f:
                for bc in counter.keys():
                    temp_f.write(f">{bc}\n")
                    temp_f.write(f"{bc[-16:] if position_only else bc}\n")

            print("Performing mmseqs2 alignment")
            # Perform mmseqs2 alignment using easy-search
            cmd = (f"mmseqs easy-search "
                   f"{temp_fasta} "  # Query 
                   f"{barcodes_file} "  # Database
                   f"{temp_dir}/results.m8 "  # Output file
                   f"{temp_dir}/tmp "  # Temp directory
                   f"--comp-bias-corr 0 "  # Disable amino acid composition bias correction
                   f"--mask 0 "  # Disable low-complexity masking
                   f"--mask-profile 0 "  # Disable profile-based masking
                   f"--min-aln-len {'13' if position_only else '30'} "  # Must mostly cover at least a barcode
                   f"-c 0.8 "  # Alignment must cover 80% of the query
                   f"--strand 1 "  # Disable reverse frames
                   f"-k 6 "  # Seed size
                   # f"--greedy-best-hits 1 "  # Greedy best hits 
                   f"--threads {n_processes} "  # Number of threads
                   f"--format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq' "  # Output format
                   f"--search-type 3")  # Nucleotide search
            print(f"> {cmd}")
            os.system(cmd)

            # Parse the results (From Melissa)
            mapping = pd.read_csv(f"{temp_dir}/results.m8", header=None, sep='\t')
            mapping.columns = ['query_bc',
                               'target_barcode',
                               'percent_identity',
                               'alignment_length',
                               'mismatches',
                               'gapopen',  # Number of gap openings
                               'qapstart',
                               'qapend',
                               'tstart',
                               'tend',
                               'evalue',  # Expect value
                               'bitscore',  # Bit score
                               'corrected_sequence']
            print("Mapped", mapping.shape[0], "barcodes")
            return mapping, counter

    print("Mapping positional barcodes for spike ins")
    positional_mappings, positional_counter = map_barcode_seqs(twist_barcodes_fasta_path, True)
    positional_mappings['umi'] = positional_mappings['query_bc'].str[
        :16]  # We only map with positional barcodes so use the full sequence to get the umi info
    positional_mappings['variant'] = positional_mappings['target_barcode'].str.split(':').str[0]
    positional_mappings = positional_mappings.sort_values(by=['query_bc', 'evalue', 'mismatches', 'gapopen'],
                                                          ascending=[True, True, True, True])
    positional_mappings = positional_mappings.loc[~positional_mappings['query_bc'].duplicated()]
    print(f"Selected {positional_mappings.shape[0]} unique control positional barcode matches")
    all_reads = sum(positional_counter.values())
    unmapped_reads = positional_counter.keys() - set(positional_mappings['query_bc'].values)
    unmapped_reads_count = sum(positional_counter[bc] for bc in unmapped_reads)
    print(
        f"{unmapped_reads_count} barcoded reads did not map to a known barcode ({unmapped_reads_count / all_reads * 100:.2f}%)")
    positional_mappings['corrected'] = (positional_mappings['mismatches'] > 0) | (positional_mappings['gapopen'] > 0)
    print(
        f"Corrected {positional_mappings['corrected'].sum()} unique barcodes ({positional_mappings['corrected'].sum() / positional_mappings.shape[0] * 100:.2f}%)")
    print(
        f"This represents {sum(positional_counter[bc] for bc in positional_mappings[positional_mappings['corrected']]['query_bc']) / all_reads * 100:.2f}% of reads")

    print(f"Excluding {positional_mappings.shape[0]} barcodes ({all_reads} reads) from clonal mapping")
    print("Mapping clonal barcodes")
    clonal_mappings, clonal_counter = map_barcode_seqs(barcodes_fasta_path, False,
                                                       set(positional_mappings.query_bc.values))
    clonal_mappings['umi'] = clonal_mappings['corrected_sequence'].str[:16]
    clonal_mappings['variant'] = clonal_mappings['target_barcode'].str.split(':').str[0]
    clonal_mappings['map_score'] = clonal_mappings['target_barcode'].str.split(':').str[1].astype(float)
    clonal_mappings = clonal_mappings.sort_values(by=['query_bc', 'evalue', 'mismatches', 'gapopen', 'map_score'],
                                                  ascending=[True, True, True, True, False])
    clonal_mappings = clonal_mappings.loc[~clonal_mappings['query_bc'].duplicated()]
    print(f"Selected {clonal_mappings.shape[0]} best unique barcode matches")
    all_reads = sum(clonal_counter.values())
    unmapped_reads = clonal_counter.keys() - set(clonal_mappings['query_bc'].values)
    unmapped_reads_count = sum(clonal_counter[bc] for bc in unmapped_reads)
    print(
        f"{unmapped_reads_count} barcoded reads did not map to a known barcode ({unmapped_reads_count / all_reads * 100:.2f}%)")
    clonal_mappings['corrected'] = (clonal_mappings['mismatches'] > 0) | (clonal_mappings['gapopen'] > 0)
    print(
        f"Corrected {clonal_mappings['corrected'].sum()} unique barcodes ({clonal_mappings['corrected'].sum() / clonal_mappings.shape[0] * 100:.2f}%)")
    print(
        f"This represents {sum(clonal_counter[bc] for bc in clonal_mappings[clonal_mappings['corrected']]['query_bc']) / all_reads * 100:.2f}% of reads")

    # Finally, any remaining unmapped barcodes are checked again to map to X using TWIST position barcodes
    print("Mapping unclassified clonal barcodes to TWIST positional barcodes")
    position_clonal_mappings, position_clonal_counter = map_barcode_seqs(twist_position_only_barcodes_fasta_path, True,
                                                                         set(positional_mappings.query_bc.values) | set(
                                                                             clonal_mappings.query_bc.values))
    position_clonal_mappings['umi'] = position_clonal_mappings['corrected_sequence'].str[:16]
    position_clonal_mappings['variant'] = position_clonal_mappings['target_barcode'].str.split(':').str[0]
    position_clonal_mappings = position_clonal_mappings.sort_values(by=['query_bc', 'evalue', 'mismatches', 'gapopen'],
                                                                    ascending=[True, True, True, True])
    position_clonal_mappings = position_clonal_mappings.loc[~position_clonal_mappings['query_bc'].duplicated()]
    print(f"Selected {position_clonal_mappings.shape[0]} best unique barcode matches")
    all_reads = sum(position_clonal_counter.values())
    unmapped_reads = position_clonal_counter.keys() - set(position_clonal_mappings['query_bc'].values)
    unmapped_reads_count = sum(position_clonal_counter[bc] for bc in unmapped_reads)
    print(
        f"{unmapped_reads_count} barcoded reads did not map to a known barcode ({unmapped_reads_count / all_reads * 100:.2f}%)")
    position_clonal_mappings['corrected'] = (position_clonal_mappings['mismatches'] > 0) | (
                position_clonal_mappings['gapopen'] > 0)
    print(
        f"Corrected {position_clonal_mappings['corrected'].sum()} unique barcodes ({position_clonal_mappings['corrected'].sum() / position_clonal_mappings.shape[0] * 100:.2f}%)")
    print(
        f"This represents {sum(position_clonal_counter[bc] for bc in position_clonal_mappings[position_clonal_mappings['corrected']]['query_bc']) / all_reads * 100:.2f}% of reads")

    # Map barcode -> variant
    print("Collecting counts...")
    bc2variant = dict()
    for i, row in itertools.chain(clonal_mappings.iterrows(), positional_mappings.iterrows(),
                                  position_clonal_mappings.iterrows()):
        bc2variant[row['query_bc']] = (row['umi'], row['variant'])
    # Create a dataframe
    rows = []
    for bc, (umi, variant) in bc2variant.items():
        if bc in positional_counter:
            rows.append((umi, variant, positional_counter[bc]))
        else:
            rows.append((umi, variant, clonal_counter[bc]))

    mapping = pd.DataFrame(rows, columns=['umi', 'variant', 'count'])

    print(f"Total counts: {mapping['count'].sum()}")

    print("Processing barcodes")
    # Duplicate the rows for each variant if there are multiple
    mapping['variant'] = mapping['variant'].str.split(',')
    mapping = mapping.explode('variant')
    # Drop low-quality STOP variants
    is_longread_stop = mapping['variant'].str.endswith('_')
    # mapping.loc[is_stop, 'variant'] = 'STOP'
    mapping = mapping[~is_longread_stop]
    # Remove any duplicated entries in case there are multiple mappings to STOP
    mapping = mapping.drop_duplicates()
    # Group by variant and UMI
    mapping = mapping.groupby(['variant', 'umi']).agg({'count': 'sum'}).reset_index()

    return mapping


def make_counts(experiment_names, file_names, barcodes_file, twist_file, spike_ins=[],
                capture_seq="CGCTTTAAGGCCGGTCCTAGCAA", barcode_offset=16, n_processes=3,
                special_vars=['STOP', 'mSTING_WT', 'RRDD', 'RRAA', 'dCTT', 'dUPR']):
    print("Compiling barcodes")
    longread_bcs = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
    df = pd.read_table(barcodes_file)
    df = df[~df['aa_mut'].isin(special_vars)]
    with open(longread_bcs, 'w') as f:
        for i, row in df.iterrows():
            f.write(f">{row['aa_mut']}:{row['score']}:{row['bc']}\n{row['bc']}\n")

    twist_bcs = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
    df = pd.read_csv(twist_file)
    df = df[df.residue.isin(spike_ins)]
    with open(twist_bcs, 'w') as f:
        for i, row in df.iterrows():
            f.write(f">{row['residue']}:{row['barcode']}\n{row['barcode']}\n")

    non_spiked_twist_bcs = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
    df = pd.read_csv(twist_file)
    df = df[~df.residue.isin(spike_ins)]
    with open(non_spiked_twist_bcs, 'w') as f:
        for i, row in df.iterrows():
            f.write(f">{row['residue']}:{row['barcode']}\n{row['barcode']}\n")

    frames = []
    for files, name in zip(file_names, experiment_names):
        mapping = map_variants(longread_bcs, twist_bcs, non_spiked_twist_bcs, files, capture_seq, barcode_offset,
                               n_processes)
        mapping['experiment'] = name
        frames.append(mapping)

    df = pd.concat(frames)
    return df


# def make_counts(experiment_names, file_names, barcode_df, capture_seq="CGCTTTAAGGCCGGTCCTAGCAA", barcode_offset=36, n_processes=3):
#     bc2variants = defaultdict(list)
#     barcode_df = barcode_df.reset_index(drop=True)
#     for i, row in barcode_df.iterrows():
#         bc = row['bc']
#         bc2variants[bc].append(i)
#     # Convert keys to tuple
#     possible_bcs = tuple(bc2variants.keys())
#
#     rev_capture_seq = reverse_complement(capture_seq)
#     variant_location = infer_position_and_variant_locations(file_names[0][0], capture_seq, rev_capture_seq, barcode_offset)
#     variant_start, variant_end, is_forward = variant_location
#     frames = []
#
#     with Pool(n_processes) as pool:
#         results = pool.imap_unordered(count_job, [
#             (bc2variants, possible_bcs, files, name, is_forward, variant_start, variant_end)
#             for files, name in zip(file_names, experiment_names)
#         ])
#
#         for name, variants_counted, valid, corrected in tqdm(results, total=len(file_names), desc="Processing Reads", unit="experiments"):
#             print(f"Processed {name}")
#             print(f"Valid reads: {valid}, Corrected reads: {corrected} ({corrected/valid*100:.2f}%)")
#             # Now that its done, replace the indices with the actual variant names
#             df = pd.DataFrame({
#                 'variant': [barcode_df.loc[i, 'variant'] for i in variants_counted],
#                 'umi': [barcode_df.loc[i, 'var_bc'] for i in variants_counted]
#             }).value_counts().reset_index()
#             df.columns = ["variant",
#                             "umi",
#                           'count']
#             df["experiment"] = name
#
#             # # Drop control variants for experiments that don't have them since there may be contamination by setting the values to 0
#             # if drop_controls:
#             #     df.loc[df["position"].isin(controls), "count"] = 0
#                 #df = df[~df["position"].isin(controls)]
#
#             frames.append(df)
#
#     df = pd.concat(frames)
#     return df


def _compute_fitness(treatment_counts, treatment_WT_counts, initial_counts, initial_WT_counts, dt, normalize):
    """
    Compute fitness, related to: https://www.science.org/doi/full/10.1126/sciadv.abb2236
    :param treatment_counts: The counts of the treatment sample. The first dim must indicate variant.
    :param treatment_WT_counts: The counts of the WT in the treatment sample. The first dim must indicate variant.
    :param initial_counts: The counts of the initial sample. The first dim must indicate variant.
    :param initial_WT_counts: The counts of the WT in the initial sample. The first dim must indicate variant.
    :param dt: The time between the two samples.
    :param normalize: Whether to normalize the counts to the WT counts.
    :return: The fitness of each variant.
    """
    treatment_counts = np.array(treatment_counts)
    treatment_WT_counts = np.array(treatment_WT_counts)
    initial_counts = np.array(initial_counts)
    initial_WT_counts = np.array(initial_WT_counts)
    # Fitness = (ln(treatment_counts / sum(treatment_WT_counts)) - ln(initital_counts / sum(initial_WT_counts)) / dt
    # Note this can create NaNs.
    # If treatment counts are 0 and WT counts are 0, the fitness should be 0
    # If treatment counts are 0 and WT counts are not 0, the fitness should be min non-inf value observed minus 1e-6
    # If treatment counts are not 0 and WT counts are 0, the fitness should be NaN, we did not observe it so we can't quantify it properly
    with np.errstate(divide='ignore', invalid='ignore'):
        if normalize:
            ref = treatment_WT_counts
            ref_init = initial_WT_counts
        else:
            ref = treatment_counts
            ref_init = initial_counts
        if ref.ndim == 1:
            ref = ref[:, np.newaxis]
        if ref_init.ndim == 1:
            ref_init = ref_init[:, np.newaxis]
        ref = np.nan_to_num(ref)
        ref_init = np.nan_to_num(ref_init)
        fitness = (np.log(treatment_counts / np.sum(ref, axis=0)) - np.log(
            initial_counts / np.sum(ref_init, axis=0))) / dt
        # Compute relative abundance and take the log ratio of relative abundances
        if np.sum(np.isfinite(fitness)) == 0:
            return fitness  # Can't adjust
        # If treatment counts are 0 and ref counts are 0, 0
        fitness[(treatment_counts == 0) & (initial_counts == 0)] = 0
        # If treatment counts are 0 and ref counts are not 0, min observed non-inf value - 1e-6
        fitness[(treatment_counts == 0) & (initial_counts != 0)] = np.nanmin(fitness[np.isfinite(fitness)]) - 1e-6
        # If treatment counts are not 0 and ref counts are 0, max observed non-inf value + 1e-6
        fitness[(treatment_counts != 0) & (initial_counts == 0)] = np.nanmax(fitness[np.isfinite(fitness)]) + 1e-6

        return fitness


def bootstrap_initial_distribution(df, n=2500):
    # Bootstrap the initial distribution for relative abundance
    return np.random.multinomial(df['count'].sum(), df['count'] / df['count'].sum(), size=n).T


def compute_fitness(experiment_name, initial_df, group_df, dt, normalize, n=2500, counts_cutoff=8):
    # # Require Intersection of position and variant column pairs
    # initial_df = initial_df[initial_df['position'].isin(group_df['position']) & initial_df['variant'].isin(group_df['variant'])]
    # group_df = group_df[group_df['position'].isin(initial_df['position']) & group_df['variant'].isin(initial_df['variant'])]
    # Fill in missing position/variant pairs in either dataframe with 0
    # Collapse the initial_df if there are multiple experiments

    # Sum the counts for each position/variant pair
    initial_df = initial_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control']).agg(
        {'count': 'sum'}).reset_index()
    group_df = group_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control']).agg(
        {'count': 'sum'}).reset_index()

    initial_df = pd.merge(initial_df, group_df[['position', 'variant', 'residue', 'umi', 'is_control']],
                          on=['position', 'variant', 'residue', 'umi', 'is_control'], how='inner').fillna(0)
    group_df = pd.merge(group_df, initial_df[['position', 'variant', 'residue', 'umi', 'is_control']],
                        on=['position', 'variant', 'residue', 'umi', 'is_control'], how='inner').fillna(0)

    # Sort the dataframes by position and variant
    initial_df = initial_df.sort_values(['position', 'variant', 'residue', 'umi', 'is_control'])
    # Re order by initial_df
    initial_df = initial_df.set_index(['position', 'variant', 'residue', 'umi', 'is_control'])
    group_df = group_df.set_index(['position', 'variant', 'residue', 'umi', 'is_control']).reindex(initial_df.index)
    initial_df = initial_df.reset_index()
    group_df = group_df.reset_index()

    # Require at least 8 counts for each position/variant pair in the initial population
    count_filter = (initial_df['count'].values >= counts_cutoff)
    initial_df = initial_df[count_filter]
    group_df = group_df[count_filter]

    # Compute the IQR filter for WT log-ratios
    # First, add matching initial WT counts to the group_df
    wt_df = group_df[group_df.variant == "REF"]
    if wt_df.shape[0] > 2 and normalize:
        wt_df = pd.merge(wt_df, initial_df[initial_df.position == "REF"], on='umi', suffixes=('', '_initial'))
        # Compute the log-ratio of the WT counts
        wt_df['log_ratio_wt'] = np.log(wt_df['count'] / wt_df['count_initial'])
        # Replace infinities with nan
        wt_df['log_ratio_wt'] = wt_df['log_ratio_wt'].replace([np.inf, -np.inf], np.nan)
        # Compute the IQR of the log-ratio
        q3 = np.nanpercentile(wt_df['log_ratio_wt'], 75)
        q1 = np.nanpercentile(wt_df['log_ratio_wt'], 25)
        iqr_log_ratio_wt = q3 - q1
        # Compute the IQR filter (1.5x IQR) for the log-ratio
        wt_df['iqr_filter'] = (wt_df['log_ratio_wt'] >= q1 - 1.5 * iqr_log_ratio_wt) & (
                    wt_df['log_ratio_wt'] <= q3 + 1.5 * iqr_log_ratio_wt)
        # Add the filter to the group_df
        group_df = pd.merge(group_df, wt_df[['position', 'umi', 'variant', 'iqr_filter']], how='left',
                            on=['position', 'variant', 'umi'])
        group_df['iqr_filter'] = group_df['iqr_filter'].fillna(True).astype(bool)  # Non-WT variants are not filtered

        # Filter by IQR
        group_df = group_df[group_df['iqr_filter']]

        # Re-intersect the dataframes
        initial_df = pd.merge(initial_df, group_df[['position', 'variant', 'umi']], on=['position', 'variant', 'umi'],
                              how='inner').fillna(0)
        group_df = pd.merge(group_df, initial_df[['position', 'variant', 'umi']], on=['position', 'variant', 'umi'],
                            how='inner').fillna(0)

    # Next, collapse all the umis to variants
    initial_df = initial_df.groupby(['position', 'variant', 'residue', 'is_control']).agg(
        {'count': 'sum'}).reset_index()
    group_df = group_df.groupby(['position', 'variant', 'residue', 'is_control']).agg({'count': 'sum'}).reset_index()

    # Compute the count filter for the treatment group (at least 8 counts)
    # after collapsing the umis because of how degenerate they are
    # count_filter = (group_df['count'].values >= counts_cutoff)
    # initial_df = initial_df[count_filter]
    # group_df = group_df[count_filter]

    # Next, bootstrap the counts
    initial_dist = bootstrap_initial_distribution(initial_df, n=n)
    group_dist = bootstrap_initial_distribution(group_df, n=n)

    # Compute the normality of the distributions using the Shapiro-Wilk test
    # initial_normality = np.array([scipy.stats.shapiro(dist)[1] for dist in initial_dist])
    # group_normality = np.array([scipy.stats.shapiro(dist)[1] for dist in group_dist])

    fitness_dist = _compute_fitness(group_dist,
                                    group_dist[group_df.variant == "REF", :],
                                    initial_dist,
                                    initial_dist[initial_df.variant == "REF", :],
                                    dt, normalize)
    observed_fitness = _compute_fitness(group_df['count'],
                                        group_df[group_df.variant == "REF"]['count'],
                                        initial_df['count'],
                                        initial_df[initial_df.variant == "REF"]['count'],
                                        dt, normalize)

    # Compute the 95% Confidence interval and see if 0 is within it
    # I.e. significant if 0 is not within the 95% CI
    lower_bound = np.nanpercentile(fitness_dist, 2.5, axis=1)
    upper_bound = np.nanpercentile(fitness_dist, 97.5, axis=1)
    is_significant = (0 < lower_bound) | (0 > upper_bound)
    # Set the fitness to anything non significant to 0
    fitness_dist[~is_significant] = 0
    # Now we have the bootstrapped fitness distributions (n_variants x n_bootstraps)
    fitness_standard_error = np.nanstd(fitness_dist, axis=1, ddof=1)
    # fitness_mean = np.nanmean(fitness_dist, axis=1)
    fitness_median = np.nanmedian(fitness_dist, axis=1)

    # Get the p-value of the observed fitness
    p_value = 1 - np.nanmean(np.abs(fitness_dist) < np.abs(observed_fitness[:, np.newaxis]), axis=1)

    fitness_df = pd.DataFrame({
        'position': group_df['position'],
        'variant': group_df['variant'],
        'residue': group_df['residue'],
        'is_control': group_df['is_control'],
        'count': group_df['count'],
        # 'bootstrap_normality_p': group_normality,
        'initial_count': initial_df['count'],
        # 'initial_bootstrap_normality_p': initial_normality,
        'observed_fitness': observed_fitness,
        'median_fitness': fitness_median,
        'standard_error': fitness_standard_error,
        'p_value': p_value,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
    })
    fitness_df['experiment'] = experiment_name
    return fitness_df


### Plotting
def shannon_diversity_index(counts, weights=None):
    """
    Calculate the Shannon Diversity Index for a given list of species counts. This is normalized by the number of species.

    Parameters:
    counts (list or np.array): A list or array of species counts.
    weights (list or np.array): A list or array of weights for each species.

    Returns:
    float: The Shannon Diversity Index.
    """
    if weights is not None:
        # Normalize the weights
        weights = weights / np.sum(weights)
    # Convert to proportions
    proportions = counts / np.sum(counts)
    # Compute the weighted Shannon Diversity Index
    if weights is not None:
        return entropy(proportions, base=2, qk=weights) / (np.log2(len(counts)) + 1e-6)
    else:
        return entropy(proportions, base=2) / (np.log2(len(counts)) + 1e-6)


def weighted_std(values, weights):
    # Set nans to 0 weight
    weights = np.nan_to_num(weights)
    weights = weights / np.sum(weights)
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)
    return np.sqrt(variance)


def weighted_avg(values, weights):
    # Set nans to 0 weight
    weights = np.nan_to_num(weights)
    weights = weights / np.sum(weights)
    return np.average(values, weights=weights, axis=0)


def summarize_data(total_output: pd.DataFrame,
                   ordering: list[str] = None,
                   weighted: bool = False,
                   infill_positions: list[str] = list(),
                   keep_variants_independent: bool = False,
                   keep_residues_independent: bool = False,
                   compute_sdi: bool = False,
                   significance_threshold: float = 0.05):  ## Polish data

    grouping_cols = ['position', 'is_control']
    if keep_residues_independent or keep_variants_independent:
        grouping_cols.append('residue')

    # Set minimum p-value
    # total_output['p_value'] = total_output['p_value'].clip(lower=1e-4)
    # Force bootstrapping to NaN if observed is nan
    total_output['standard_error'] = np.where(total_output['observed_fitness'].isna(), np.nan,
                                              total_output['standard_error'])
    total_output['lower_bound'] = np.where(total_output['observed_fitness'].isna(), np.nan, total_output['lower_bound'])
    total_output['upper_bound'] = np.where(total_output['observed_fitness'].isna(), np.nan, total_output['upper_bound'])
    total_output['p_value'] = np.where(total_output['observed_fitness'].isna(), np.nan, total_output['p_value'])
    # Compute number of variants for each experiment/position
    total_output['n'] = total_output.groupby(grouping_cols + ['experiment'])['count'].transform('count')

    # Weights
    if 'weight' not in total_output.columns and weighted:
        total_output['weight'] = total_output['n']

    # First, compute the average wrt the mouse
    # Note, to aggregate the p-values, we will use Fisher's method

    def aggregate(group):
        # Aggregate the p-values using Fisher's method
        # p_values = combine_pvalues(group['p_value'], method='fisher')
        # NOTE: The p-values are not independent, so we can't use Fisher's method, instead we will use the Harmonic Mean P-value
        # https://en.wikipedia.org/wiki/Extensions_of_Fisher%27s_method

        if weighted:
            # Note that we are using median in unweighted and mean in weighted
            fitness_std = weighted_std(group['observed_fitness'], group['weight'])
            fitness_std = np.nan_to_num(fitness_std, nan=0.01)
            avg_fitness = weighted_avg(group['observed_fitness'], group['weight'])
            # avg_standard_error = weighted_avg(group['standard_error'], group['weight'])
            standard_error = fitness_std / np.sqrt(group.shape[0])  # SEM
            # p_values = hmean(group['p_value'], weights=group['weight'])
            # Compute SDI wrt observed counts
            sdi = shannon_diversity_index(group['count'], weights=group['weight']) if compute_sdi else None
        else:
            fitness_std = group['observed_fitness'].std()
            fitness_std = np.nan_to_num(fitness_std, nan=0.01)
            avg_fitness = np.nanmedian(group['observed_fitness'].values)
            # avg_standard_error = group['standard_error'].median()
            standard_error = fitness_std / np.sqrt(group.shape[0])  # SEM
            # p_values = hmean(group['p_value'])
            # Compute SDI wrt observed counts
            sdi = shannon_diversity_index(group['count']) if compute_sdi else None

        # Compute a p-value to test if the fitness is significantly different from 0
        # if np.all(group['observed_fitness'] == 0.):
        #     p = 1.0
        # else:
        #     p = scipy.stats.wilcoxon(group['observed_fitness'].values, alternative='two-sided')[1]

        to_append = {
            'mouse': group['experiment'].iloc[0],
            'observed_fitness': avg_fitness,
            'standard_error': standard_error,
            # 'p_value': p_values,
            'observed_fitness_std': fitness_std,
            'weight': group['weight'].sum() if weighted else 1.,
            'sdi': sdi if compute_sdi else 0.0,
            'n': group.shape[0],
            'count': group['count'].sum(),
            'initial_count': group['initial_count'].sum(),
            'updated_residue': group['residue'].iloc[0] if group['residue'].nunique() == 1 else 'X',
            # 'p_value': p,
        }

        if keep_variants_independent:
            to_append['variant'] = group['variant'].iloc[0]
        else:
            to_append['variant'] = group['variant'].iloc[0] if group['is_control'].iloc[0] else (
                        group['variant'].iloc[0][:-1] + 'X')

        return pd.Series(to_append)

    # Aggregate variants to positions/experiments
    # if groupby_mouse:
    #     total_output_by_experiment = total_output.groupby(['position', 'umi']).apply(aggregate).reset_index()
    # else:
    total_output_by_experiment = total_output.groupby(grouping_cols + ['experiment']).apply(aggregate).reset_index()
    total_output_by_experiment['residue'] = total_output_by_experiment['updated_residue']
    # Aggregate to positions
    if keep_variants_independent:
        total_output_by_position = total_output_by_experiment.groupby(grouping_cols).apply(aggregate).reset_index()
    else:
        total_output_by_position = total_output_by_experiment.groupby(['position', 'is_control']).apply(
            aggregate).reset_index()
    total_output_by_position['residue'] = total_output_by_position['updated_residue']

    # total_output_by_experiment['fdr'] = multipletests(total_output_by_experiment['p_value'].values, method='fdr_bh')[1]
    # total_output_by_position['fdr'] = multipletests(total_output_by_position['p_value'].values, method='fdr_bh')[1]

    # Count variants with a significant p-value per position
    n_significant = (
        total_output[total_output['p_value'] < significance_threshold]
        .groupby('position')['variant']
        .nunique()
        .rename('n_significant')
    )
    total_output_by_position = total_output_by_position.merge(
        n_significant, on='position', how='left'
    )
    total_output_by_position['n_significant'] = total_output_by_position['n_significant'].fillna(0).astype(int)

    min_fitness = total_output_by_position['observed_fitness'].abs().min() * 1.01
    # In-fill missing data from various sources to capture complete elimination of a variant
    missing_positions_in_total_output = set(infill_positions) - set(total_output_by_position['variant'])
    # Add missing positions to total_output_by_position (min observed fitness - 0.1)
    missing_positions_df = pd.DataFrame({
        'position': list(missing_positions_in_total_output),
        'variant': list(missing_positions_in_total_output),
        'residue': [p if not p.endswith('X') else 'X' for p in missing_positions_in_total_output],
        'updated_residue': [p if not p.endswith('X') else 'X' for p in missing_positions_in_total_output],
        'is_control': [not m.endswith('X') for m in missing_positions_in_total_output],
        'observed_fitness': min_fitness,
        'standard_error': 0.0,
        'observed_fitness_std': 0.0,
        'weight': 1.,
        'sdi': 0.0,
        'n': 0,
        'count': 0,
        # 'p_value': 1.0,
        # 'fdr': 1.0,
    })
    total_output_by_position = pd.concat([total_output_by_position, missing_positions_df])

    # # Replace standard error by position with average of standard error by experiment
    # error = total_output_by_experiment.groupby('position')['standard_error'].median().reset_index()
    # total_output_by_position = pd.merge(total_output_by_position, error, on='position', suffixes=('', '_experiment'))
    # total_output_by_position['standard_error'] = total_output_by_position['standard_error_experiment']

    if ordering is not None:
        def map_var(var):
            if var in ordering:
                return ordering.index(var)
            elif (var[:-1] + 'X') in ordering:
                return ordering.index(var[:-1] + 'X')
            else:
                raise KeyError(var)

        total_output_by_experiment['position_i'] = total_output_by_experiment['variant'].map(map_var)
        total_output_by_position['position_i'] = total_output_by_position['variant'].map(map_var)
        total_output_by_experiment = total_output_by_experiment.sort_values('position_i', ascending=True)
        total_output_by_position = total_output_by_position.sort_values('position_i', ascending=True)
    else:
        # Sort the data by mean fitness
        total_output_by_experiment = total_output_by_experiment.sort_values('observed_fitness', ascending=True)
        total_output_by_position = total_output_by_position.sort_values('observed_fitness', ascending=True)

    return total_output_by_experiment, total_output_by_position


def make_fitness_heatmap(fitness_tests_df: pd.DataFrame,
                         wt_sequence: str,
                         ordering: list[str] = None,
                         outline_variants: list[Optional[tuple[str, str]]] = None,
                         annotate_significant_fitness=True,
                         annotated: list[str] = None,
                         ax: plt.Axes = None):
    """
    Make a heatmap of fitness data (N_AAs x N_residues).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(100, 20))
    else:
        fig = None

    amino_acids = [
        "R",  # Arginine
        "H",  # Histidine
        "K",  # Lysine
        "D",  # Aspartic Acid
        "E",  # Glutamic Acid
        "S",  # Serine
        "T",  # Threonine
        "N",  # Asparagine
        "Q",  # Glutamine
        "C",  # Cysteine
        "G",  # Glycine
        "P",  # Proline
        "A",  # Alanine
        "V",  # Valine
        "I",  # Isoleucine
        "L",  # Leucine
        "M",  # Methionine
        "F",  # Phenylalanine
        "Y",  # Tyrosine
        "W",  # Tryptophan
    ]

    # Extract the fitness matrix
    N_positions = fitness_tests_df['position'].max()
    N_AAs = len(amino_acids)
    fitness_matrix = np.full((N_AAs, N_positions), np.nan)
    p_value_matrix = np.full((N_AAs, N_positions), np.nan)
    significant_fitness_matrix = np.full((N_AAs, N_positions), np.nan)
    wt_mask = np.zeros((N_AAs, N_positions))
    outline_mask = np.full((N_AAs, N_positions), "")
    positions = []
    for i, row in fitness_tests_df.iterrows():
        position = row['position']
        positions.append(position)
        alt = row['alt']
        if alt not in amino_acids:
            continue
        fitness = row['observed_fitness']
        p_value = row['t_fdr'] if not np.isnan(row['t_fdr']) else 1.0
        fitness_p_value = row['fdr'] if not np.isnan(row['fdr']) else 1.0
        fitness_matrix[amino_acids.index(alt), position - 1] = fitness
        p_value_matrix[amino_acids.index(alt), position - 1] = p_value
        significant_fitness_matrix[amino_acids.index(alt), position - 1] = fitness_p_value < 0.05
        if outline_variants is not None:
            for outline in outline_variants:
                if row['ref'] + str(position) + alt == outline[0]:
                    outline_mask[amino_acids.index(alt), position - 1] = outline[1]

    for pos, aa in enumerate(wt_sequence):
        wt_mask[amino_acids.index(aa), pos] = 1

    not_detected_mat = np.isnan(fitness_matrix)

    # Plot the labelled heatmap of fitnesses (y-axis is residue, x-axis is position)
    if ordering is not None:
        new_fitness_mat = np.full((N_AAs, len(ordering)), np.nan)
        new_p_value_mat = np.full((N_AAs, len(ordering)), np.nan)
        new_not_detected_mat = np.full((N_AAs, len(ordering)), True)
        new_wt_mask = np.full((N_AAs, len(ordering)), 0.0)
        new_significant_fitness_mat = np.full((N_AAs, len(ordering)), np.nan)
        outline_mask = np.full((N_AAs, len(ordering)), "")
        new_positions = []
        for i, var in enumerate(ordering):
            try:
                position = int(var) - 1
                new_positions.append(position + 1)
            except:
                try:
                    position = int(var[1:-1]) - 1
                    new_positions.append(position + 1)
                except:
                    new_positions.append(var)
                    continue
            if position >= N_positions:
                continue
            new_fitness_mat[:, i] = fitness_matrix[:, position]
            new_p_value_mat[:, i] = p_value_matrix[:, position]
            new_not_detected_mat[:, i] = not_detected_mat[:, position]
            new_wt_mask[:, i] = wt_mask[:, position]
            new_significant_fitness_mat[:, i] = significant_fitness_matrix[:, position]
            outline_mask[:, i] = outline_mask[:, position]

        fitness_matrix = new_fitness_mat
        p_value_matrix = new_p_value_mat
        not_detected_mat = new_not_detected_mat
        wt_mask = new_wt_mask
        significant_fitness_matrix = new_significant_fitness_mat
        positions = new_positions
        N_positions = len(ordering)
    else:
        ordering = list(wt_sequence)

    fitness_cmap = mpl.colormaps.get_cmap('coolwarm')
    # fitness_norm = mpl.colors.TwoSlopeNorm(vmin=np.nanmin(fitness_matrix.flatten()),
    #                                        vmax=np.nanmax(fitness_matrix.flatten()),
    #                                        vcenter=0)
    fitness_norm = mpl.colors.TwoSlopeNorm(vmin=-.1,
                                           vmax=.1,
                                           vcenter=0)
    allwhite_grey_cmap = mpl.colors.ListedColormap([[1, 1, 1], [0.9, 0.9, 0.9]])
    allwhite_grey_norm = mpl.colors.BoundaryNorm(boundaries=[0, 0.5, 1], ncolors=2)
    all_darkgrey_cmap = mpl.colors.ListedColormap([[1, 1, 1], [0.2, 0.2, 0.2]])
    all_darkgrey_norm = mpl.colors.BoundaryNorm(boundaries=[0, 0.5, 1], ncolors=2)
    # colors = fitness_cmap(fitness_norm(fitness_mat.flatten()))

    # Plot the heatmap
    sns.heatmap(fitness_matrix, ax=ax, cmap=fitness_cmap, norm=fitness_norm,
                yticklabels=amino_acids, xticklabels=ordering, cbar=False,
                linewidths=0.5, linecolor='black', annot=False, mask=not_detected_mat == 1)
    # Add a colorbar by making a fake plot
    # plt.colorbar(plt.cm.ScalarMappable(norm=fitness_norm, cmap=fitness_cmap),
    #              ax=ax, pad=0.01, shrink=0.5, label='Observed Fitness')
    # Annotate not detected with gray color (I.e. only color the boxes that are not detected)
    sns.heatmap(not_detected_mat, cbar=False, ax=ax, linewidths=0.5, linecolor='black',
                yticklabels=amino_acids, xticklabels=ordering, mask=not_detected_mat == 0, cmap=allwhite_grey_cmap,
                norm=allwhite_grey_norm)
    # Annotate WT with dark grey color
    sns.heatmap(wt_mask, cbar=False, ax=ax, linewidths=0.5, linecolor='black',
                yticklabels=amino_acids, xticklabels=ordering, mask=wt_mask == 0, cmap=all_darkgrey_cmap,
                norm=all_darkgrey_norm)

    if outline_variants is not None:
        # Transform the p-values to a color
        # Reverse plasma cmap
        to_add = []
        for i in range(N_positions):
            for j in range(N_AAs):
                if outline_mask[j, i] == "":
                    continue
                rect = plt.Rectangle((i, j), 1, 1, linewidth=2,
                                     edgecolor=outline_mask[j, i],
                                     facecolor='none', alpha=0.75)
                to_add.append(rect)
        for rect in to_add:  # Draw insignificant first
            ax.add_patch(rect)

    if annotated is not None:
        annotated = [var[1:] for var in annotated]  # Remove the starting AA since it doesn't matter

    # Annotate Positions
    for i in range(N_AAs):
        for j in range(N_positions):
            # if wt_mask[i, j] == 1:
            #     ax.text(j+0.5, i+0.5, "+", ha='center', va='center', color='black', fontsize=10)
            if annotate_significant_fitness:
                if significant_fitness_matrix[i, j] == 1:
                    if annotated is not None:
                        # Construct the variant name
                        possible_var_name = f"{positions[j]}{amino_acids[i]}"
                        if possible_var_name in annotated:
                            # Note we need to specify the default serif font because archivo doesn't contain this glyph
                            ax.text(j + 0.5, i + 0.5, "\u272a", ha='center', va='center', color='black', fontsize=18,
                                    fontdict={'family': 'sans-serif'}, rotation=360 - 15)
                            continue
                    ax.text(j + 0.5, i + 0.5, "\u2605", ha='center', va='center', color='black', fontsize=18,
                            fontdict={'family': 'sans-serif'}, rotation=360 - 15)

    # Shift x-limits by -0.5 so heatmap cells (spanning [i, i+1]) center at integer i,
    # matching matplotlib's categorical barplot which centers bars at integer positions.
    ax.set_xlim(-0.5, N_positions - 0.5)

    # Move the x-axis labels to the top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    # Rotate the x-axis labels
    ax.set_xticks(np.arange(len(ordering)) + 0.5)
    ax.set_xticklabels(ordering, rotation=90 if any([len(var) > 2 for var in ordering]) else 0)
    ax.set_yticklabels(amino_acids, rotation=0)

    return ax


def make_fitness_df_matrix(total_output: pd.DataFrame) -> tuple[np.array, list[str], list[str]]:
    no_control_output = total_output[~total_output['variant'].str.contains('_')].copy()
    no_control_output['position'] = pd.to_numeric(no_control_output['position'], errors='coerce')
    no_control_output = no_control_output.dropna(subset=['position'])
    n_positions = int(no_control_output['position'].max())
    residues = no_control_output['residue'].unique().tolist()
    n_residues = len(residues)
    experiments = no_control_output.experiment.unique().tolist()
    n_experiments = len(experiments)

    fitness = np.full((n_positions, n_residues, n_experiments), np.nan)
    for i in range(n_positions):
        position_df = no_control_output[no_control_output['position'] == (i + 1)]
        for j, residue in enumerate(residues):
            residue_df = position_df[position_df['variant'].str.endswith(residue)]
            if residue_df.shape[0] > 0:
                for k, experiment in enumerate(experiments):
                    experiment_df = residue_df[residue_df['experiment'] == experiment]
                    if experiment_df.shape[0] > 0:
                        fitness[i, j, k] = np.nanmedian(experiment_df['observed_fitness'])

    return fitness, residues, experiments


def test_residues(
        fitness_matrix: np.array,  # N_positions x N_residues x N_replicates
        residues: list[str],  # Residue names
        compare_to_ref: bool = False  # If given, compare to ref (0 fitness)
) -> pd.DataFrame:
    """
    Compute the significance of each residue at each position by first running an ANOVA test,
    then a post-hoc series of t-tests.
    other residues if not given.
    :param fitness_matrix: The matrix of fitnesses. N_positions x N_residues x N_replicates
    :param residues: The list of residues names.
    :param compare_to_ref: If given, compare to the reference fitness (0).
    :return: A dataframe with the p-values of the ANOVA test and the post-hoc tests.
    """
    N_positions, N_residues, N_replicates = fitness_matrix.shape

    results = defaultdict(list)
    for i in range(N_positions):
        pos = i + 1
        position_data = fitness_matrix[i, :, :]
        comparison_data = []
        for res_idx in range(N_residues):
            res_data = position_data[res_idx, :]
            filtered_res = res_data[~np.isnan(res_data)]
            if filtered_res.size > 0:
                comparison_data.append(filtered_res)
            else:
                comparison_data.append([])

        if compare_to_ref:
            for res_idx in range(N_residues):
                if len(comparison_data) == 0 or len(comparison_data[res_idx]) < 3:
                    results['position'].append(pos)
                    results['residue'].append(residues[res_idx])
                    results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                    results['statistic'].append(np.nan)
                    results['p_value'].append(np.nan)
                    continue

                # Test to see if the fitnesses are significantly different from 0
                res = scipy.stats.ttest_1samp(comparison_data[res_idx], 0, alternative='two-sided')
                statistic = res.statistic
                p_val = res.pvalue
                results['position'].append(pos)
                results['residue'].append(residues[res_idx])
                results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                results['statistic'].append(statistic)
                results['p_value'].append(p_val)

        else:
            # Remove all empty comparisons
            comparison_data = [data for data in comparison_data if len(data) > 0]
            if len(comparison_data) <= 1:
                for res_idx in range(N_residues):
                    results['position'].append(pos)
                    results['residue'].append(residues[res_idx])
                    results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                    results['anova_stat'].append(np.nan)
                    results['anova_pvalue'].append(np.nan)
                    results['t_stat'].append(np.nan)
                    results['t_pvalue'].append(np.nan)
                continue

            # ANOVA
            stat, pvalue = scipy.stats.f_oneway(*comparison_data)

            if pvalue < 0.05:
                # One vs rest t-tests
                for res_idx in range(N_residues):
                    values = position_data[res_idx, :]
                    other_values = np.delete(position_data, res_idx, axis=0).flatten()
                    # Filter NaNs
                    values = values[~np.isnan(values)]
                    other_values = other_values[~np.isnan(other_values)]
                    if values.size <= 1 or other_values.size <= 1:
                        results['position'].append(pos)
                        results['residue'].append(residues[res_idx])
                        results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                        results['anova_stat'].append(stat)
                        results['anova_pvalue'].append(pvalue)
                        results['t_stat'].append(np.nan)
                        results['t_pvalue'].append(np.nan)
                        continue

                    t_stat, t_p = scipy.stats.ttest_ind(values, other_values, equal_var=False, alternative='two-sided')
                    results['position'].append(pos)
                    results['residue'].append(residues[res_idx])
                    results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                    results['anova_stat'].append(stat)
                    results['anova_pvalue'].append(pvalue)
                    results['t_stat'].append(t_stat)
                    results['t_pvalue'].append(t_p)
            else:
                for res_idx, res in enumerate(residues):
                    results['position'].append(pos)
                    results['residue'].append(res)
                    results['observed_fitness'].append(np.nanmedian(position_data[res_idx, :]))
                    results['anova_stat'].append(stat)
                    results['anova_pvalue'].append(pvalue)
                    results['t_stat'].append(np.nan)
                    results['t_pvalue'].append(np.nan)

    results = pd.DataFrame(results)
    # Correct for multiple testing
    if compare_to_ref:
        pvals = results[~results['p_value'].isna()]['p_value'].values
    else:
        pvals = results[~results['t_pvalue'].isna()]['t_pvalue'].values

    adjusted = multipletests(pvals, method='fdr_bh')[1]
    if compare_to_ref:
        results['fdr'] = np.nan
        results.loc[~results['p_value'].isna(), 'fdr'] = adjusted
    else:
        results['t_fdr'] = np.nan
        results.loc[~results['t_pvalue'].isna(), 't_fdr'] = adjusted

    return results


def make_fitness_barplot(total_output: pd.DataFrame,  # Data to plot
                         secondary_output: pd.DataFrame = None,
                         # If given, plot the secondary data as well using the same order as the primary output
                         additional_outputs: list[tuple[str, pd.DataFrame]] = None,
                         # If given, plot the additional data as well using the same order as the primary output Displayed as a heatmap
                         title: Optional[str] = "Observed Fitnesses",  # Title of the plot
                         secondary_output_title: Optional[str] = None,  # Title of the secondary plot
                         output_name: str = None,  # If given, save plot to this file (in svg and pdf)
                         intermediate_output_file: str = None,
                         # If given, save all the raw data used to make the plot into an excel file
                         weighted: bool = False,  # If given, weighted average using these weights per position
                         ordering: list[str] = None,  # If given, order the positions according to this list
                         tracks: list[tuple[str, Union[dict[str, str], Callable[[Any, list[str]], None]]]] = None,
                         # If given, add tracks below the x axis for each position (with a given name per track)
                         plot_points: bool = True,
                         plot_by_residues: bool = False,
                         # If true, tests are performed to see if the fitnesses are significantly different from expectation
                         split_by_residue: bool = False,  # Have each residue as a separate point
                         scale_points_by_sem: bool = False,  # If true, scale point size by the standard error
                         color_points: bool = False,
                         color_by_sdi: bool = False,
                         color_by_residue: bool = False,
                         wt_sequence: str = None,
                         # If given, use this sequence to determine the WT amino acid at each position
                         include_heatmap: bool = True,  # If true, include a heatmap of fitness by residue
                         highlighted_positions: list[str] = None,
                         ax: plt.Axes = None,
                         ylim: Optional[tuple[float, float]] = None,
                         symmetrical_y_axis: bool = True,
                         important_variants: list[str] = None,
                         # List of variants to annotate if there is significant fitness change
                         significance_threshold: float = 0.05,
                         plot_singletons_as_letters: bool = False,
                         n_significant_by_variants: bool = True
                         ):
    # Positions that may have dropped out and should be added
    additional_positions = set()
    if secondary_output is not None:
        additional_positions.update(set(secondary_output['variant']))
    if highlighted_positions is not None:
        additional_positions.update(set(highlighted_positions))
    if ordering is not None:
        additional_positions.update(set(ordering))
    if additional_outputs is not None:
        # Drop positions not in the total output or the additional positions
        updated_additional_outputs = []
        for name, output in additional_outputs:
            positions_to_remove = set(output.position.values) - set(total_output.position.values) - additional_positions
            output = output[~output.position.isin(positions_to_remove)].copy()
            updated_additional_outputs.append((name, output))
        additional_outputs = updated_additional_outputs
        # print([o[1].shape for o in additional_outputs])

    additional_positions = list(additional_positions)

    if plot_by_residues or include_heatmap:
        fitness_matrix, residues, experiments = make_fitness_df_matrix(total_output)
        fitness_test_df = test_residues(fitness_matrix, residues, compare_to_ref=False)  # Test for outlier residues
        wt_fitness_test_df = test_residues(fitness_matrix, residues, compare_to_ref=True)  # Test for difference to WT
        fitness_test_df = pd.merge(fitness_test_df, wt_fitness_test_df, on=['position', 'residue', 'observed_fitness'])
        # Add the reference residue to the fitness test
        if wt_sequence:
            ref_residue_df = pd.DataFrame({
                'position': [i + 1 for i in range(len(wt_sequence))],
                'ref': list(wt_sequence),
            })
        else:
            # Infer the WT residue from the dataframe
            for i in range(fitness_matrix.shape[0]):
                pos_df = total_output[total_output['position'] == (i + 1)]
                if pos_df.shape[0] > 0:
                    ref_residue = pos_df['variant'].str[0].values[0]
                else:
                    ref_residue = 'X'
                residues[i] = ref_residue
            ref_residue_df = pd.DataFrame({
                'position': [i + 1 for i in range(fitness_matrix.shape[0])],
                'ref': residues,
            })
        fitness_test_df = fitness_test_df.rename(columns={
            'residue': 'alt'
        })
        fitness_test_df = pd.merge(fitness_test_df, ref_residue_df, on='position', how='left')

        fitness_test_df['variant'] = fitness_test_df['ref'] + fitness_test_df.position.astype(str) + \
                                     fitness_test_df['alt']
    else:
        fitness_test_df = None
    total_output_by_experiment, total_output_by_position = summarize_data(total_output,
                                                                          ordering=ordering,
                                                                          weighted=weighted, compute_sdi=color_by_sdi,
                                                                          infill_positions=additional_positions,
                                                                          # keep_variants_independent=True,
                                                                          keep_residues_independent=plot_by_residues,
                                                                          significance_threshold=significance_threshold)

    if n_significant_by_variants and fitness_test_df is not None:
        _heatmap_aas = {"K","R","H","Y","F","W","L","I","V","P","G","A","M","C","S","T","N","Q","E","D"}
        _n_sig_by_pos = {}
        for _, _row in fitness_test_df.iterrows():
            if _row['alt'] not in _heatmap_aas:
                continue
            _fdr = _row['fdr'] if not pd.isna(_row['fdr']) else 1.0
            if _fdr < significance_threshold:
                _pos = int(_row['position'])
                _n_sig_by_pos[_pos] = _n_sig_by_pos.get(_pos, 0) + 1
        def _lookup_n_sig(p):
            try:
                return _n_sig_by_pos.get(int(float(p)), 0)
            except (ValueError, TypeError):
                return 0
        total_output_by_position['n_significant'] = total_output_by_position['position'].apply(_lookup_n_sig).astype(int)

    # total_output_by_position = total_output_by_position.sort_values('observed_fitness', ascending=True)
    # total_output_by_position = pd.concat([total_output_by_position.head(250), total_output_by_position.tail(250)])
    # total_output_by_experiment = total_output_by_experiment[total_output_by_experiment.variant.isin(total_output_by_position.variant)]

    ## Get constants
    N_experiments = total_output_by_experiment['experiment'].nunique()
    n_positions = total_output_by_position['variant'].nunique()

    if intermediate_output_file is not None:
        writer = pd.ExcelWriter(intermediate_output_file)
        total_output_by_experiment.to_excel(writer, sheet_name='Data Points', index=False)
        total_output_by_position.to_excel(writer, sheet_name='Fully Aggregated Data', index=False)
        if fitness_test_df is not None:
            fitness_test_df.to_excel(writer, sheet_name='Fitness Tests', index=False)
        writer.close()

    ## Plotting
    # Bar plot of fitnesses with each positions as dots as well
    if ax is None:
        fig, ax = plt.subplots(figsize=(100, 60 if secondary_output is not None else 40))
    else:
        fig = None

    # Add title
    if title is not None:
        ax.set_title(title)
    if ylim is None:
        if symmetrical_y_axis:
            ylim = [-total_output_by_position.observed_fitness.abs().max() * 1.25,
                    total_output_by_position.observed_fitness.abs().max() * 1.25]
        else:
            ylim = [total_output_by_experiment.observed_fitness.min() * 1.1,
                    total_output_by_experiment.observed_fitness.max() * 1.1]
    # Set the ylimit
    ax.set_ylim(ylim)
    # Add labels
    ax.set_xlabel("Amino Acid Substitution")
    ax.set_ylabel("Log2 Fold-Change in Relative Abundance")

    # First plot the the barplot of observed mean fitness colored by the values themselves
    mean_obs_fitness_cmap = mpl.colormaps.get_cmap('coolwarm')
    mean_obs_fitness_range = mpl.colors.TwoSlopeNorm(vmin=-0.1, vmax=0.1, vcenter=0)
    # mean_obs_fitness_range = mpl.colors.TwoSlopeNorm(vmin=min(total_output_by_position['observed_fitness'].min(), -.1), vmax=max(total_output_by_position['observed_fitness'].max(), 0.1), vcenter=0)
    mean_obs_fitness_colors = mean_obs_fitness_cmap(
        mean_obs_fitness_range(np.nan_to_num(total_output_by_position['observed_fitness'].values)))
    bar_max_n_sig = total_output_by_position['n_significant'].max()
    mean_obs_fitness_colors[:, 3] = total_output_by_position['n_significant'].fillna(0).apply(
        lambda n: 0.25 + 0.75 * (n / bar_max_n_sig) if bar_max_n_sig > 0 else 0.25
    ).values
    ax.bar(total_output_by_position['variant'], total_output_by_position['observed_fitness'],
           color=mean_obs_fitness_colors, edgecolor='black', zorder=1)
    # Next add error bars
    ax.errorbar(total_output_by_position['variant'], total_output_by_position['observed_fitness'],
                yerr=total_output_by_position['standard_error'],
                elinewidth=1, capsize=1.5, ls='none', color='black',
                alpha=.75, zorder=2)  # Error bar for the mean fitness
    # Add x-axis tick labels
    ax.set_xticks(total_output_by_position['variant'])
    ax.set_xticklabels(total_output_by_position['variant'], rotation=90, fontsize=6)

    if plot_points:
        # ax.errorbar(total_output['position'], total_output['observed_fitness'],
        #             # Use lower and upper bound
        #             yerr=total_output['standard_error'],
        #             elinewidth=.5, capsize=.25, ls='none', color='black', zorder=3)  # Error bar for the observed fitnesses
        # Add the individual points
        # point_scale = (np.abs(np.log10(total_output_by_experiment['p_value'].values + 1.))+1)*5  # Scale the points by the p-value
        # Scale by number of variants
        # point_scale = (total_output_by_experiment['n'] / total_output_by_experiment['n'].max()) * 10
        if scale_points_by_sem:
            # Scale by inverse of SEM
            total_output_by_experiment['point_scale'] = (1 / (total_output_by_experiment['standard_error'] + .05))
        else:
            # Scale by number of variants within a position
            total_output_by_experiment['point_scale'] = 25
            # total_output_by_experiment['point_scale'] = total_output_by_experiment.groupby(['variant'])['n'].transform(lambda x: x / x.max()) * 20

        AA_COLORMAP = dict(
            T='#0E019D',
            Y='#9B9C00',
            W='#CF6100',
            V='#9C0057',
            S='#008247',
            R='#CFCF8E',
            N='#948ECF',
            M='#8FB5A4',
            L='#CECF12',
            Q='#CFA784',
            P='#CF84AE',
            K='#D06507',
            I='#CF0878',
            H='#2212CF',
            G='#23B574',
            F='#FDFF47',
            E='#FF973C',
            D='#FF3BAA',
            C='#5948FF',
            A='#59E6A6',
            X='#454544'
        )

        default = '#333333'

        def aa_color_mapping(row):
            return mcolors.to_rgba(AA_COLORMAP.get(row['residue'], default))

        if color_points:
            if color_by_sdi:
                sdi_cmap = mpl.colormaps.get_cmap('viridis')
                sdi_range = mpl.colors.Normalize(vmin=total_output_by_experiment['sdi'].min(),
                                                 vmax=total_output_by_experiment['sdi'].max())
                total_colors = lambda df: sdi_cmap(sdi_range(df['sdi']))
            elif color_by_residue or split_by_residue:
                total_colors = lambda df: df.apply(aa_color_mapping, axis=1)
            else:
                total_colors = lambda df: mean_obs_fitness_cmap(mean_obs_fitness_range(df['observed_fitness']))
        else:
            total_colors = lambda df: mcolors.to_rgba("#333333")
        # if position_weights:
        #     total_colors[:, 3] = total_output[position_weights + '_x'].values * .75
        # else:
        #     total_colors[:, 3] = .75
        if fitness_test_df is not None and plot_by_residues or split_by_residue:
            control_data = total_output_by_experiment[total_output_by_experiment['is_control']]
            ax.scatter(control_data['variant'], control_data['observed_fitness'], s=15,
                       alpha=.5, edgecolors='black', linewidths=.15,
                       marker='o', color=total_colors(control_data), zorder=2)

            if split_by_residue:
                for res in total_output_by_experiment[~total_output_by_experiment['is_control']]['residue'].unique():
                    filtered_variant = total_output_by_experiment[total_output_by_experiment['residue'] == res]
                    ax.scatter(filtered_variant['variant'].values,
                               filtered_variant['observed_fitness'].values,
                               s=15, alpha=.5, edgecolors='black', linewidths=.15,
                               marker="$\\text{" + res + "}$", color=total_colors(filtered_variant), zorder=2)
            else:
                fitness_test_df['is_significant'] = (~fitness_test_df['t_fdr'].isna()) & (
                            fitness_test_df['t_fdr'] < significance_threshold)
                # filtered_test_df['point_scale'] = -np.log10(filtered_test_df['t_fdr'])
                # filtered_test_df['point_scale'] = filtered_test_df['point_scale'].fillna(np.nanmax(filtered_test_df['point_scale'].values))
                fitness_test_df['point_scale'] = 20
                fitness_test_df['residue'] = fitness_test_df['alt']
                # Replace variant with X to map correctly
                fitness_test_df['variant'] = fitness_test_df['ref'] + fitness_test_df['position'].astype(str) + "X"

                insignificant_color = (*mcolors.to_rgba(default)[:3], .25)
                for alt in fitness_test_df['alt'].unique():
                    filtered_alt = fitness_test_df[fitness_test_df['alt'] == alt]
                    color = mcolors.to_rgba(AA_COLORMAP.get(alt, default))
                    ax.scatter(filtered_alt['variant'].values,
                               filtered_alt['observed_fitness'].values,
                               s=20,
                               edgecolors='black', linewidths=.15,  # Add the text
                               color=[(color if row['is_significant'] else insignificant_color) for i, row in
                                      filtered_alt.iterrows()],
                               zorder=3,
                               marker="$\\text{" + alt + "}$")
        else:
            max_n_sig = total_output_by_position['n_significant'].max()
            sig_alpha_map = total_output_by_position.set_index('variant')['n_significant'].apply(
                lambda n: 0.25 + 0.75 * (n / max_n_sig) if max_n_sig > 0 else 0.25
            ).to_dict()

            def colors_with_sig_alpha(df, colors):
                n = len(df)
                if n == 0:
                    return np.empty((0, 4))
                if isinstance(colors, tuple) and len(colors) == 4:
                    rgba = np.tile(np.array(colors, dtype=float), (n, 1))
                else:
                    rgba = np.array([mcolors.to_rgba(c) for c in colors], dtype=float)
                rgba[:, 3] = df['variant'].map(sig_alpha_map).fillna(0.25).values
                return rgba

            singleton_data = total_output_by_experiment[total_output_by_experiment['n'] <= 1]
            nonsingleton_data = total_output_by_experiment[total_output_by_experiment['n'] > 1]
            ax.scatter(nonsingleton_data['variant'], nonsingleton_data['observed_fitness'], s=20,
                       # nonsingleton_data['point_scale'],
                       edgecolors='black', linewidths=.15,
                       marker='o',
                       color=colors_with_sig_alpha(nonsingleton_data, total_colors(nonsingleton_data)), zorder=2)
            for res, group in singleton_data.groupby(['residue']):
                res = res[0]
                ax.scatter(group['variant'], group['observed_fitness'], s=20,  # singleton_data['point_scale'],
                           edgecolors='black', linewidths=.15,
                           marker='o' if (not plot_singletons_as_letters or len(res) > 1) else "$\\text{" + res + "}$",
                           color=colors_with_sig_alpha(group, total_colors(group)), zorder=2)

    if fig is not None:
        # Add the colorbar to the axis
        sm = plt.cm.ScalarMappable(cmap=mean_obs_fitness_cmap, norm=mean_obs_fitness_range)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Fitness')

        if color_by_sdi and plot_points:
            # Add a colorbar for SDI values
            sm = plt.cm.ScalarMappable(cmap=sdi_cmap, norm=sdi_range)
            sm.set_array([])
            fig.colorbar(sm, ax=ax, label='Normalized Shannon Diversity Index')

        # Add a colorbar for alpha values
        # if position_weights:
        #     sm = plt.cm.ScalarMappable(cmap=plt.cm.Greys, norm=plt.Normalize(vmin=0, vmax=1))
        #     sm.set_array([])
        #     fig.colorbar(sm, ax=ax, label='Probability of Inclusion')

    # Disable top and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Disable grid
    ax.grid(False)
    # Make xticks invisible
    # Ensure the original x axis is preserved
    ax.set_xticks(total_output_by_position['variant'])
    ax.set_xticklabels(total_output_by_position['variant'], rotation=90, fontsize=6)

    # Highlight positions by setting x labels bold
    if highlighted_positions is not None:
        for tick in ax.get_xticklabels():
            if tick.get_text() in highlighted_positions:
                tick.set_fontweight('bold')

    # Add tracks
    if tracks is not None or secondary_output is not None or include_heatmap:
        if tracks is None:
            tracks = []

        divider = make_axes_locatable(ax)
        if secondary_output is not None:
            new_ax = divider.append_axes("top", size="100%", pad=.65, sharex=ax)
            make_fitness_barplot(
                total_output=secondary_output,
                secondary_output=None,
                additional_outputs=None,
                title=secondary_output_title,
                secondary_output_title=None,
                output_name=None,
                intermediate_output_file=None,
                weighted=weighted,
                ordering=total_output_by_position['variant'].values.tolist(),
                tracks=None,
                plot_points=plot_points,
                scale_points_by_sem=scale_points_by_sem,
                color_points=color_points,
                color_by_sdi=color_by_sdi,
                highlighted_positions=highlighted_positions,
                ax=new_ax,
                ylim=ylim
            )
            # Disable xaxis for the secondary plot
            # new_ax.set_xticklabels([])
            # new_ax.set_xlabel('')
            # new_ax.set_xticks([])
            # new_ax.spines['bottom'].set_visible(False)
            plt.setp(new_ax.get_xticklabels(), visible=False)
        first = True

        if include_heatmap and fitness_test_df is not None:
            new_ax = divider.append_axes("top", size="50%", pad=1.0, sharex=ax)
            make_fitness_heatmap(fitness_test_df, wt_sequence,
                                 ordering=total_output_by_position['variant'].values.tolist(),
                                 annotated=important_variants,
                                 ax=new_ax)

        if additional_outputs is not None and len(
                additional_outputs) > 0:  # Display additional data as a heatmap track below the barplot
            def heatmap_maker(axis, positions):
                # First, collapse data
                # label2data = dict()
                matrix = np.zeros((len(additional_outputs), len(positions)))
                labels = []
                for i, (name, data) in enumerate(additional_outputs):
                    # Collapse the data
                    collapsed_by_experiment, collapsed_by_position = summarize_data(data, ordering=positions,
                                                                                    weighted=weighted,
                                                                                    infill_positions=additional_positions)
                    # label2data[name] = collapsed_by_position
                    matrix[i, :] = collapsed_by_position['observed_fitness'].values
                    labels.append(name)

                # Set y labels
                im = axis.imshow(matrix, aspect='auto', cmap=mean_obs_fitness_cmap, norm=mean_obs_fitness_range)
                axis.set_yticks(np.arange(len(labels)))
                axis.set_yticklabels(labels)
                # Disable x ticks
                axis.set_xticks([])
                axis.grid(False)
                return im

            tracks = [('', heatmap_maker)] + tracks

        if tracks is not None:
            for track_name, data_per_position in tracks:
                is_colormap = isinstance(data_per_position, dict)
                new_ax = divider.append_axes("bottom", size="5%" if is_colormap else "10%", pad=.65 if first else 0.2,
                                             sharex=ax)
                # new_ax.set_title(track_name, y=-.4)
                # Add title to the right of each bar
                # new_ax.text(.96, .5, track_name, transform=new_ax.transAxes, rotation=0, ha='left', va='center')
                # Add title to the left of each bar
                new_ax.text(-0.03, .5, track_name, transform=new_ax.transAxes, rotation=90, ha='center', va='center')
                if is_colormap:
                    # Add the track data
                    # new_ax.axis('off')
                    # Disable axis for the new ax
                    new_ax.spines['top'].set_visible(False)
                    new_ax.spines['right'].set_visible(False)
                    new_ax.spines['left'].set_visible(False)
                    new_ax.spines['bottom'].set_visible(False)
                    new_ax.tick_params(left=False, right=False, labelleft=False, labelbottom=False, top=False,
                                       bottom=False)
                    new_ax.bar(list(range(n_positions)), np.ones(n_positions),
                               color=[data_per_position[position] for position in
                                      total_output_by_position.variant.values.tolist()], width=1, linewidth=0)
                    plt.setp(new_ax.get_xticklabels(), visible=False)
                else:
                    # Pass the axis and positions
                    data_per_position(new_ax, total_output_by_position.variant.values.tolist())
                first = False

        ax.set_xticklabels(total_output_by_position['variant'], rotation=90, fontsize=6)
        plt.setp(ax.get_xticklabels(), visible=True)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, right=False, labelbottom=True)

    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % 25 != 0:
            label.set_visible(False)  # Show every 25th label

    if fig is not None:
        # Finally, prepare the figure for saving
        fig.tight_layout()

        if output_name is not None:
            plt.savefig(output_name + ".pdf", bbox_inches='tight')
            plt.savefig(output_name + ".svg", bbox_inches='tight')
            plt.savefig(output_name + ".png", bbox_inches='tight')
            plt.clf()
        else:
            return fig
    else:
        return ax


mSTING = "MPYSNLHPAIPRPRGHRSKYVALIFLVASLMILWVAKDPPNHTLKYLALHLASHELGLLLKNLCCLAEELCHVQSRYQGSYWKAVRACLGCPIHCMAMILLSSYFYFLQNTADIYLSWMFGLLVLYKSLSMLLGLQSLTPAEVSAVCEEKKLNVAHGLAWSYYIGYLRLILPGLQARIRMFNQLHNNMLSGAGSRRLYILFPLDCGVPDNLSVVDPNIRFRDMLPQQNIDRAGIKNRVYSNSVYEILENGQPAGVCILEYATPLQTLFAMSQDAKAGFSREDRLEQAKLFCRTLEEILEDVPESRNNCRLIVYQEPTDGNSFSLSQEVLRHIRQEEKEEVTMNAPMTSVAPPPSVLSQEPRLLISGMDQPLPLRTDLI"


def generate_figures(total_output: pd.DataFrame, output_dir: str, name: str):
    print("Generating figures...")
    fig_dir = os.path.join(output_dir, "figures")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # fitness_mat, residues, experiments = make_fitness_df_matrix(total_output)
    # results_df = test_residues(fitness_mat, residues)

    # Make the base fitness barplot
    make_fitness_barplot(total_output, None, None, name,
                         os.path.join(fig_dir, "fitness_barplot"),
                         os.path.join(fig_dir, "fitness_barplot_data.xlsx"),
                         color_points=True,
                         # color_by_residue=False,
                         color_by_residue=True,
                         split_by_residue=True,
                         plot_by_residues=True,
                         symmetrical_y_axis=True, wt_sequence=mSTING)


def generate_abundances(experiment_prefix, barcodes_file, twist_file, reference_names, alt_names, data_dirs,
                        capture_seq, barcode_offset, n_processes):
    if not os.path.exists(experiment_prefix + 'abundance_with_variant.csv'):
        # Find all the fastq file prefixes
        data_files = []
        experiment_names = reference_names + alt_names
        for name in experiment_names:
            all_files = []
            for data_dir in data_dirs:
                files = glob(os.path.join(data_dir, f"{name}*R2*.fastq.gz"))
                all_files.extend(files)
            print(f"Found {len(all_files)} files for {name}")
            data_files += [all_files]
        control_variants = [
            "V146L", "N153S", "V154M", "S161A", "I229A", "Y239S", "N241A", "S365A", "I229A", "RRAA", "RRDD", "ΔCTT",
            "ΔUPR", "mSTING_WT", "STOP"
        ]
        abundance_df = make_counts(
            experiment_names=experiment_names,
            file_names=data_files,
            barcodes_file=barcodes_file,
            twist_file=twist_file,
            spike_ins=control_variants,
            capture_seq=capture_seq,
            barcode_offset=barcode_offset,
            n_processes=n_processes
        )
        print("Number of parsed reads:", abundance_df['count'].sum())
        # Add columns for position and residue
        # Transform the variant to a more consistent format
        abundance_df['variant'] = abundance_df['variant'].str.replace('mSTING_WT', 'REF').str.replace('dUPR',
                                                                                                      'ΔUPR').str.replace(
            'dCTT', 'ΔCTT')
        # Compute which variants are controls
        abundance_df['is_control'] = abundance_df['variant'].isin(control_variants + ['REF'])
        # For position, extract the number between the first and last letter, else do nothing
        abundance_df['position'] = abundance_df['variant'].str.extract(r'(\d+)')
        abundance_df.loc[abundance_df['position'].isna(), 'position'] = abundance_df.loc[
            abundance_df['position'].isna(), 'variant']
        # For residue, extract the last letter if and only if there is a number, else do nothing
        abundance_df['residue'] = abundance_df['variant'].str.extract(r'\d+([A-Z])')
        abundance_df.loc[abundance_df['residue'].isna(), 'residue'] = abundance_df.loc[
            abundance_df['residue'].isna(), 'variant']
        abundance_df['is_control'] = abundance_df['variant'].isin(control_variants + ['REF'])
        # Save the corrected umis
        abundance_df.to_csv(experiment_prefix + 'abundance_with_variant.csv', index=False)

        # Collapse umis
        # abundance_df = abundance_df.drop(columns=['residue', 'variant']).groupby(['experiment', 'position']).sum().reset_index()
        # abundance_df.columns = ['experiment','position', 'count']  # Rename the count column
        #
        # for name in experiment_names:
        #     sub_abundance = abundance_df[abundance_df['experiment'] == name]
        #     # Add 0s for variants that are missing
        #     missing_aas = [aa for aa in abundance_df.position.unique() if aa not in sub_abundance["position"].unique()]
        #     if len(missing_aas) > 0:
        #         abundance_df = pd.concat([abundance_df, pd.DataFrame({"position": missing_aas, "count": [0]*len(missing_aas), "experiment": name})], ignore_index=True)
        #
        # abundance_df.to_csv(experiment_prefix + 'abundance.csv', index=False)


def bootstrap_fitness(abundance_df, total_boot_output, reference_names, alt_names, paired, overwrite_boots, dt,
                      normalize, counts_cutoff, N, n_processes, collapse_control_variants: bool = True):
    if not os.path.exists(total_boot_output) or overwrite_boots:
        # Compute all the comparisons we must do:
        if paired:
            comparisons = [([ref], [alt]) for ref, alt in zip(reference_names, alt_names)]
        else:  # We will compute each alt against all references
            comparisons = [(reference_names, [alt]) for alt in alt_names]

        jobs = []
        for ref, alt in comparisons:
            ref_df = abundance_df[abundance_df['experiment'].isin(ref)]
            alt_df = abundance_df[abundance_df['experiment'].isin(alt)]

            ref_df = ref_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control', 'experiment']).agg(
                {'count': 'sum'}).reset_index()
            alt_df = alt_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control', 'experiment']).agg(
                {'count': 'sum'}).reset_index()

            # Collapse the counts if there are more than one being grouped
            if len(ref) > 1:
                ref_df = ref_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control']).agg(
                    {'count': 'sum'}).reset_index()
                ref_df['experiment'] = ref[0]
            if len(alt) > 1:
                alt_df = alt_df.groupby(['position', 'variant', 'residue', 'umi', 'is_control']).agg(
                    {'count': 'sum'}).reset_index()
                alt_df['experiment'] = alt[0]

            jobs.append((alt[0], ref_df.copy(), alt_df.copy(), dt, normalize, N, counts_cutoff))

        if n_processes == 1:
            import itertools
            total_output = itertools.starmap(compute_fitness,
                                             tqdm(jobs, "Computing fitness", total=len(jobs), unit='comparisons'))
        else:
            with Pool(n_processes) as p:
                total_output = p.starmap(compute_fitness,
                                         tqdm(jobs, "Computing fitness", total=len(jobs), unit='comparisons'))

        # for initial_pop, (ref, alt) in tqdm(zip(initial_population, comparisons), "Computing fitness", total=len(comparisons), unit='comparisons'):
        #     boot_results = compute_fitness(ref_df.copy(), alt_df.copy(), dt=dt, n=N)
        #
        #     # Append the results to the total boot output
        #     total_output.append(boot_results)
        total_output = pd.concat(total_output).sort_values(
            ['position', 'variant', 'residue', 'is_control', 'experiment'])
        total_output.to_csv(total_boot_output, index=False)


def run_pipeline_from_toml(config_path: str, overwrite: bool = False):
    """
    Run the pipeline using a TOML configuration file.

    Parameters specified in the TOML file will override defaults.
    Missing parameters will use the run_pipeline defaults.

    :param config_path: Path to the TOML configuration file.
    :param overwrite: Whether to overwrite existing bootstrap results.
    """
    with open(config_path, 'rb') as f:
        config = tomllib.load(f)

    # Call run_pipeline with the config, letting Python handle default values
    run_pipeline(**config, overwrite_boots=overwrite)


def run_pipeline(barcodes_file: str, twist_file: str, experiment_name: str, output_dir: str,
                 data_dirs, reference_names, alt_names, dt=1, paired=False,
                 counts_cutoff=8, normalize=True, barcode_length=16, N=2500,
                 n_processes=15, overwrite_boots=False, capture_seq="CGCTTTAAGGCCGGTCCTAGCAA",
                 make_plots: bool = True):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    experiment_prefix = os.path.join(output_dir, f"{experiment_name}_")
    generate_abundances(experiment_prefix, barcodes_file, twist_file, reference_names, alt_names, data_dirs,
                        capture_seq, barcode_length, n_processes)
    abundance_df = pd.read_csv(experiment_prefix + 'abundance_with_variant.csv')

    # Given the abundance, we can now compute the fitness and bootstraps
    total_boot_output = experiment_prefix + 'boot.csv'
    bootstrap_fitness(abundance_df, total_boot_output, reference_names, alt_names, paired, overwrite_boots, dt,
                      normalize, counts_cutoff, N, n_processes)
    total_output = pd.read_csv(total_boot_output)

    # Make figures
    if make_plots:
        generate_figures(total_output, output_dir, experiment_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the SatSeq analysis pipeline.')
    parser.add_argument('config', help='Path to a TOML configuration file or directory containing TOML files')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing bootstrap results')
    args = parser.parse_args()

    if os.path.isdir(args.config):
        # Run all TOML files in the directory
        toml_files = sorted(glob(os.path.join(args.config, '*.toml')))
        if not toml_files:
            print(f"No TOML files found in directory: {args.config}")
            sys.exit(1)
        print(f"Found {len(toml_files)} TOML files in {args.config}")
        for i, config_path in enumerate(toml_files, 1):
            print(f"\n[{i}/{len(toml_files)}] Running pipeline from config: {config_path}")
            run_pipeline_from_toml(config_path, overwrite=args.overwrite)
        print(f"\nCompleted all {len(toml_files)} experiments.")
    else:
        # Run single TOML file
        print(f"Running pipeline from config: {args.config}")
        run_pipeline_from_toml(args.config, overwrite=args.overwrite)
