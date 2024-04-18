import ribopy
from ribopy import Ribo
from functions import get_sequence, get_cds_range_lookup, get_psite_offset
import numpy as np
import multiprocessing
import time
import pickle

ribo_path = f'./mouse/all.ribo'
experiments = ['WT_control_A', 'WT_10min_A', 'WT_30min_A', 'WT_1hr_A']
exp = experiments[0]
min_len = 25
max_len = 31
alias = True
ribo_object = Ribo(ribo_path, alias = ribopy.api.alias.apris_human_alias)

def get_filtered_transcripts(ribo_object, experiments, min_len, max_len, alias, threshold):
    region_counts = ribo_object.get_region_counts(region_name = "CDS", range_lower = min_len, range_upper = max_len,
                                                sum_lengths = True, sum_references = False, alias = alias)
    filtered_region_counts = region_counts.loc[:, experiments]

    reads_per_nt_dict = {} 
    for transcript, (start, stop) in get_cds_range_lookup(ribo_object).items():
        cds_length = stop - start
        cds_counts = filtered_region_counts.loc[transcript]
        reads_per_nt = cds_counts.sum() / len(experiments) / cds_length
        reads_per_nt_dict[transcript] = reads_per_nt
        
    filtered_transcripts = {transcript: value for transcript, value in reads_per_nt_dict.items() if value >= threshold}

    return filtered_transcripts


def get_adj_coverage(ribo_object, transcript, exp, min_len, max_len, alias):
    start, stop = get_cds_range_lookup(ribo_object)[transcript]
    offset = get_psite_offset(ribo_object, exp, min_len, max_len)

    if start < max(offset.values()):
        return None
    
    coverages = [
        ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)
        [transcript][start - offset[i] : stop - offset[i]]
        for i in range(min_len, max_len + 1)
    ]

    coverage = sum(coverages, np.zeros_like(coverages[0]))
    return coverage


def process_transcript(transcript):
    coverage = get_adj_coverage(ribo_object, transcript, exp, min_len, max_len, alias)
    return coverage

if __name__ == '__main__':
    start_time = time.time()

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    results = []
    transcripts = get_filtered_transcripts(ribo_object, experiments, min_len, max_len, alias, threshold=1.0)
    for transcript in transcripts:
        result = pool.apply_async(process_transcript, args=(transcript,))
        results.append((transcript, result))  # Store both transcript and its result

    pool.close()
    pool.join()

    adj_coverage_dict = {}
    for transcript, result in results:
        coverage = result.get()  # Get the actual result
        adj_coverage_dict[transcript] = coverage

    with open(f'adj_coverage_filtered_{exp}.pkl', 'wb') as f:
        pickle.dump(adj_coverage_dict, f)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")
