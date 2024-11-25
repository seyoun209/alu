#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer

# Read and process samplesheet
samples = pd.read_csv(config["update_samplesheet"], sep='\t')
samples = samples.astype(str)
samples = samples[~samples['Donor'].str.contains('|'.join(config['samples_to_omit']))]
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

# Group samples by condition and timepoint
conditions = samples[config['condition_column']].unique()
timepoints = config['timepoints']
control_condition = config['control_condition']

# Create sample lists by condition and timepoint
sample_groups = {}
for condition in conditions:
    for timepoint in timepoints:
        key = f"{condition}_{timepoint}"
        sample_groups[key] = samples[
            (samples[config['condition_column']] == condition) & 
            (samples[config['time_column']] == timepoint)
        ]['mn'].tolist()


rule all:
    input:

rule samplelist:
    output:
        txt = "splicing_output/samplelists/{condition}_{timepoint}.txt"
    params:
        bam_pattern = "rna_output/align/{sample}_sorted.bam"
    run:
        group_key = f"{wildcards.condition}_{wildcards.timepoint}"
        bam_files = [params.bam_pattern.format(sample=sample) 
                    for sample in sample_groups[group_key]]
        with open(output.txt, 'w') as f:
            f.write(','.join(bam_files))

rule rmats:
    input:
        b1 = "splicing_output/samplelists/{control_timepoint}.txt",
        b2 = "splicing_output/samplelists/{condition}_{timepoint}.txt"
    output:
        directory("splicing_output/rmats_{condition}")
    params:
        rmats_ver = config['rmats_turbo'],
        gtf = config['gtf'],
        strand = config['t'],
        readlen = config['readlength']
    threads: 4
    shell:
        """
        module load rmats-turbo/{params.rmats_ver}
        run_rmats --b1 {input.b1} --b2 {input.b2} \
            --gtf {params.gtf} -t {params.strand} \
            --readLength {params.readlen} \
            --nthread {threads} \
            --od {output} \
            --tmp {output}/tmp
        """
