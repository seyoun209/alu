#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer

# Read and process samplesheet
samples = pd.read_csv(config["update_samplesheet"], sep='\t')
#samples = pd.read_csv("bam_samplesheet.txt",sep='\t')
samples = samples.astype(str)
#samples = samples[~samples['Donor'].str.contains('|'.join(config['samples_to_omit']))]

# Create mn column based on mergeBy config
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy= ['Project','Cell_Type','Condition','Time','Bio_Rep','Tech_Rep','Seq_Rep']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

# Get unique groups based on mergeCondition
group_columns = config['mergeCondition']
#mergeCondition= ['Time']
group_columns = config['mergeCondition']
if isinstance(group_columns, list):
    all_values = samples[group_columns].agg('_'.join, axis=1).unique()
else:
    all_values = samples[group_columns].unique()


##selected_conditions = ['0', '360', '4320']
# Filter to only selected conditions
group_values = [x for x in all_values if x in config['selected_conditions']]
print("Analyzing conditions:", group_values)

# Create sample groups dictionary
sample_groups = {}
for group in group_values:
    if isinstance(group_columns, list):
        mask = samples[group_columns].agg('_'.join, axis=1) == group
    else:
        mask = samples[group_columns] == group
    sample_groups[group] = samples[mask]['mn'].tolist()


rule all:
    input:
        # Sample lists for each condition
        expand("rna_output/samplelists/{condition}.txt",
               condition=group_values),
        # rMATS comparisons
        expand("rna_output/rmats_{ctrl}_vs_{cond}",
               ctrl=config['control_condition'],
               cond=[c for c in group_values if c != config['control_condition']])

rule samplelist:
    output:
        txt = "rna_output/samplelists/{condition}.txt"
    params:
        bam_pattern = "rna_output/align/{sample}_sorted.bam"
    run:
        if wildcards.condition not in group_values:
            raise ValueError(f"Condition {wildcards.condition} not in selected conditions")
        bam_files = [params.bam_pattern.format(sample=sample)
                    for sample in sample_groups[wildcards.condition]]
        with open(output.txt, 'w') as f:
            f.write(','.join(bam_files))

rule rmats:
    input:
        b1 = "splicing_output/samplelists/{control}.txt",
        b2 = "splicing_output/samplelists/{condition}.txt"
    output:
        directory("splicing_output/rmats_{control}_vs_{condition}")
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

