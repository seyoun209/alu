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

# Group by mn and bam
bamfile = samples.groupby('mn')['alignbam'].apply(list).to_dict()

# Build dictionary of merged BAM files using config
mergeSample = samples.groupby(config['mergeCondition'])['mn'].apply(list).to_dict()


## Define actions on success
onsuccess:
        ## Success message
        print("signal track completed successfully! Wahoo!")

##### Define rules #####
rule all:
    input:
        expand("signal_output/signals/{cond}_norm.bw", cond=mergeSample.keys()),
        expand("signal_output/signals/strand/{cond}_{direction}.bw",
               cond=mergeSample.keys(), direction=['fwd', 'rev'])

rule signal:
    input:
        bam = lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        signal = "signal_output/signals/{sampleName}.bw"
    log:
        err = 'signal_output/logs/signal_{sampleName}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p signal_output/signals
        bamCoverage -b {input.bam} -o {output.signal} > {log.err} 2>&1
        """

rule mergeAlign:
    input:
        bam_files = lambda wildcards: ["signal_output/align/{sampleName}_sorted.bam".format(sampleName=value)
                                     for value in mergeSample[wildcards.cond]]
    output:
        bam = "signal_output/merged/{cond}_sorted.bam",
        bai = "signal_output/merged/{cond}_sorted.bam.bai",
        stats = "signal_output/merged/{cond}_stats.txt"
    log:
        err = 'signal_output/logs/mergeAlign_{cond}.err'
    params:
        samtoolsVersion = config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        mkdir -p signal_output/merged
        samtools merge {output.bam} {input.bam_files} 2> {log.err}
        samtools flagstat {output.bam} > {output.stats} 2>> {log.err}
        samtools index {output.bam} 2>> {log.err}
        """

rule mergeSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal = "signal_output/signals/{cond}.bw"
    log:
        err = 'signal_output/logs/mergeSignal_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p signal_output/signals
        bamCoverage -b {input.bam} -o {output.signal} 2> {log.err}
        """

rule mergeForwardSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_fw = "signal_output/signals/strand/{cond}_fwd.bw"
    log:
        err = 'signal_output/logs/mergeSignal_forward_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p signal_output/signals/strand
        bamCoverage --filterRNAstrand forward --normalizeUsing BPM --binSize 10 \
        --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_fw} 2> {log.err}
        """

rule mergeReverseSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_rv = "signal_output/signals/strand/{cond}_rev.bw"
    log:
        err = 'signal_output/logs/mergeSignal_reverse_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p signal_output/signals/strand
        bamCoverage --filterRNAstrand reverse --normalizeUsing BPM --binSize 10 \
        --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_rv} 2> {log.err}
        """

rule mergeSignal_norm:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal = "signal_output/signals/{cond}_norm.bw"
    log:
        err = 'signal_output/logs/mergeSignal_normalized_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p signal_output/signals
        bamCoverage --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 \
        --bam {input.bam} -o {output.signal} 2> {log.err}
        """
