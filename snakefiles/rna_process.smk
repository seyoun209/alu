#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakemake.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
#samples = pd.read_csv("rna_samplesheet.txt",sep='\t')

## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy= ['Project','Cell_Type','Condition','Time','Bio_Rep','Tech_Rep','Seq_Rep']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

## Group by mn and extract Read1 & Read2
read1 = samples.groupby('mn')['Read1'].apply(list).to_dict()
read2 = samples.groupby('mn')['Read2'].apply(list).to_dict()

## Set run summary name using helper script
runName = namer(samples, config['mergeBy'])
#runName = namer(samples, samples[mergeBy])

#align output format

## Define actions on success
onsuccess:
	## Success message
	print("RNA_Splicing completed successfully! Wahoo!")


##### Define rules #####
rule all:
    input:
#        [expand("rna_output/fastq/{sampleName}_{read}.fastq.gz",
#                sampleName=key, read=['R1','R2'])
#         for key in read1],

        # FastQC outputs
        [expand("rna_output/QC/{sampleName}_{read}_fastqc.{ext}",
                sampleName=key, read=['R1', 'R2'], ext=['zip', 'html'])
         for key in read1],

        # Trimming outputs
#        [expand('rna_output/trim/{sampleName}_{read}{ext}',
#                sampleName=key,
#                read=['R1', 'R2'],
#                ext=['_val_1.fq.gz', '_val_2.fq.gz', '.fastq.gz_trimming_report.txt'])
#         for key in read1],

        # Alignment outputs
#        [expand("rna_output/align/{sampleName}_{ext}",
#                sampleName=key,
#                ext=["Aligned.sortedByCoord.out.bam",
#                     "Log.out",
#                     "Log.progress.out",
#                     "Log.final.out"])
#         for key in read1],

        # BAM indexes
#        [expand('rna_output/align/{sampleName}_sorted.{ext}',
#                sampleName=key,
#                ext=['bam', 'bai'])
#         for key in read1],

        # Quantification
        [expand('rna_output/quant/{sampleName}',
                sampleName=key)
         for key in read1]


rule catReads:
    input:
        R1 = lambda wildcards: read1.get(wildcards.sampleName),
        R2 = lambda wildcards: read2.get(wildcards.sampleName)
    output:
        R1 = 'rna_output/fastq/{sampleName}_R1.fastq.gz',
        R2 = 'rna_output/fastq/{sampleName}_R2.fastq.gz'
    benchmark:
        'rna_output/benchmarks/{sampleName}_catReads.tsv'
    log:
        errR1 = 'rna_output/logs/{sampleName}_R1_catReads.err',
        errR2 = 'rna_output/logs/{sampleName}_R2_catReads.err'
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.errR1}
        cat {input.R2} > {output.R2} 2> {log.errR2}
        """

rule fastqc:
    input:
        QC1 = rules.catReads.output.R1,
        QC2 = rules.catReads.output.R2
    output:
        zip1 = "rna_output/QC/{sampleName}_R1_fastqc.zip",
        zip2 = "rna_output/QC/{sampleName}_R2_fastqc.zip",
        html1 = "rna_output/QC/{sampleName}_R1_fastqc.html",
        html2 = "rna_output/QC/{sampleName}_R2_fastqc.html"
    log:
        err = 'rna_output/logs/fastqc_{sampleName}.err',
        out = 'rna_output/logs/fastqc_{sampleName}.out'
    params:
        dir = "rna_output/QC",
        version = config['fastqcVers']
    benchmark:
        'rna_output/benchmarks/fastqc_{sampleName}.tsv'
    shell:
        """
        module load fastqc/{params.version};
        fastqc -o {params.dir} {input.QC1} {input.QC2} 1> {log.out} 2> {log.err};
        """

rule trim:
    input:
        R1 = rules.catReads.output.R1,
        R2 = rules.catReads.output.R2
    output:
        trim1 = 'rna_output/trim/{sampleName}_R1_val_1.fq.gz',
        trim2 = 'rna_output/trim/{sampleName}_R2_val_2.fq.gz',
        report1 = 'rna_output/trim/{sampleName}_R1.fastq.gz_trimming_report.txt',
        report2 = 'rna_output/trim/{sampleName}_R2.fastq.gz_trimming_report.txt'
    threads: 4
    params:
        trim_galore_ver = config['trim_galore'],
        python_ver = config['pythonVers'],
        pigz_ver = config['pgizVers']
    log:
        err = 'rna_output/logs/trim_{sampleName}.err',
        out = 'rna_output/logs/trim_{sampleName}.out'
    shell:
        """
        module load trim_galore/{params.trim_galore_ver}
        module load python/{params.python_ver}
        module load pigz/{params.pigz_ver}
        mkdir -p rna_output/trim
        trim_galore -o rna_output/trim --cores {threads} --paired {input.R1} {input.R2} 2> {log.err}
        """

rule align:
    input:
        R1 = rules.trim.output.trim1,
        R2 = rules.trim.output.trim2
    output:
        bam = "rna_output/align/{sampleName}_Aligned.sortedByCoord.out.bam",
        log_out = "rna_output/align/{sampleName}_Log.out",
        progress_log = "rna_output/align/{sampleName}_Log.progress.out",
        final_log = "rna_output/align/{sampleName}_Log.final.out"
    threads: 8
    log:
        err = 'rna_output/logs/align_{sampleName}.err',
        out = 'rna_output/logs/align_{sampleName}.out'
    params:
        index = config['star'],
        sjdb = config['starsjdb'],
        starVer=config['starVers'],
        dir_align = "rna_output/align/{sampleName}_"
    shell:
        """
        ml star/{params.starVer};
        mkdir -p rna_output/align
        
        STAR --genomeDir {params.index} \
                --runThreadN {threads} \
                --sjdbFileChrStartEnd {params.sjdb} \
                --outFileNamePrefix {params.dir_align} \
                --outSAMstrandField intronMotif \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesIn {input.R1} {input.R2} \
                --outFilterType BySJout \
                --outFilterMultimapNmax 20 \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 1> {log.out} 2> {log.err}
        """

rule index:
    input:
        rules.align.output.bam
    output:
        temp='rna_output/align/{sampleName}_Aligned.sortedByCoord.out.bam.bai',
        bam= 'rna_output/align/{sampleName}_sorted.bam',
        bai = 'rna_output/align/{sampleName}_sorted.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        out = 'rna_output/logs/index_{sampleName}.out',
        err = 'rna_output/logs/index_{sampleName}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools index -@ {threads} {input} {output.temp} 
        samtools view -b {input} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY |samtools sort -o {output.bam} 
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err}
        """

rule quant:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        'rna_output/quant/{sampleName}'
    params:
        salmonVer = config['salmonVers'],
        index=config['salmon'],
        gcFlag = config['gcBias'],
        seqFlag = config['seqBias']
    log:
        out='rna_output/logs/salmon_{sampleName}.out',
        err='rna_output/logs/salmon_{sampleName}.err'
    shell:
        """
        ml salmon/{params.salmonVer}
        mkdir -p rna_output/quant
        salmon quant -i {params.index} -l A -1 {input.trim1} -2 {input.trim2} -o {output} --seqBias --gcBias --validateMappings
        """
