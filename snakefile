#! /usr/bin/env bash
## Snakefile
####################
Samples = ['1','2','3','4','5','6','7','8','9','10','11','12','NCM','NCR']
Numbers = ['0']
BOWTIE2_INDEX = "/storage/zhangyanxiaoLab/share/bowtie2_index/hg38"
MARKDUP="/storage/zhangyanxiaoLab/share/Pipelines/atac-seq-pipeline-snakemake/dependencies/picard.jar MarkDuplicates"

rule all:
  input: 
    expand("bam/NPCD{number}_{sample}_1.nodup.bam",number=Numbers,sample=Samples),
    expand("bam/NPCD{number}_{sample}_1.nodup.bam.bai",number=Numbers,sample=Samples),
    expand("bigWig/NPCD{number}_{sample}_1.nodup.bw",number=Numbers,sample=Samples),

    expand("bam/NPCD{number}_{sample}_2.nodup.bam",number=Numbers,sample=Samples),
    expand("bam/NPCD{number}_{sample}_2.nodup.bam.bai",number=Numbers,sample=Samples),
    expand("bigWig/NPCD{number}_{sample}_2.nodup.bw",number=Numbers,sample=Samples),

    expand("bam/NPCD{number}_{sample}_merge.nodup.bam",number=Numbers,sample=Samples),
    expand("bam/NPCD{number}_{sample}_merge.nodup.bam.bai",number=Numbers,sample=Samples),
    expand("bigWig/NPCD{number}_{sample}_merge.nodup.bw",number=Numbers,sample=Samples),

    expand("peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_peaks.narrowPeak",number=Numbers,sample=Samples),
    expand("peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_peaks.xls",number=Numbers,sample=Samples),
    expand("peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_summits.bed",number=Numbers,sample=Samples)

rule bowtie2_align1:
  output: 
    bam=temp("bam/NPCD{number}_{sample}_1.sorted.bam"),
    raw_qc = "qc/NPCD{number}_{sample}_1.raw.flagstat.qc",
    log="logs/NPCD{number}_{sample}_1.bowtie2.log"
  input:
#    lambda wildcards: FASTQ_DICT[wildcards.sample]
    "D{number}_{sample}_1/D{number}_{sample}_1_1.clean.fq.gz",
    "D{number}_{sample}_1/D{number}_{sample}_1_2.clean.fq.gz"
  threads: 10 
  run:
    print(input)
    if len(input) == 1:
      middle = "-U " + str(input)
    elif len(input) == 2:
#     middle = "-1 " + str(input[0]).replace("fastq","fastq_trim",1).replace(".fastq","_val_1.fq") + " -2 " + str(input[1]).replace("fastq","fastq_trim",1).replace(".fastq","_val_2.fq") + " -X 2000"
      middle = "-1 " + str(input[0]) + " -2 " + str(input[1])  + " -X 2000"
    shell(
    "bowtie2 -x {BOWTIE2_INDEX} "
    "{middle} "
    "-p {threads} 2> logs/NPCD{wildcards.number}_{wildcards.sample}_1.bowtie2.log|"
    "samtools view -bS |"
    "samtools sort -@ {threads} -m 4G > {output.bam};"
    "samtools flagstat {output.bam} > {output.raw_qc};"
    )

rule bowtie2_align2:
  output: 
    bam=temp("bam/NPCD{number}_{sample}_2.sorted.bam"),
    raw_qc = "qc/NPCD{number}_{sample}_2.raw.flagstat.qc",
    log="logs/NPCD{number}_{sample}_2.bowtie2.log"
  input:
#    lambda wildcards: FASTQ_DICT[wildcards.sample]
    "D{number}_{sample}_2/D{number}_{sample}_2_1.clean.fq.gz",
    "D{number}_{sample}_2/D{number}_{sample}_2_2.clean.fq.gz"
  threads: 10 
  run:
    print(input)
    if len(input) == 1:
      middle = "-U " + str(input)
    elif len(input) == 2:
#     middle = "-1 " + str(input[0]).replace("fastq","fastq_trim",1).replace(".fastq","_val_1.fq") + " -2 " + str(input[1]).replace("fastq","fastq_trim",1).replace(".fastq","_val_2.fq") + " -X 2000"
      middle = "-1 " + str(input[0]) + " -2 " + str(input[1])  + " -X 2000"
    shell(
    "bowtie2 -x {BOWTIE2_INDEX} "
    "{middle} "
    "-p {threads} 2> logs/NPCD{wildcards.number}_{wildcards.sample}_2.bowtie2.log|"
    "samtools view -bS |"
    "samtools sort -@ {threads} -m 4G > {output.bam};"
    "samtools flagstat {output.bam} > {output.raw_qc};"
    )

rule bam_rmdup1:
  input:
    bam = "bam/NPCD{number}_{sample}_1.sorted.bam",
  output:
    bam = "bam/NPCD{number}_{sample}_1.nodup.bam",
    bai = "bam/NPCD{number}_{sample}_1.nodup.bam.bai",
    qc = "qc/NPCD{number}_{sample}_1.dup.qc"
  log:
    "logs/markdup/NPCD{number}_{sample}_1.markdup.log"
  threads: 3
  shell:
    "java -Xmx12G -XX:ParallelGCThreads=3 -jar {MARKDUP} TMP_DIR=tmp/NPCD{wildcards.number}_{wildcards.sample}_1 INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log};"
    "samtools index {output.bam}"

rule bam_rmdup2:
  input:
    bam = "bam/NPCD{number}_{sample}_2.sorted.bam",
  output:
    bam = "bam/NPCD{number}_{sample}_2.nodup.bam",
    bai = "bam/NPCD{number}_{sample}_2.nodup.bam.bai",
    qc = "qc/NPCD{number}_{sample}_2.dup.qc"
  log:
    "logs/markdup/NPCD{number}_{sample}_2.markdup.log"
  threads: 3
  shell:
    "java -Xmx12G -XX:ParallelGCThreads=3 -jar {MARKDUP} TMP_DIR=tmp/NPCD{wildcards.number}_{wildcards.sample}_2 INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log};"
    "samtools index {output.bam}"

rule bam2bigwig1:
  input:
    bam = "bam/NPCD{number}_{sample}_1.nodup.bam",
  output: 
    bw = "bigWig/NPCD{number}_{sample}_1.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule bam2bigwig2:
  input:
    bam = "bam/NPCD{number}_{sample}_2.nodup.bam"
  output: 
    bw = "bigWig/NPCD{number}_{sample}_2.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule merge:
    input:
        bam1 = "bam/NPCD{number}_{sample}_1.nodup.bam",
        bam2 = "bam/NPCD{number}_{sample}_2.nodup.bam"
    output:
        bam = "bam/NPCD{number}_{sample}_merge.nodup.bam",
        bai = "bam/NPCD{number}_{sample}_merge.nodup.bam.bai"
    shell:
        "samtools merge {output.bam} {input.bam1} {input.bam2}; samtools index {output.bam}"

rule bam2bigwig_merge:
  input:
    bam = "bam/NPCD{number}_{sample}_merge.nodup.bam"
  output: 
    bw = "bigWig/NPCD{number}_{sample}_merge.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule macs2callpeaks:
    input:
        "bam/NPCD{number}_{sample}_merge.nodup.bam"
    output:
        "peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_peaks.narrowPeak",
        "peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_peaks.xls",
        "peaks/NPCD{number}_{sample}_merge/NPCD{number}_{sample}_summits.bed"
    log:
        "logs/NPCD{number}_{sample}_merge_peaks.log"
    shell:
        "macs2 callpeak -t {input} -f BAMPE -n NPCD{wildcards.number}_{wildcards.sample} --outdir peaks/NPCD{wildcards.number}_{wildcards.sample}_merge -g hs --keep-dup all"
        
