#!/usr/bin/env nextflow

params.reads = "./data/reads.fast5"  // Input FAST5 file(s)
params.output = "./results"          // Output directory

process basecalling {
    input:
    path reads from params.reads

    output:
    path "${params.output}/basecalled" into basecalled

    script:
    """
    dorado basecaller --model dna_r9.4.1 ${reads} -o ${params.output}/basecalled
    """
}

process quality_control {
    input:
    path basecalled

    output:
    path "${params.output}/qc" into qc_results

    script:
    """
    fastqc ${basecalled} -o ${params.output}/qc
    pycoQC -i ${basecalled} -o ${params.output}/qc/pycoqc.html
    """
}

process read_filtering {
    input:
    path basecalled

    output:
    path "${params.output}/filtered" into filtered_reads

    script:
    """
    NanoFilt --quality 10 --length 500 < ${basecalled} > ${params.output}/filtered/filtered_reads.fasta
    """
}

process genome_assembly {
    input:
    path filtered_reads

    output:
    path "${params.output}/assembly" into assembly

    script:
    """
    flye --nano-raw ${filtered_reads} --out-dir ${params.output}/assembly
    """
}

process polishing {
    input:
    path assembly

    output:
    path "${params.output}/polished" into polished

    script:
    """
    racon ${filtered_reads} ${assembly}/assembly.fasta > ${params.output}/polished/polished.fasta
    medaka_consensus -i ${params.output}/polished/polished.fasta -d ${assembly}/assembly.fasta -o ${params.output}/polished
    """
}

process variant_calling {
    input:
    path polished

    output:
    path "${params.output}/variants" into variants

    script:
    """
    sniffles -m ${polished}/polished.fasta -v ${params.output}/variants/variants.vcf
    """
}

workflow {
    basecalling | quality_control | read_filtering | genome_assembly | polishing | variant_calling
}
