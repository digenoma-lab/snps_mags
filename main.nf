#!/usr/bin/env nextflow




process FASTQC {
    tag "$sample"
    //label 'process_medium'
    publishDir "$params.outdir/QC/FASTQC", mode: "copy"

    container "oras://community.wave.seqera.io/library/fastqc_mosdepth:14f170d6650930e0"

    input:
    tuple val(sample), path(reads)

    output:
    path("${sample}.fastqc"), emit: fqc 

    script:
    if(params.debug == true){
    """
    echo fastqc -o ${sample}.fastqc ${reads[0]} ${reads[1]}
    mkdir -p ${sample}.fastqc
    touch ${sample}.fastqc/report.fastqc
    
    """
    } else{
    """
    mkdir -p ${sample}.fastqc
    fastqc -t $task.cpus -o ${sample}.fastqc ${reads[0]} ${reads[1]}
    """
    }
    
}

process bwa_index {
    tag "Index Ref BWA"
    publishDir "${params.outdir}/bwa_index", mode: 'copy'
    
    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"


    input:
    path reference

    output:
    path "${reference}.*", emit: index

    script:
    if(params.debug){
    """
    echo bwa-mem2 index ${reference}
    touch ${reference}.ref
    touch ${reference}.bwt
    touch ${reference}.sai
    """
   }else{
    """
    bwa-mem2 index ${reference}
    """
   }
}

process samtools_index {
    tag "Index ref Sam"
    publishDir "${params.outdir}/samtools_index", mode: 'copy'


    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai

    script:
    if(params.debug){
 	"""
	echo samtools faidx ${reference}
        touch ${reference}.fai
 	"""
    }else{
    """
    samtools faidx ${reference}
    """
    }
}



process aln_pipe{

    tag "aln-${sample}"
    publishDir "${params.outdir}/alignment_results", mode: 'copy'
    
    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    tuple val(sample), path(reads)
    path index
    path reference

    output:
    tuple val(sample), path("${sample}.marked.bam"), emit: marked_bam
    path "${sample}.marked.bam.bai" , emit: bai
    path "${sample}.log.bwamem" , emit: log

    script:
    if(params.debug){
	"""
         echo "bwa-mem2 mem -R '@RG\\tID:${sample}\\tSM:${sample}\\tLB:${params.rlibrary}\\tPL:${params.rplat}' -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} 2> ${sample}.log.bwamem"
    echo "samtools sort -@ ${task.cpus} - |"
    echo "samtools markdup -@ ${task.cpus} - ${sample}.marked.bam"
    
    # Index final BAM
    echo samtools index ${sample}.marked.bam

	touch ${sample}.marked.bam
        touch ${sample}.marked.bam.bai
        touch ${sample}.log.bwamem
        """
   }else{
	"""
   bwa-mem2 mem -R '@RG\\tID:${sample}\\tSM:${sample}\\tLB:${params.rlibrary}\\tPL:${params.rplat}' -t ${task.cpus} ${index[0]} ${reads[0]} ${reads[1]} 2> ${sample}.log.bwamem | \
    samtools sort -@ ${task.cpus} - | \
    samtools markdup -@ ${task.cpus} - ${sample}.marked.bam
    
    # Index final BAM
    samtools index ${sample}.marked.bam
 
	"""
   }
}


process qualimap {
    tag "${sample}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'


    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    tuple val(sample), path(bam)
    path reference

    output:
    path "${sample}", emit: qualimap_results

    script:
    if(params.debug){
    """
     echo qualimap bamqc \
        -bam ${bam} \
        -outdir qualimap_results/${sample} \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
         mkdir ${sample}
    """
    }else{
    """
    qualimap bamqc \
        -bam ${bam} \
        -outdir ${sample} \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
    """
    }
}

process instrain_variant_calling {
    tag "${sample}"
    publishDir "${params.outdir}/instrain", pattern: "*.{html,stb,tsv}", mode: 'copy'


    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    tuple val(sample), path(marked_bam)
    path reference
    path fai

    output:
    path "${sample}_instrain", emit: instrain_out

    script:
    if(params.debug){
     """
    echo inStrain profile \
        ${marked_bam} \
        ${reference} \
        -o ${sample}_instrain \
        -p ${task.cpus}
       mkdir ${sample}_instrain
    """
  
    }else{
    """
    inStrain profile \
        ${marked_bam} \
        ${reference} \
        -o ${sample}_instrain \
        -p ${task.cpus} 
    """
    }
}


process DEPTH{
    tag "$sample-depth"
    //label 'process_medium'
    
    publishDir "$params.outdir/DEPTH", mode: "copy"

    container "oras://community.wave.seqera.io/library/fastqc_mosdepth:14f170d6650930e0"
    input:
    tuple val(sample), path(marked_bam)
    path reference
    

    output:
    path("${sample}.depth.*"), emit: depth
    
    script:
    if(params.debug == true){
    """
    echo mosdepth -f ${reference} -t $task.cpus ${sample}.depth $marked_bam
    touch ${sample}.depth.mosdepth.dist.txt
    touch ${sample}.depth.mosdepth.summary.txt
    """
    }else{
    """
    mosdepth -f ${reference} -t $task.cpus ${sample}.depth $marked_bam
    """
    }   
}


process multiqc {
    tag "Generating MultiQC report"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    path "*"

    output:
    path "multiqc*"

    script:
    if(params.debug){
       """
	echo multiqc . -f -o . --interactive
        touch multiqc_report.html
        
       """
    }else{
    """
    multiqc . -f -o . --interactive
    """
    }
}

workflow {

    reference_file = file(params.reference)
    read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)

    bwa_index(reference_file)
    samtools_index(reference_file)
   
   FASTQC(read_pairs) 
   aln_pipe(read_pairs, bwa_index.out.index, reference_file)

    //alignment_pipeline(read_pairs, bwa_index.out.index)

    qualimap(aln_pipe.out.marked_bam, reference_file)
    instrain_variant_calling(aln_pipe.out.marked_bam, reference_file, samtools_index.out.fai)
    DEPTH(aln_pipe.out.marked_bam, reference_file)
   //we create a report of the alignments and mapping
   inmul=qualimap.out.qualimap_results.collect().mix(instrain_variant_calling.out.instrain_out.collect()).flatten().collect()
   inmul=inmul.mix(DEPTH.out.depth.collect()).flatten().collect()
   // inmul.view()

    multiqc(inmul)

}
