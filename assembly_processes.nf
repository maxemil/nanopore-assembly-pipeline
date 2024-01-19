params.min_read_length = "500"
params.min_read_quality = "0"
params.min_contig_length = "2000"
params.flye_options = "--nano-raw"
params.medaka_model = "r1041_e82_400bps_sup_v4.3.0"
params.input_format = "fastq.gz"

process remove_short_reads {
  input:
    path 'fastq/*'
    val draft_name
  output:
    path "${draft_name}.l${params.min_read_length}_q${params.min_read_quality}.${params.input_format}", emit: min500

  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    seqkit seq --min-qual ${params.min_read_quality} \
               --min-len ${params.min_read_length} \
               --out-file ${draft_name}.l${params.min_read_length}_q${params.min_read_quality}.${params.input_format} \
               fastq/* 
    """
}

process flye_assembly {
  input:
    file fastq 
  output:
    tuple val("${fastq.simpleName}"), path("${fastq.simpleName}-flye/"), emit: flye

  label "high_cpu"
  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    flye -o ${fastq.simpleName}-flye -t ${task.cpus} --meta ${params.flye_options} ${fastq} &> ${fastq.simpleName}-flye.log
    """
}

process remove_short_contigs {
  input:
    tuple val(pre), path(flye_dir)
  output:
    path "${pre}.flye.fasta", emit: draft
    path "${flye_dir}", emit: flye_dir

  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    flye-post.py ${flye_dir} -m ${params.min_contig_length} -p ${pre}
    cp ${flye_dir}/${pre}.fasta ${pre}.flye.fasta
    """
    // seqkit seq -m 2000 ${assembly} -o ${assembly.simpleName}.min2000.fasta.gz
}

process racon {
  input:
    file fastq
    file assembly
    val iteration
  output:
    path "${assembly.simpleName}.racon${iteration}.fasta", emit: racon

  label "high_cpu"
  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    minimap2 -t ${task.cpus} -x map-ont ${assembly} ${fastq} > ${fastq.simpleName}_vs_${assembly.simpleName}.paf
    racon -m 8 -x -6 -g -8 -w 500 -t ${task.cpus} ${fastq} ${fastq.simpleName}_vs_${assembly.simpleName}.paf ${assembly} > ${assembly.simpleName}.racon${iteration}.fasta
    """
}


process medaka {
  input:
    file fastq
    file assembly
  output:
    path "${assembly.simpleName}.medaka.fasta", emit: medaka
    path "${assembly.simpleName}-medaka/", emit: medaka_dir

  label "high_cpu"
  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    # medaka_consensus -i ${fastq} -d ${assembly} -o ${assembly.simpleName}-medaka -t ${task.cpus} -m ${params.medaka_model}
    run_medaka_parallel.py -a ${assembly} -f ${fastq} -p ${assembly.simpleName}-medaka/${assembly.simpleName} -t ${task.cpus} -m ${params.medaka_model}
    # cp ${assembly.simpleName}-medaka/consensus.fasta ${assembly.simpleName}.medaka.fasta
    cp consensus.fasta ${assembly.simpleName}.medaka.fasta
    """
}

process proovframe {
  input:
    file assembly
    file prot_db
  output:
    path "${assembly.simpleName}.proovframe.fasta", emit: proovframe
    path "*.log", emit: proovframe_logs
    path "${assembly.simpleName}-proovframe.o6", emit: proovframe_map

  label "high_cpu"
  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    proovframe map -a ${prot_db} -t ${task.cpus} -o ${assembly.simpleName}-proovframe.o6 ${assembly} &> ${assembly.simpleName}-proovframe-map.log
    proovframe fix ${assembly} ${assembly.simpleName}-proovframe.o6 -o ${assembly.simpleName}.proovframe.fasta &> ${assembly.simpleName}-proovframe-fix.log
    """
}
