process remove_short_reads {
  input:
    file fastq
    val draft_name
  output:
    path "${draft_name}.min500.fastq.gz", emit: min500

  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    seqkit seq -m 500 ${fastq} -o ${draft_name}.min500.fastq.gz
    """
}

process flye_assembly {
  input:
    file fastq
  output:
    path "${fastq.simpleName}.flye.fasta", emit: flye
    path "${fastq.simpleName}-flye/", emit: flye_dir

  label "high_cpu"
  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    flye -o ${fastq.simpleName}-flye -t ${task.cpus} --meta --nano-raw ${fastq} &> ${fastq.simpleName}-flye.log
    cp ${fastq.simpleName}-flye/assembly.fasta ${fastq.simpleName}.flye.fasta
    """
}

process remove_short_contigs {
  input:
    file assembly
  output:
    path "${assembly.simpleName}.min2000.fasta.gz", emit: min2000

  publishDir "${params.output_folder}", mode: 'copy'

  script:
    """
    seqkit seq -m 2000 ${assembly} -o ${assembly.simpleName}.min2000.fasta.gz
    """
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
    medaka_consensus -i ${fastq} -d ${assembly} -o ${assembly.simpleName}-medaka -t ${task.cpus} -m r941_min_hac_g507
    cp ${assembly.simpleName}-medaka/consensus.fasta ${assembly.simpleName}.medaka.fasta
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
