params.fastq = ""
params.draft_name = ""
params.prot_db = ""
params.min_read_length = "500"
params.min_read_quality = "0"
params.min_contig_length = "2000"
params.output_folder = params.draft_name
params.flye_options = "--nano-raw"
params.medaka_model = "r941_min_sup_g507"
params.input_format = "fastq.gz"

include { remove_short_reads; flye_assembly; racon; medaka; proovframe; remove_short_contigs} from './assembly_processes.nf'
include {racon as racon02} from './assembly_processes.nf'

workflow {
    fastq = Channel.fromPath(params.fastq)
    prot_db = Channel.fromPath(params.prot_db)

    remove_short_reads(fastq, params.draft_name)
    flye_assembly(remove_short_reads.out.min500)
    remove_short_contigs(flye_assembly.out.flye)
    racon(remove_short_reads.out.min500, remove_short_contigs.out.draft, "01")
    racon02(remove_short_reads.out.min500, racon.out.racon, "02")
    medaka(remove_short_reads.out.min500, racon02.out.racon)
    proovframe(medaka.out.medaka, prot_db)
}

workflow.onComplete {
    File file = new File("$params.output_folder/${workflow.start}.log")
    file.append("Pipeline $workflow.scriptName started at $workflow.start \n")
    file.append("Pipeline $workflow.scriptName with hash $workflow.scriptId \n")
    file.append("Pipeline $workflow.scriptName was launched at $workflow.launchDir \n")
    file.append("Pipeline $workflow.scriptName was launched as $workflow.commandLine \n")
    file.append("Pipeline $workflow.scriptName completed at $workflow.complete \n")
    file.append("Execution status: ${ workflow.success ? 'OK' : 'failed' } \n")
}