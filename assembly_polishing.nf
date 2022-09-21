params.fastq = ""
params.draft_name = ""
params.prot_db = ""

params.output_folder = params.draft_name

include { remove_short_reads; flye_assembly; racon; medaka; proovframe; remove_short_contigs} from './assembly_processes.nf'
include {racon as racon02} from './assembly_processes.nf'

workflow {
    fastq = Channel.fromPath(params.fastq)
    prot_db = Channel.fromPath(params.prot_db)

    remove_short_reads(fastq, params.draft_name)
    flye_assembly(remove_short_reads.out.min500)
    remove_short_contigs(flye_assembly.out.flye)
    racon(remove_short_reads.out.min500, remove_short_contigs.out.min2000, "01")
    racon02(remove_short_reads.out.min500, racon.out.racon, "02")
    medaka(remove_short_reads.out.min500, racon02.out.racon)
    proovframe(medaka.out.medaka, prot_db)
}
