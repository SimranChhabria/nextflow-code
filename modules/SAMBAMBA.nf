process SAMBAMBA_DEDUP {
      publishDir "${params.dedup_dirname}", mode: 'copy' 

     input:
     tuple val(meta), path(bam)

     output:
     tuple val(meta), path('*dd.bam'), emit: dedup_bam
     //tuple val(meta), path('*dd.bam.bai'), emit: dedup_bam_bai
     tuple val(meta), path('*log') , emit: dedup_log

    script:
    def prefix = task.ext.prefix ?: "${meta.seqid}"
    """
    sambamba markdup --remove-duplicates "${bam}" "${prefix}.dd.bam"
    cat .command.err > "${prefix}.log"
    """
  

}

process SAMBAMBA_DEDUP_INDEX {
    publishDir "${params.dedup_dirname}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
 
    output:
    tuple val(meta), path('*dd.bam.bai'), emit: dedup_bam_bai

    script:
    def prefix = task.ext.prefix ?: "${meta.seqid}"

    """
    samtools index "${bam}" > "${prefix}.dd.bam.bai" 
    """
}

process SAMBAMBA_DEDUP_LOG_TABLE {
    publishDir "${params.stats_dirname}", mode: 'copy' 

     input:
     tuple val(meta), path(dedupstat)

     output:
     tuple val(meta), path('*dd.tsv'), emit: dedup_stats

     script:
     def prefix = task.ext.prefix ?: "${meta.seqid}"
     """
     dedup-log2table.R "${dedupstat}" tmp.tsv
     paste_col.py -i tmp.tsv --header "Sample" -v "${prefix}" > "${prefix}.dd.tsv"
    """


}
