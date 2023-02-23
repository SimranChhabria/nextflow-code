process SAMTOOLS_FLAGSTATS {
    //tag "${meta.id}"
    label 'samtools_flagstat'
    publishDir "${params.stats_dirname}", mode: 'copy' , pattern : '*tsv'
   
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*tsv') , emit: flagstats

    script:
    def prefix = task.ext.prefix ?: "${meta.seqid}"
 
    """
    samtools flagstat ${bam} > ${prefix}_flagstats.tsv
    """

}

process SAMTOOLS_FLAGSTATS_TABLE {
     publishDir "${params.stats_dirname}", mode: 'copy' 

     input:
     tuple val(meta), path(flagstat)

     output:
     tuple val(meta), path('*sampleID.tsv'), emit: flagstat_sampleID

     script:
     def prefix = task.ext.prefix ?: "${meta.seqid}"
     """
     flagstat2table.R "${flagstat}" tmp.tsv
     paste_col.py -i tmp.tsv --header "Sample" -v "${prefix}" > "${prefix}_flagstat_sampleID.tsv"
     """
}
   
    
