process BAMTOBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir"${params.bedfile}",mode: 'copy'

   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools \\
        bamtobed \\
        -i $bam \\
        > ${prefix}.bed
    
    """

}

process SORT_BED {
    tag "$meta.id"
    publishDir"${params.bedfile}",mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*sorted.bed"), emit: sorted_bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sort_python.py -i ${bam} -o ${prefix}_sorted.bed 
    """


    
}
