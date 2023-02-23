process BSDOWNLOAD {
    label 'process_medium'

    
    input:
    val (run_id) //Takes the run_id (params.run_id)

    output:
    path "${params.run_name}_raw"  
          

    script:
    """
    mkdir -p ${params.run_name}_raw
    bs download run -i ${run_id} -o ${params.run_name}_raw --extension=bcl.gz,bcl.bgzf,cbcl,bci,stats,locs,filter,bin,xml
    """


}


process BCL_CONVERT {
    label 'process_high'
    container "nfcore/bclconvert:4.0.3"
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BCLCONVERT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path (samplesheet)
    path (run_dir)

    output:
    path("**[!Undetermined]_S*_R?_00?.fastq.gz")   ,emit: fastq
    path("**[!Undetermined]_S*_I?_00?.fastq.gz")   ,optional:true ,emit: fastq_idx
    path("**Undetermined_S0*_R?_00?.fastq.gz")     ,optional:true ,emit: undetermined
    path("**Undetermined_S0*_I?_00?.fastq.gz")     ,optional:true, emit: undetermined_idx
    path("Reports/Demultiplex_Stats.csv")          ,emit: demultiplex
    path("Reports/*Top_Unknown_Barcodes.csv")      ,emit: unknown_barcodes
    path("Reports")                                ,emit: reports
    path("Logs")                                   ,emit: logs
    path("**/InterOp/*.bin")                       ,emit: interop
    path("versions.yml")                                            ,emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''

    """
    bcl-convert \\
        $args \\
        --output-directory . \\
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --no-lane-splitting=true \\
         --force \\
        --bcl-sampleproject-subdirectories=true\\
        --bcl-num-parallel-tiles ${task.cpus}

    mv Reports/Top_Unknown_Barcodes.csv Reports/${params.run_name}_Top_Unknown_Barcodes.csv
    scp Reports/${params.run_name}_Top_Unknown_Barcodes.csv ${params.QC_dirname}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """
}
