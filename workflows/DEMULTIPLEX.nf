#!/usr/bin/env nextflow

//
// Download raw data from basespace and demultiplex Illumina BCL data using bcl-convert 
//

include { BSDOWNLOAD }     from '../modules/BCLCONVERT.nf'
include { BCL_CONVERT }     from '../modules/BCLCONVERT.nf'



workflow DEMULTIPLEX {

       main:
       ch_versions = Channel.empty()
       ch_fastq    = Channel.empty()
       ch_reports  = Channel.empty()
       ch_stats    = Channel.empty()
          
    /*  -- * Downloading raw data (bs download) using run.id -- */

        //BS_DOWNLOAD  (params.runid)

    /* --Collects the output from BS_DOWNLOAD for BCL CONVERT (Demultiplexing) -- */
                     
       //BCL_CONVERT  (params.samplesheetBCL, BS_DOWNLOAD.out.collect())
        BCL_CONVERT  (params.samplesheetBCL,params.raw_dir)
        ch_fastq    = ch_fastq.mix(BCLCONVERT.out.fastq)
        ch_interop  = ch_interop.mix(BCLCONVERT.out.interop)
        ch_reports  = ch_reports.mix(BCLCONVERT.out.reports)
        ch_versions = ch_versions.mix(BCLCONVERT.out.versions)

        // Generate meta for each fastq
        ch_fastq_with_meta = generate_fastq_meta(ch_fastq)

        emit:
        fastq    = ch_fastq_with_meta
        reports  = ch_reports
        stats    = ch_stats
        interop  = ch_interop
        versions = ch_version

        


}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def generate_fastq_meta(ch_reads){


            ch_reads 
            .flatten() // Seperate all the fastq files
            .map { file -> 
            //project =   file.toString().tokenize('/').get(8);
            projectID = file.getParent()
            project = projectID.getSimpleName().toString() //Get the ProjectID
            file =  file;
            return [project, file] } //Return ProjectID and fastq
}

