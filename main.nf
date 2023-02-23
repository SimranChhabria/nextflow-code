#!/usr/bin/env nextflow


/* This is the main nextflow script which calls all the subworkflows */

nextflow.enable.dsl = 2

/* Include all the workflows */

include {  ILLUMINA }   from './workflows/ILLUMINA.nf'


workflow {
    
     ILLUMINA()
}