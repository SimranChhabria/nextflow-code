#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Check input path parameters to see if they exist
def checkPathParamList = [
     params.samplesheetBCL, params.samplesheetmeta
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


if (params.samplesheetBCL)        { ch_samplesheetBCL      = file(params.samplesheetBCL)      } else { exit 1, 'BCL samplesheet file not specified!' }
if (params.samplesheetmeta)       { ch_samplesheetmeta     = file(params.samplesheetmeta)     } else { exit 1, 'Meta samplesheet file not specified!' } 

include { DEMULTIPLEX } from './DEMULTIPLEX.nf'
include { COVID }       from './COVID.nf'

workflow ILLUMINA {
     
    //Download the raw data from basespcae and demultiplexing 
    //DEMULTIPLEX()
    //ch_input = DEMULTIPLEX.out.ch_fastq

    //Input the reads from the folder
    ch_input = Channel.fromFilePairs(params.inputfastq)
    
    //Combine the meta data and fastq reads into one channel
    ch_input_sample = join_csv((file(params.samplesheetmeta)), ch_input)
    

    //Subset the fastq files based on project name:
     

    //---- COVID ANALYSIS -----

    COVID(ch_input_sample)
    //ONCOLOGY()     

    

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def join_csv(csv_file, input_fastq) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
           def seqID  = row['SEQ_ID']
           def client = row['Client']
           def projID = row['Sample_Project']
           def clientID =row['Client_Sample_ID']
           def ref   = row['Ref']
           def primer_set = row['Primer_Set']
           def output_format = row['Output_Format']
           return [seqID,projID,client,clientID,ref,primer_set,output_format]
        }.join(input_fastq)
         .map { seqID,proID,client,clientID,ref,primer_set,output_format,fastq ->
          def meta = [
                    "id"      : seqID,
                    "projectID"  : proID,
                    "Client"     : client.toString(),
                    "ClientID"   : clientID.toString(),
                    "Ref"        : ref.toString(),
                    "primer"     : primer_set.toString(),
                    "output"     : output_format.toString(),
                    ] 
                    return[meta, fastq]
                    }
         

}
