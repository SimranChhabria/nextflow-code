nextflow.enable.dsl = 2

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)


 /* -- Include all the modules -- */
include { FASTQC }               from '../modules/FASTQC.nf'
include { FASTP }                from '../modules/FASTP.nf'
include { BOWTIE2_BUILD }        from '../modules/BOWTIE2.nf'
include { BOWTIE2_ALIGN }        from '../modules/BOWTIE2.nf'
include { BAM_SORT_SAMTOOLS }    from '../subworkflow/BAM_SORT_SAMTOOLS.nf'
//include { BWAMEM2_INDEX }        from '../modules/BWA.nf' 
//include { BWAMEM2_MEM }          from '../modules/BWA.nf'
include { BAMTOBED }             from '../modules/BAMTOBED.nf'
include { SORT_BED }             from '../modules/BAMTOBED.nf'

include { MULTIQC }                     from '../modules/MULTIQC.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/custom/CUSTOM_DUMPSOFTWAREVERSIONS.nf'

def multiqc_report    = []
def pass_mapped_reads = [:]
def fail_mapped_reads = [:]

workflow COVID {
    
  take:
  ch_reads_input


  main:
// To gather all QC reports for MultiQC
  ch_reports  = Channel.empty()
// To gather used softwares versions for MultiQC
  ch_versions = Channel.empty()

//QC
  if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
      FASTQC(ch_reads_input)
      ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
      ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

// Trimming and/or splitting
   if (params.trim_fastq || params.split_fastq > 0) {
      save_trimmed_fail = false
      save_merged = false
      FASTP(ch_reads_input,
             [], // we are not using any adapter fastas at the moment
             save_trimmed_fail,
             save_merged)

      ch_reports = ch_reports.mix(
                    FASTP.out.json.collect{meta, json -> json},
                    FASTP.out.html.collect{meta, html -> html}
                  )
      ch_versions = ch_versions.mix(FASTP.out.versions)
      ch_reads_to_map = FASTP.out.reads
      }
  
    // Alignment with bowtie2
    ch_bam                      = Channel.empty()
    ch_bai                      = Channel.empty()
    ch_bowtie2_multiqc          = Channel.empty()
    ch_bowtie2_flagstat_multiqc = Channel.empty()
    
    if (!params.skip_variants) {
      //Indexing wuhan genome
      BOWTIE2_BUILD(file(params.wuhan_ref))
      ch_index    = BOWTIE2_BUILD.out.index
      ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

      //Mapping with wuhan genome BOWTIE2
      sort_bam = false
      BOWTIE2_ALIGN(ch_reads_to_map,ch_index,params.save_unaligned,sort_bam)
      ch_bam                      = BOWTIE2_ALIGN.out.bam
      //ch_bowtie2_multiqc          = BOWTIE2_ALIGN.out.log_out
      ch_reports                  = ch_reports.mix(BOWTIE2_ALIGN.out.log_out.collect{meta, report -> report})
      ch_versions                 = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
      BAM_SORT_SAMTOOLS (
           BOWTIE2_ALIGN.out.bam
       )
      ch_versions                 = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)
      ch_bai                      = BAM_SORT_SAMTOOLS.out.bai
      ch_bowtie2_flagstat_multiqc = BAM_SORT_SAMTOOLS.out.flagstat
      ch_reports                  = ch_reports.mix(BAM_SORT_SAMTOOLS.out.flagstat.collect{meta, report -> report})  

      //CRAM_QC_NO_MD.out.qc.collect{meta, report -> report}
      
      ch_fail_mapping_multiqc = Channel.empty()
      ch_bowtie2_flagstat_multiqc
            .map { meta, flagstat -> [ meta ] + getFlagstatMappedReads(flagstat, params.min_mapped_reads) }
            .set { ch_mapped_reads }
      
      ch_bam
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_bam }

      
      ch_mapped_reads
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }.set { ch_pass_fail_mapped }

      

      
    } 

   ch_version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }
    
     if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        // workflow_summary    = WorkflowPractice.paramsSummaryMultiqc(workflow, summary_params)
        // ch_workflow_summary = Channel.value(workflow_summary)
        

        ch_multiqc_files = Channel.empty()
        ch_reports.view()
        ch_multiqc_files = ch_multiqc_files.mix(ch_version_yaml)
        //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        //ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))
        //ch_multiqc_files.view()

        //ch_multiqc_configs = ch_multiqc_config.mix(ch_multiqc_custom_config).ifEmpty([])

        MULTIQC (
            ch_multiqc_files.collect()
            //ch_multiqc_config.collect().ifEmpty([])
            //ch_multiqc_custom_config.collect().ifEmpty([]),
            //ch_multiqc_logo.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }
}

    


  // // BWA Alignment
  // if (!(params.skip_tools && params.skip_tools.split(',').contains('bwamem2'))) {
  //   //Indexing wuhan genome
  //   BWAMEM2_INDEX(file(params.wuhan_ref))
  //   ch_index    = BWAMEM2_INDEX.out.index
  //   ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
  
  //   //Mapping with wuhan genome BWA MEM
  //    sort_bam = false
  //    BWAMEM2_MEM(ch_reads_to_map,ch_index,sort_bam)
  //    ch_bam_mapped = BWAMEM2_MEM.out.bam
  //    ch_bam_mapped.view()
  //    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)


  // IVAR CONSENSES AND VARIANT CALLING
   

  


//-----------------------
//      FUNCTIONS
// ----------------------


 def getFlagstatMappedReads(flagstat_file, min_mapped_reads) {
        def mapped_reads = 0
        flagstat_file.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }

        def pass = false
        def logname = flagstat_file.getBaseName() - 'flagstat'
        if (mapped_reads > min_mapped_reads.toInteger()) {
            pass = true
        }
        return [ mapped_reads, pass ]
    }





