/*
 * Pipeline default parameters
 */
 params.reference = "/ceph/cbio/users/lindo/REFERENCE/ref2.fasta"
 params.dict = "/ceph/cbio/users/lindo/REFERENCE/ref2.dict"
 params.refindex = "/ceph/cbio/users/lindo/REFERENCE/ref2.fasta.fai"
 params.bam = "$baseDir/*.bam"
 params.gatk = "/usr/local/pipeline/Tools/gatk-4.1.0.0/gatk"
 params.outdir = "$baseDir"

/*
 * Parse the input parameters
 */
 genome_file = file(params.reference)
 dictionary = file(params.dict)
 reference_index = file(params.refindex)
 bam_files = Channel.fromPath(params.bam)
 gatk = params.gatk

 println """\
          V A R I A N T - C A L L I N G    P I P E L I N E    W I T H    G A T K
          ===================================
          GENOME        : ${params.reference}
          GENOME INDEX  : ${params.refindex}
          DICTIONARY    : ${params.dict}
          BAM FILES     : ${bam_files}
          OUTDIR        : ${params.outdir}
          """
          .stripIndent()

/**********
 * PART 1: Variant Calling
 */

process '1A_variant_calling' {
  input:
      file reference from genome_file
      file dict from dictionary
      file refindex from reference_index
      file bam from bam_files

  output:
      file "GATK_${bam.baseName}.vcf"

  script:
  """
    $gatk HaplotypeCaller \
        -R $reference \
        -I $bam \
        -O GATK_${bam.baseName}.vcf
  """
}

workflow.onComplete = {
   // any workflow property can be used here
   println "Pipeline complete"
   println "Command line: $workflow.commandLine"
}
