/*
 * Pipeline default parameters
 */
 params.reference = "$baseDir/reference.fa"
 params.refindex = "$baseDir/reference.fa.fai"
 params.dict = "$baseDir/reference.dict"
 params.bam = "$baseDir/*.bam"
 params.bamindex = "$baseDir/*.bai"
 params.gatk = "/usr/local/pipeline/Tools/gatk-4.1.0.0/gatk"
 params.outdir = "$baseDir"

/*
 * Parse the input parameters
 */
 genome_file = file(params.reference)
 dictionary = file(params.dict)
 reference_index = file(params.refindex)
 bam_files = Channel.fromPath(params.bam)
 bam_index = Channel.fromPath(params.bamindex)
 GATK = params.gatk

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
  publishDir '/create/a/new/directory/GATK'

  input:
      file reference from genome_file
      file dict from dictionary
      file refindex from reference_index
      file bam from bam_files
      file bamindex from bam_index

  output:
      file "${bam.baseName}_gatk.vcf" into VCF
      file "${bam.baseName}.g.vcf" into GVCF

  script:
  """
    $GATK HaplotypeCaller \
        -R $reference \
        -I $bam \
        -O ${bam.baseName}_gatk.vcf &
    $GATK HaplotypeCaller \
        -R $reference \
        -I $bam \
        -O ${bam.baseName}.g.vcf \
        -ERC GVCF
  """
}

workflow.onComplete = {
   // any workflow property can be used here
   println "Pipeline complete"
   println "Command line: $workflow.commandLine"
}
