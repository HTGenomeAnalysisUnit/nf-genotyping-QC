nextflow.enable.dsl=2

outdir = params.outdir ? params.outdir : "QC"
file(outdir).mkdirs()
outdir_samples = "${outdir}/qc_samples"
outdir_vars = "${outdir}/qc_variants"

workflow {
    log.info """\
    =============================================
      GENOTYPE QC    -    N F   P I P E L I N E    
    =============================================
    input file   : ${params.input}
    input format : ${params.format}
    =============================================
    """
    .stripIndent()

    if (params.format == "bed") {
        genotype_files = Channel
            .fromFilePairs("*.{bed,bim,fam}", size: 3)
    } else if (params.format == "vcf") {
        genotype_files = Channel
            .fromFilePath("*.vcf.gz")
            .map { tuple(it.name, it) }
    }

    PLINK_SAMPLE_MISS(genotype_files)
    PLINK_SAMPLE_HET(genotype_files)
    PLINK_SEX_CHECK(genotype_files)
    PLINK_IBD(PLINK_PRUNE_VARS.out)

    metrics_ch = PLINK_SAMPLE_MISS.out.mix(
        PLINK_SAMPLE_HET.out,
        PLINK_SEX_CHECK.out,
        PLINK_IBD.out
    ).collect()

    GET_SAMPLES_LIST(filters_ch)
    FILTER_SAMPLES(genotype_files, filters_ch)

    PLINK_VARIANTS()

    REPORT()
}

process GET_SAMPLES_LIST {
    publishDir "${outdir_samples}", mode: 'copy'

    input:
    file(filtered_samples_files)

    output:
    file('fail_qc_samples.txt')

    script:
    """
    
    """
}

process FILTER_SAMPLES {
    input:
    tuple val(file_basename), file(genotype_data)
    file(filter_lists)

    output:
    tuple val(file_basename), file("${file_basename}.QCed_samples")

    script:

    """
    plink --bfile $file_basename \
        --remove fail_qc.txt \
        --allow-no-sex --make-bed \
        --out ${file_basename}.QCed_samples
    """
}

process PLINK_IBD {    
    publishDir "${outdir_samples}", mode: 'copy'

    input:
    tuple val(file_basename), file(genotype_data)  

    output:
    file('QC_samples.ibd.genome')

    script:
    """
    #create a subset of frequent SNPs, with low missing rate
    plink --bfile $file_basename  \
        --chr 1-22 --maf 0.25 --geno 0.05 --hwe 1e-6 \
        --allow-no-sex --make-bed --out QC_step2_frequent

    #perform pruning
    plink --bfile QC_step2_frequent \
        --indep-pairwise 50 5 0.2 \
        --allow-no-sex --out QC_step2_prunedSNPs

    #generate IBD report from pruned SNPs
    plink --bfile QC_step2_frequent \
        --extract QC_step2_prunedSNPs.prune.in \
        --genome --allow-no-sex --out QC_samples.ibd
    """
}

process PLINK_SEX_CHECK {
    publishDir "${outdir_samples}", mode: 'copy'
    
    input:
    tuple val(file_basename), file(genotype_data)  

    output:
    file('QC_samples.sex_check.sexcheck')

    script:
    """
    #Rename pseudo autosomale region SNPs
    plink --bfile $file_basename \
        --split-x hg38 no-fail \
        --make-bed --out temp

    #Use plink to perform the sex check
    plink --bfile tem --sex-check 0.5 0.8 --out QC_samples.sex_check
    """
}

process PLINK_SAMPLE_MISS {
    publishDir "${outdir_samples}", mode: 'copy'
    
    input:
    tuple val(file_basename), file(genotype_data)  

    output:
    file('QC_samples.miss.imiss')

    script:
    """
    plink --bfile $file_basename \
        --missing --allow-no-sex \
        --out QC_samples.miss
    """
}

process PLINK_SAMPLE_HET {
    publishDir "${outdir_samples}", mode: 'copy'
    
    input:
    tuple val(file_basename), file(genotype_data)  

    output:
    file('QC_samples.het.het')

    script:
    """
    plink --bfile $file_basename \
        --het --allow-no-sex \
        --out QC_samples.het
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? """\
    
    PIPELINE COMPLETED!    
    ============================================
    QC folder    --> ${params.outdir}
    Alignements  --> ${outdir_align}
    Variants     --> ${outdir_vars}   
    """.stripIndent() : "Oops .. something went wrong" )
}