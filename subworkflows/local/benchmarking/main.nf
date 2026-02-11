/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BENCHMARKING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Compares variant calls against a truth set using hap.py.
    Produces a JSON metrics file consumed by the stimulus optimization framework.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PREPROCESS_QUERY_VCF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_fai), path(fai)

    output:
    tuple val(meta), path("*.norm.vcf.gz"), path("*.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools norm \\
        -m -both \\
        -f ${fasta} \\
        ${vcf} \\
    | bcftools view \\
        -f PASS,. \\
        -e 'ALT="<*>"' \\
        -Oz \\
        -o ${prefix}.norm.vcf.gz

    bcftools index -t ${prefix}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm.vcf.gz
    touch ${prefix}.norm.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process BENCHMARK_HAPPY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hap.py=0.3.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.15--py27h5c5a3ab_0' :
        'community.wave.seqera.io/library/hap.py:0.3.15--76851d3e8624e5b7' }"

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_bed
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_fai), path(fai)

    output:
    tuple val(meta), path("*.summary.csv"), emit: summary_csv
    tuple val(meta), path("*.extended.csv"), optional: true, emit: extended_csv
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_arg = truth_bed ? "-f ${truth_bed}" : ""
    """
    export HGREF=${fasta}

    hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        ${bed_arg} \\
        -r ${fasta} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        --engine vcfeval \\
    || true

    # Generate fallback summary if hap.py did not produce one
    # (e.g. VCF incompatibility).  All metrics are zero so the
    # optimizer learns to avoid this parameter combination.
    if [ ! -f ${prefix}.summary.csv ]; then
        {
            echo "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score"
            echo "INDEL,ALL,0,0,0,0,0,0,0,0,0.000000,0.000000,0.000000,0.000000"
            echo "INDEL,PASS,0,0,0,0,0,0,0,0,0.000000,0.000000,0.000000,0.000000"
            echo "SNP,ALL,0,0,0,0,0,0,0,0,0.000000,0.000000,0.000000,0.000000"
            echo "SNP,PASS,0,0,0,0,0,0,0,0,0.000000,0.000000,0.000000,0.000000"
        } > ${prefix}.summary.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: \$(hap.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo '0.3.15')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_CSV > ${prefix}.summary.csv
    Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio
    INDEL,ALL,100,90,10,95,5,0,0,0,0.900000,0.947368,0.000000,0.923077,,,1.000000,1.000000
    INDEL,PASS,100,90,10,95,5,0,0,0,0.900000,0.947368,0.000000,0.923077,,,1.000000,1.000000
    SNP,ALL,1000,980,20,990,10,0,0,0,0.980000,0.989899,0.000000,0.984925,2.000000,2.000000,1.000000,1.000000
    SNP,PASS,1000,980,20,990,10,0,0,0,0.980000,0.989899,0.000000,0.984925,2.000000,2.000000,1.000000,1.000000
    END_CSV

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: 0.3.15
    END_VERSIONS
    """
}

process EXTRACT_METRICS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(summary_csv)
    val variant_caller
    val aligner

    output:
    tuple val(meta), path("metrics.json"), emit: metrics
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_metrics.py \\
        --summary-csv ${summary_csv} \\
        --sample ${meta.id} \\
        --variant-caller ${variant_caller} \\
        --aligner ${aligner} \\
        --output metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        extract_metrics: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_JSON > metrics.json
    {
        "sample": "${meta.id}",
        "variant_caller": "${variant_caller}",
        "aligner": "${aligner}",
        "snp": {"truth_total": 0, "tp": 0, "fn": 0, "fp": 0, "recall": 0.0, "precision": 0.0, "f1_score": 0.0},
        "indel": {"truth_total": 0, "tp": 0, "fn": 0, "fp": 0, "recall": 0.0, "precision": 0.0, "f1_score": 0.0},
        "summary": {"weighted_f1": 0.0, "total_errors": 0, "snp_f1": 0.0, "indel_f1": 0.0}
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        extract_metrics: 1.0.0
    END_VERSIONS
    """
}

workflow BENCHMARKING {
    take:
    vcf_to_benchmark   // channel: [ val(meta), path(vcf) ]
    truth_vcf           // path: truth VCF file
    truth_bed           // path: truth high-confidence BED regions
    fasta               // channel: [ val(meta), path(fasta) ]
    fasta_fai           // channel: [ val(meta), path(fai) ]
    variant_caller      // val: name of variant caller used
    aligner             // val: name of aligner used

    main:
    ch_versions = Channel.empty()

    // Step 1: Normalize and filter query VCF
    PREPROCESS_QUERY_VCF(
        vcf_to_benchmark,
        fasta,
        fasta_fai,
    )
    ch_versions = ch_versions.mix(PREPROCESS_QUERY_VCF.out.versions.first())

    // Step 2: Run hap.py comparison against truth set
    BENCHMARK_HAPPY(
        PREPROCESS_QUERY_VCF.out.vcf,
        truth_vcf,
        truth_bed,
        fasta,
        fasta_fai,
    )
    ch_versions = ch_versions.mix(BENCHMARK_HAPPY.out.versions.first())

    // Step 3: Extract metrics from hap.py CSV output
    EXTRACT_METRICS(
        BENCHMARK_HAPPY.out.summary_csv,
        variant_caller,
        aligner,
    )
    ch_versions = ch_versions.mix(EXTRACT_METRICS.out.versions.first())

    emit:
    metrics  = EXTRACT_METRICS.out.metrics   // channel: [ val(meta), path(metrics.json) ]
    versions = ch_versions                    // channel: [ path(versions.yml) ]
}
