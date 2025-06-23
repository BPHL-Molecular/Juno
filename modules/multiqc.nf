process MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'

    input:
    file summary_report

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data

    script:
    """
    multiqc ${params.output_dir} --interactive
    """
}
