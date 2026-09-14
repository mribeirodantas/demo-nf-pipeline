nextflow.enable.dsl = 2

process RAW_QC {
    tag "$meta.id"
    publishDir "${params.outdir}/qc/raw", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_raw_stats.tsv"), emit: stats

    script:
    """
    awk 'NR%4==2 {
        n++; len=length(\$0); total+=len;
        gc+=gsub(/[GCgc]/,"",\$0);
    } END {
        printf "num_reads\\t%d\\n", n;
        printf "total_bases\\t%d\\n", total;
        printf "mean_length\\t%.2f\\n", (n>0 ? total/n : 0);
        printf "gc_percent\\t%.2f\\n", (total>0 ? (gc/total)*100 : 0);
    }' ${fastq} > ${meta.id}_raw_stats.tsv
    """
}

process TRIM_FILTER {
    tag "$meta.id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq"), emit: trimmed

    script:
    """
    awk -v minlen=${params.min_length} -v trim=${params.trim_bases} '
    BEGIN { RS="@"; ORS="" }
    NR > 1 {
        split(\$0, lines, "\\n")
        header = lines[1]; seq = lines[2]; qual = lines[4]
        len = length(seq)
        if (len - trim >= minlen) {
            newseq  = substr(seq, 1, len - trim)
            newqual = substr(qual, 1, len - trim)
            print "@" header "\\n" newseq "\\n+\\n" newqual "\\n"
        }
    }' ${fastq} > ${meta.id}_trimmed.fastq
    """
}

process TRIMMED_QC {
    tag "$meta.id"
    publishDir "${params.outdir}/qc/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_trimmed_stats.tsv"), emit: stats

    script:
    """
    awk 'NR%4==2 {
        n++; len=length(\$0); total+=len;
        gc+=gsub(/[GCgc]/,"",\$0);
    } END {
        printf "num_reads\\t%d\\n", n;
        printf "total_bases\\t%d\\n", total;
        printf "mean_length\\t%.2f\\n", (n>0 ? total/n : 0);
        printf "gc_percent\\t%.2f\\n", (total>0 ? (gc/total)*100 : 0);
    }' ${fastq} > ${meta.id}_trimmed_stats.tsv
    """
}

process SAMPLE_REPORT {
    tag "$meta.id"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(raw_stats), path(trimmed_stats)

    output:
    path("${meta.id}_report.tsv"), emit: row

    script:
    """
    raw_reads=\$(awk -F'\\t' '\$1=="num_reads"{print \$2}' ${raw_stats})
    trimmed_reads=\$(awk -F'\\t' '\$1=="num_reads"{print \$2}' ${trimmed_stats})
    raw_len=\$(awk -F'\\t' '\$1=="mean_length"{print \$2}' ${raw_stats})
    trimmed_len=\$(awk -F'\\t' '\$1=="mean_length"{print \$2}' ${trimmed_stats})
    raw_gc=\$(awk -F'\\t' '\$1=="gc_percent"{print \$2}' ${raw_stats})
    retained=\$(awk -v r="\$raw_reads" -v t="\$trimmed_reads" 'BEGIN { printf "%.2f", (r>0 ? (t/r)*100 : 0) }')

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "\$raw_reads" "\$trimmed_reads" "\$retained" "\$raw_len" "\$trimmed_len" "\$raw_gc" \\
        > ${meta.id}_report.tsv
    """
}

process AGGREGATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(rows)

    output:
    path("summary_report.tsv")

    script:
    """
    printf "sample\\traw_reads\\ttrimmed_reads\\tretained_pct\\traw_mean_len\\ttrimmed_mean_len\\traw_gc_pct\\n" > summary_report.tsv
    sort ${rows} >> summary_report.tsv
    """
}

workflow {
    samples_ch = channel.fromList(
        params.samples.collect { id, path -> [ [id: id], file(path) ] }
    )

    RAW_QC(samples_ch)
    TRIM_FILTER(samples_ch)
    TRIMMED_QC(TRIM_FILTER.out.trimmed)

    reports_input = RAW_QC.out.stats.join(TRIMMED_QC.out.stats)

    SAMPLE_REPORT(reports_input)
    AGGREGATE_REPORT(SAMPLE_REPORT.out.row.collect())
}
