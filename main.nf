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

process LENGTH_HISTOGRAM {
    tag "$meta.id"
    publishDir "${params.outdir}/qc/histogram", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_length_histogram.tsv"), emit: histogram

    script:
    """
    awk 'NR%4==2 {
        len=length(\$0);
        if (len < 30) short++;
        else if (len < 60) medium++;
        else long++;
    } END {
        printf "short\\t%d\\n", short+0;
        printf "medium\\t%d\\n", medium+0;
        printf "long\\t%d\\n", long+0;
    }' ${fastq} > ${meta.id}_length_histogram.tsv
    """
}

process TRIM_AND_QC {
    tag "$meta.id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq"), emit: trimmed
    tuple val(meta), path("${meta.id}_trimmed_stats.tsv"), emit: stats

    script:
    """
    awk -v minlen=${params.min_length} -v trim3=${params.trim_bases} -v trim5=${params.trim_bases5} '
    BEGIN { RS="@"; ORS="" }
    NR > 1 {
        split(\$0, lines, "\\n")
        header = lines[1]; seq = lines[2]; qual = lines[4]
        len = length(seq)
        keep_len = len - trim3 - trim5
        if (keep_len >= minlen) {
            newseq  = substr(seq, trim5 + 1, keep_len)
            newqual = substr(qual, trim5 + 1, keep_len)
            print "@" header "\\n" newseq "\\n+\\n" newqual "\\n"
        }
    }' ${fastq} > ${meta.id}_trimmed.fastq

    awk 'NR%4==2 {
        n++; len=length(\$0); total+=len;
        gc+=gsub(/[GCgc]/,"",\$0);
    } END {
        printf "num_reads\\t%d\\n", n;
        printf "total_bases\\t%d\\n", total;
        printf "mean_length\\t%.2f\\n", (n>0 ? total/n : 0);
        printf "gc_percent\\t%.2f\\n", (total>0 ? (gc/total)*100 : 0);
    }' ${meta.id}_trimmed.fastq > ${meta.id}_trimmed_stats.tsv
    """
}

process SAMPLE_REPORT {
    tag "$meta.id"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(raw_stats), path(trimmed_stats), path(histogram)

    output:
    path("${meta.id}_report.tsv"), emit: row

    script:
    """
    raw_reads=\$(awk -F'\\t' '\$1=="num_reads"{print \$2}' ${raw_stats})
    trimmed_reads=\$(awk -F'\\t' '\$1=="num_reads"{print \$2}' ${trimmed_stats})
    raw_len=\$(awk -F'\\t' '\$1=="mean_length"{print \$2}' ${raw_stats})
    trimmed_len=\$(awk -F'\\t' '\$1=="mean_length"{print \$2}' ${trimmed_stats})
    raw_gc=\$(awk -F'\\t' '\$1=="gc_percent"{print \$2}' ${raw_stats})
    short=\$(awk -F'\\t' '\$1=="short"{print \$2}' ${histogram})
    medium=\$(awk -F'\\t' '\$1=="medium"{print \$2}' ${histogram})
    long=\$(awk -F'\\t' '\$1=="long"{print \$2}' ${histogram})
    retained=\$(awk -v r="\$raw_reads" -v t="\$trimmed_reads" 'BEGIN { printf "%.2f", (r>0 ? (t/r)*100 : 0) }')

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "\$raw_reads" "\$trimmed_reads" "\$retained" "\$raw_len" "\$trimmed_len" "\$raw_gc" "\$short" "\$medium" "\$long" \\
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
    printf "sample\\traw_reads\\ttrimmed_reads\\tretained_pct\\traw_mean_len\\ttrimmed_mean_len\\traw_gc_pct\\tshort_reads\\tmedium_reads\\tlong_reads\\n" > summary_report.tsv
    sort ${rows} >> summary_report.tsv
    """
}

workflow {
    samples_ch = channel.fromList(
        params.samples.collect { id, path -> [ [id: id], file(path) ] }
    )

    RAW_QC(samples_ch)
    LENGTH_HISTOGRAM(samples_ch)
    TRIM_AND_QC(samples_ch)

    reports_input = RAW_QC.out.stats
        .join(TRIM_AND_QC.out.stats)
        .join(LENGTH_HISTOGRAM.out.histogram)

    SAMPLE_REPORT(reports_input)
    AGGREGATE_REPORT(SAMPLE_REPORT.out.row.collect())
}
