nextflow.enable.dsl = 2

process SPLIT_WORDS {
    tag "$meta.id"

    input:
    tuple val(meta), val(text)

    output:
    tuple val(meta), path("words.txt"), emit: words

    script:
    """
    echo "${text}" | tr ' ' '\\n' | grep -v '^\$' > words.txt
    """
}

process COUNT_WORDS {
    tag "$meta.id"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(meta), path(words)

    output:
    tuple val(meta), path("${meta.id}_summary.txt"), emit: summary

    script:
    """
    wc -w < ${words} | tr -d ' ' > count.txt
    echo "Sample: ${meta.id}" > ${meta.id}_summary.txt
    echo "Total words: \$(cat count.txt)" >> ${meta.id}_summary.txt
    """
}

process UPPERCASE {
    tag "$meta.id"
    publishDir "${params.outdir}/transformed", mode: 'copy'

    input:
    tuple val(meta), path(words)

    output:
    tuple val(meta), path("${meta.id}_upper.txt"), emit: upper

    script:
    """
    tr '[:lower:]' '[:upper:]' < ${words} > ${meta.id}_upper.txt
    """
}

workflow {
    samples = channel.of(
        [ [id: 'sample1'], params.text1 ],
        [ [id: 'sample2'], params.text2 ]
    )

    SPLIT_WORDS(samples)
    COUNT_WORDS(SPLIT_WORDS.out.words)
    UPPERCASE(SPLIT_WORDS.out.words)
}
