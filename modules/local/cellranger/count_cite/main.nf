process CELLRANGER_COUNT_CITE {

    // CITE-seq: gene expression + antibody capture from ONE GEM well, counted
    // together by `cellranger count --libraries --feature-ref`. A separate
    // process rather than a branch inside CELLRANGER_COUNT so that GEX-only
    // samples keep the exact task hash they had before -- a resumed run must
    // not realign every library because this was added.
    //
    // Output layout is identical to CELLRANGER_COUNT (<sample>/outs/*), so QC
    // groups it by the same rule. filtered_feature_bc_matrix.h5 then holds both
    // feature types; the QC notebook splits them (RNA QC on "Gene Expression",
    // "Antibody Capture" carried as an ADT assay).

    tag "Running GEX+ADT on ${sample}"
    label 'process_high'

    input:
        tuple val(sample), path(reads), path(adt_reads, stageAs: 'adt_fastqs/*')
        path(reference)
        path(feature_ref)

    output:
        tuple val(sample), path("${sample}/outs/*"), emit: outs
        path "versions.yml",                         emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
            if [ "${feature_ref.name}" = "NO_FILE_feature_ref" ]; then
                echo "ERROR: sample ${sample} has ADT (CITE-seq) FASTQs but no --feature_reference was given." >&2
                echo "  Pass the antibody feature reference CSV (id,name,read,pattern,sequence,feature_type)." >&2
                exit 1
            fi
            if ! grep -q "Antibody Capture" "${feature_ref}"; then
                echo "ERROR: --feature_reference ${feature_ref.name} has no 'Antibody Capture' row." >&2
                exit 1
            fi

            # Same renaming rule as CELLRANGER_COUNT, applied to each library in
            # its own directory. The ADT library gets its own prefix so the two
            # libraries can never be mistaken for one another by cellranger's
            # <sample>_S<n>_L<nnn>_R<n>_001 matcher.
            cellranger_renaming.py "${sample}" .
            cellranger_renaming.py "${sample}_ADT" adt_fastqs

            mkdir -p gex_fastqs
            for f in ./*.fastq.gz; do ln -s "\$PWD/\${f#./}" gex_fastqs/; done

            cat > libraries.csv <<-END_LIBS
            fastqs,sample,library_type
            \$PWD/gex_fastqs,${sample},Gene Expression
            \$PWD/adt_fastqs,${sample}_ADT,Antibody Capture
            END_LIBS
            sed -i 's/^[[:space:]]*//' libraries.csv
            cat libraries.csv

            cellranger \\
                count \\
                --id="${sample}" \\
                --transcriptome="${reference.name}" \\
                --libraries=libraries.csv \\
                --feature-ref="${feature_ref}" \\
                --localcores=${task.cpus} \\
                --localmem=${task.memory.toGiga()} \\
                --create-bam true ${args}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
                reference: ${reference.name}
                feature_reference: ${feature_ref.name}
            END_VERSIONS
        """

    stub:
        """
            cellranger_renaming.py "${sample}" .
            cellranger_renaming.py "${sample}_ADT" adt_fastqs

            mkdir -p ${sample}/outs/filtered_feature_bc_matrix
            touch ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
            touch ${sample}/outs/filtered_feature_bc_matrix/features.tsv.gz
            touch ${sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
            touch ${sample}/outs/metrics_summary.csv
            touch ${sample}/outs/filtered_feature_bc_matrix.h5
            touch ${sample}/outs/raw_feature_bc_matrix.h5
            touch ${sample}/outs/possorted_genome_bam.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: stub
            END_VERSIONS
        """
}
