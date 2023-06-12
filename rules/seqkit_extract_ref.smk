
rule seqkit_extract_ref:
    input:
        reftxt = os.path.join(DATAFOLDER["mapping"],"5Refs", "{sample}", "{segment}.{sample}_best_refs.txt"),
        segmentDB= os.path.join(SEGMENT_DB,"{segment}.all_noIdent_fewAmbig_corLen.fasta")
    output:
        temp(os.path.join(DATAFOLDER["mapping"],"5Refs", "{sample}","{segment}.{sample}_best_refs.fasta"))
    conda:
        "../envs/seqkit.yaml"
    threads:
        50
    shell:
        """
            seqkit grep -f {input.reftxt} {input.segmentDB} --quiet > {output}
        """
