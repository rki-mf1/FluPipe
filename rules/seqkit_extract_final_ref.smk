rule seqkit_extract_final_ref:
    input:
        reftxt = os.path.join(DATAFOLDER["mapping"], "{sample}", "{segment}_final_ref_{sample}.txt"),
        segmentDB= os.path.join(SEGMENT_DB,"{segment}.all_noIdent_fewAmbig_corLen.fasta")
    output:
        temp(os.path.join(DATAFOLDER["mapping"], "{sample}","{segment}_final_ref_{sample}.fasta"))
    conda:
        "../envs/seqkit.yaml"
    threads:
        50
    shell:
        """
            seqkit grep -f {input.reftxt} {input.segmentDB} --quiet > {output}
        """
