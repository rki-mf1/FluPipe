
rule combine_final_refs:
    input:
        segment_refs = expand(os.path.join(DATAFOLDER["mapping"], "{{sample}}","{segment}_final_ref_{{sample}}.fasta"), segment = ALL_SEGMENTS)
    output:
        temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "temp_final_refs_{sample}.fasta"))
    shell:
        """
            cat {input.segment_refs}> {output}
            sed -i 's/kraken:taxid|//' {output}
        """