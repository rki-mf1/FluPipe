rule indexPicard:
    input:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "{ref}.fasta")
    output:
        temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{ref}.dict"))
    params:
        os.path.join(DATAFOLDER["mapping"], "{sample}")
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
            gatk CreateSequenceDictionary -R {params}{wildcards.ref}.fasta
        """
