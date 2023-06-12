singularity: "docker://rkibioinf/bwa:0.7.17--a9f152d"

if REFERENCE:
    rule indexBwa:
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE))
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".bwt")
        params:
            os.path.join(DATAFOLDER["mapping"], "{sample}")
        conda:
            "../envs/bwa.yaml"
        shell:
            r"""
                bwa index {input} &> /dev/null
            """
else:
    rule indexBwa:
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta")
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.bwt")
        params:
            os.path.join(DATAFOLDER["mapping"], "{sample}")
        conda:
            "../envs/bwa.yaml"
        shell:
            r"""
                bwa index {input} &> /dev/null
            """