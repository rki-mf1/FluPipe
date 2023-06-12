singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"
if REFERENCE:
    rule indexSamtools:
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE))
        output:
            temp(os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".fai"))
        conda:
            "../envs/samtools.yaml"
        shell:
            r"""
                samtools faidx {input}
            """
else:
    rule indexSamtools:
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta")
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.fai")
        conda:
            "../envs/samtools.yaml"
        shell:
            r"""
                samtools faidx {input}
            """
