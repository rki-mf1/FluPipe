rule indexTabix:
    input:
        "{fname}.vcf.gz"
    output:
        temp("{fname}.vcf.tbi")
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
            bcftools index -t {input} -o {output}
        """
