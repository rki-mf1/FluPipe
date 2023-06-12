rule bgzipVcf:
    input:
        "{fname}.vcf"
    output:
        "{fname}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
            bcftools view -I {input} -O z -o {output}
        """
