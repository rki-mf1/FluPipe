rule indexTabix:
    input:
        "{fname}.vcf"
    output:
        temp("{fname}.tabix")
    conda:
        "../envs/tabix.yaml"
    shell:
        r"""
            tabix {input}
        """
