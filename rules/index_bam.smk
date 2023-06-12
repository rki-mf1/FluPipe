rule indexMapping:
    input:
        "{path}/{fname}.sort.bam"
    output:
        temp("{path}/{fname}.sort.bam.bai")
    #log:
    #    os.path.join(DATAFOLDER["logs"], "mapping", "index_{fname}.log")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools index {input} &> {log}
        """
