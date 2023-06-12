rule getBamStatsRefDetec:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "Mapping_2", "{sample}", "{segment}.{sample}.minimap.bam")
    output:
        stats =temp(os.path.join(DATAFOLDER["mapping"], "Mapping_2", "{sample}", "{segment}.{sample}.minimap.bamstats.txt")),
        depth = temp(os.path.join(DATAFOLDER["mapping"], "Mapping_2", "{sample}", "{segment}.{sample}.minimap.bam_depth.txt")),
        sort = temp(os.path.join(DATAFOLDER["mapping"], "Mapping_2","{sample}", "{segment}.{sample}.minimap.sort.bam"))
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools sort {input.bam} > {output.sort}
            samtools coverage {output.sort} > {output.stats}
            samtools depth -a {output.sort} > {output.depth}
        """
