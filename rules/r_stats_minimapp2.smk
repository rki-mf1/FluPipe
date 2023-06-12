
rule Rstats_minimap_segment_refs:
    input:
        os.path.join(DATAFOLDER["mapping"],"Mapping_1", "{sample}", "{segment}.{sample}.minimap.txt")
    output:
        temp(os.path.join(DATAFOLDER["mapping"],"5Refs", "{sample}","{segment}.{sample}_best_refs.txt"))
    conda:
        "../envs/R_stats.yaml"
    params:
        script = os.path.join(workflow.basedir, "scripts", "minimap_stats.r"),
        path_out= os.path.join(DATAFOLDER["mapping"],"5Refs","{sample}")
    threads:
        50
    shell:
        """
            Rscript {params.script} {input} {params.path_out}
        """