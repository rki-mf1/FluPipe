
rule Rstats_final_refs:
    input:
        stats = expand(os.path.join(DATAFOLDER["mapping"], "Mapping_2", "{{sample}}", "{segment}.{{sample}}.minimap.bamstats.txt"), segment = ALL_SEGMENTS),
        depth = expand(os.path.join(DATAFOLDER["mapping"], "Mapping_2", "{{sample}}", "{segment}.{{sample}}.minimap.bam_depth.txt"), segment = ALL_SEGMENTS)
    output:
        stats_for_top5 = os.path.join(DATAFOLDER["mapping"], "{sample}", "Ranking_top5_Refs_{sample}.txt"),
        name_per_segment = temp(expand(os.path.join(DATAFOLDER["mapping"], "{{sample}}", "{segment}_final_ref_{{sample}}.txt"), segment = ALL_SEGMENTS))
    conda:
        "../envs/R_stats.yaml"
    params:
        script = os.path.join(workflow.basedir, "scripts", "getBestOutOf5.r"),
        path_out = os.path.join(DATAFOLDER["mapping"], "{sample}"),
        path_in = os.path.join(DATAFOLDER["mapping"],"Mapping_2", "{sample}")
    threads:
        50
    shell:
        """
            Rscript {params.script} {params.path_in} {params.path_out}
        """