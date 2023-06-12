
rule R_filter_variants_special_variant_case:
    input:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.no.filter.filtered.vcf.gz")
    output:
        os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.special_case_variant_mask.bed")
    conda:
        "../envs/R_stats.yaml"
    params:
        script = os.path.join(workflow.basedir, "scripts", "filter_good_quality_failed_variants_special_case.r"),
        path_out= os.path.join(DATAFOLDER["masking"], "{sample}")
    threads:
        10
    shell:
        """
            Rscript {params.script} {input} {params.path_out}
        """