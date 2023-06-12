# get sequence depth in single base resolution, an aggregation has to be performed at plotting level
# this data set can also be used to get the amount of bases with >0 sequence depth (genome coverage)

singularity: "docker://rkibioinf/bedtools:2.29.2--0bfe8ac"

rule getCoverage:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        counts = (os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.coverage.tsv"))
    log:
        os.path.join(DATAFOLDER["logs"], "mapping_stats", "{sample}.coverage.log")
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
            (bedtools genomecov -ibam {input.bam} -d 1> {output.counts}) 2> {log}
        """