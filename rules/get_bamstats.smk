rule getBamStats:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam")
    output:
        stats = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.txt"),
        covstats = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats_coverage.txt"),
        stats_forR = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping_stats", "{sample}.bamstats.log")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools flagstat {input.bam} 1> {output.stats} 2> {log};
            samtools coverage {input.bam} > {output.covstats};
            cat {output.stats} | sed -e 's/ + /|/' | sed -e 's/ /|/' 1> {output.stats_forR} 2> {log}
        """
