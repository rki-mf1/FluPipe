singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"

rule getInsertSize:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        sizes = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.fragsize.tsv")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping_stats", "{sample}.fragsize.log")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            echo 0 1> {output.sizes}
            samtools view -F 4 {input.bam} | cut -f9 1>> {output.sizes} 2> {log}
        """