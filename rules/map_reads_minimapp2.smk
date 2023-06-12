
rule minimap2segmentDBs:
    input:
        r1 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
        r2 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz")),
        segment = os.path.join(SEGMENT_DB,"{segment}.all_noIdent_fewAmbig_corLen.fasta")
    output:
        temp(os.path.join(DATAFOLDER["mapping"],"Mapping_1", "{sample}", "{segment}.{sample}.minimap.txt"))
    log:
        temp(os.path.join(DATAFOLDER["mapping"],"Logs", "mapping", "{segment}.{sample}.log"))
    conda:
        "../envs/minimap2.yaml"
    threads:
        50
    shell:
        r"""
            (   time \
                minimap2 -t {threads} -x sr -p 0.99 -B18 -O24,40 -E7,3 -N200 {input.segment} {input.r1} {input.r2} -o {output}
            ) &> {log}
        """