if REFERENCE:
    rule map2reference:
        # add read group tags: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        input:
            PE1 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
            PE2 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz")),
            ref = (os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE))),
            ref_index = (os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".bwt"))
        output:
            temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.bam"))
        params:
            os.path.join(DATAFOLDER["mapping"], "{sample}")
        log:
            os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.log")
        conda:
            "../envs/bwa.yaml"
        threads:
            10
        shell:
            r"""
                (   time \
                    bwa mem -t {threads} -L 10 -Y\
                        -R '@RG\tID:{wildcards.sample}\tPU:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:000' \
                        {input.ref} \
                        {input.PE1} {input.PE2} | \
                        samtools view -Sb -@ {threads} -o {output} \
                ) &> {log}
            """
else:
    rule map2reference:
        # add read group tags: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        input:
            PE1 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
            PE2 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz")),
            ref = (os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta")),
            ref_index = (os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.bwt"))
        output:
            temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.bam"))
        params:
            os.path.join(DATAFOLDER["mapping"], "{sample}")
        log:
            os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.log")
        conda:
            "../envs/bwa.yaml"
        threads:
            10
        shell:
            r"""
                (   time \
                    bwa mem -t {threads} -L 10 -Y\
                        -R '@RG\tID:{wildcards.sample}\tPU:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:000' \
                        {input.ref} \
                        {input.PE1} {input.PE2} | \
                        samtools view -Sb -@ {threads} -o {output} \
                ) &> {log}
        """