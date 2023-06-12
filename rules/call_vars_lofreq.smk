if REFERENCE:
    rule callGenomicVariants_lofreq:
        input:
            bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE)),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".bwt")
        output:
            iqbam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam"),
            tbam = temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.viterbi.indel.alnqual.bam")),
            vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.gz"))
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.log")
        params:
            cov = CNS_MIN_COV,
        conda:
            "../envs/lofreq.yaml"
        threads: 10
        shell:
            r"""
            lofreq viterbi --ref {input.ref} {input.bam} | \
                    samtools sort - | \
                    lofreq indelqual --dindel --ref {input.ref} - > {output.tbam}
            lofreq alnqual -b {output.tbam} {input.ref} 1 | samtools sort > {output.iqbam} 2> {log}
            lofreq index {output.iqbam}
            lofreq call-parallel --pp-threads {threads} --ref {input.ref} --call-indels {output.iqbam} -o {output.vcf} &> {log}
            #lofreq call --no-default-filter -A -B -a 1 -b 1  --ref {input.ref} --call-indels {output.iqbam} -o {output.vcf} &> {log}
            """
else:
    rule callGenomicVariants_lofreq:
        input:
            bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta"),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.bwt")
        output:
            iqbam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam"),
            tbam = temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.viterbi.indel.alnqual.bam")),
            vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.gz"))
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.log")
        params:
            cov = CNS_MIN_COV,
        conda:
            "../envs/lofreq.yaml"
        threads: 10
        shell:
            r"""
            lofreq viterbi --ref {input.ref} {input.bam} | \
                    samtools sort - | \
                    lofreq indelqual --dindel --ref {input.ref} - > {output.tbam}
            lofreq alnqual -b {output.tbam} {input.ref} 1 | samtools sort > {output.iqbam} 2> {log}
            lofreq index {output.iqbam}
            lofreq call-parallel --pp-threads {threads} --ref {input.ref} --call-indels {output.iqbam} -o {output.vcf} &> {log}
            #lofreq call --no-default-filter -A -B -a 1 -b 1  --ref {input.ref} --call-indels {output.iqbam} -o {output.vcf} &> {log}
            """
