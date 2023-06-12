if REFERENCE:
    rule lofreq_secial_variant_case:
        input:
            iqbam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE)),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".bwt")
        output:
            vcf = (os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.no.filter.vcf.gz")),
            vcf_filter = (os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.no.filter.filtered.vcf.gz"))
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.no.filter.log")
        conda:
            "../envs/lofreq.yaml"
        threads: 10
        shell:
            r"""
                lofreq call-parallel --pp-threads {threads} --no-default-filter -A -B -a 1 -b 1 --verbose --ref {input.ref} --call-indels {input.iqbam} -o {output.vcf} &> {log}
                lofreq filter --snvqual-alpha 0.001 --snvqual-mtc 'fdr' --print-all -i {output.vcf} -o {output.vcf_filter}
            """
else:
    rule lofreq_secial_variant_case:
        input:
            iqbam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta"),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.bwt")
        output:
            vcf = (os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.no.filter.vcf.gz")),
            vcf_filter = (os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.no.filter.filtered.vcf.gz"))
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.no.filter.log")
        conda:
            "../envs/lofreq.yaml"
        threads: 10
        shell:
            r"""
                lofreq call-parallel --pp-threads {threads} --no-default-filter -A -B -a 1 -b 1 --verbose --ref {input.ref} --call-indels {input.iqbam} -o {output.vcf} &> {log}
                lofreq filter --snvqual-alpha 0.001 --snvqual-mtc 'fdr' --print-all -i {output.vcf} -o {output.vcf_filter}
            """