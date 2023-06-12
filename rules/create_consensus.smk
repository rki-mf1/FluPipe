if REFERENCE:
    rule normVcf:
        input:
            vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.gz"),
            tbi = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.tbi"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE)),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE) + ".fai")
        output:
            vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.norm.vcf.gz")),
            index = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.norm.vcf.gz.tbi"))
        conda:
            "../envs/bcftools.yaml"
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.norm.vcf.log")
        threads: 1
        shell:
            r"""
                bcftools norm -O z -c s -f {input.ref} -o {output.vcf} {input.vcf} 2> {log}
                bcftools index -t {output.vcf}
            """
else:
    rule normVcf:
        input:
            vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.gz"),
            tbi = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf.tbi"),
            ref = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta"),
            ref_index = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta.fai")
        output:
            vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.norm.vcf.gz")),
            index = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.norm.vcf.gz.tbi"))
        conda:
            "../envs/bcftools.yaml"
        log:
            os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.norm.vcf.log")
        threads: 1
        shell:
            r"""
                bcftools norm -O z -c s -f {input.ref} -o {output.vcf} {input.vcf} 2> {log}
                bcftools index -t {output.vcf}
            """

#remove indels below freq 90 because it is not possible to write them as ambiguous bases
rule filterVarsConsensus:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.norm.vcf.gz"),
    output:
        vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.vcf.gz"))
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    shell:
        r"""
           bcftools filter -e 'TYPE="indel" && INFO/AF < 0.9' -O z -o {output.vcf} {input.vcf}
        """


rule adjustDeletionConsensus:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.vcf.gz")
    output:
        vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.del.vcf"))
    conda:
        "../envs/pyvcf.yaml"
    script:
        "../scripts/adjust_dels.py"

# low coverage hard filtering
rule createMaskConsensus:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.alnqual.bam"),
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.del.vcf")
    output:
        tmp_bed = temp(os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.lowcov.tmp")),
        final_bed = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.lowcov.bed")
    params:
        cov = CNS_MIN_COV
    conda:
        "../envs/bedtools.yaml"
    singularity:
        "docker://rkibioinf/bedtools:2.29.2--0bfe8ac"
    log:
        os.path.join(DATAFOLDER["logs"], "masking", "{sample}.consensus_mask.log")
    shell:
        r"""
            bedtools genomecov -bga -ibam {input.bam} | awk '$4 < {params.cov}' | bedtools merge > {output.tmp_bed}
            bedtools subtract -a {output.tmp_bed} -b {input.vcf} > {output.final_bed}
        """
#add the special variant case to the low coverage masking bed
rule createMaskConsensus_special_variant_case:
    input:
        bed_lowcov = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.lowcov.bed"),
        bed_special_var = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.special_case_variant_mask.bed")
    output:
        final_bed = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.final.bed")
    shell:
        r"""
            cat {input.bed_lowcov} {input.bed_special_var} > {output}
        """



## genotype adjustment
#explicit call of variants above 90, ambiguous call of SNVs below 90
rule adjustGtConsensus:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.vcf.gz")
    output:
        tmp_vcf = temp(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz.tmp")),
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz")
    conda:
        "../envs/bcftools.yaml"
    params:
        workdir = workflow.basedir
    shell:
        """
            bash {params.workdir}//scripts/add_fake_gt.sh  -i {input.vcf} -g 1 -o {output.tmp_vcf}
            bcftools index {output.tmp_vcf}
            bcftools +setGT {output.tmp_vcf} -- -t q -i 'GT="1" && INFO/AF < 0.9' -n 'c:0/1' | bcftools +setGT -o {output.vcf} -- -t q -i 'GT="1" && INFO/AF >= 0.9' -n 'c:1/1'
            bcftools index {output.vcf}
        """

## create ambig consensus

if REFERENCE:
    rule createAmbiguousConsensus:
        input:
            fasta = os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE)),
            mask = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.final.bed"),
            vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz")
            #vcf_index = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz.tbi")
        output:
            os.path.join(IUPAC_CNS_FOLDER, "{sample}.iupac_consensus.fasta")
        params:
            varlow = VAR_CALL_FRAC, 
            sample = str("{sample}" + "_"),
        log:
            os.path.join(DATAFOLDER["logs"], "consensus", "{sample}.iupac_consensus.tmp.log")
        conda:
            "../envs/bcftools.yaml"
        singularity:
            "docker://rkibioinf/bcftools:1.11--19c96f3"
        shell:
            r"""
                ( bcftools consensus \
                    -H I \
                    -o {output} \
                    -f {input.fasta} \
                    -m {input.mask} \
                    -p {params.sample} \
                    -i 'INFO/AF >= {params.varlow}' \
                    {input.vcf} ) &> {log}
            """
else:
    rule createAmbiguousConsensus:
        input:
            fasta = os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta"),
            mask = os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.final.bed"),
            vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz")
            #vcf_index = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.filtered.gt_adjust.vcf.gz.tbi")
        output:
            os.path.join(IUPAC_CNS_FOLDER, "{sample}.iupac_consensus.fasta")
        params:
            varlow = VAR_CALL_FRAC, 
            sample = str("{sample}" + "_"),
        log:
            os.path.join(DATAFOLDER["logs"], "consensus", "{sample}.iupac_consensus.tmp.log")
        conda:
            "../envs/bcftools.yaml"
        singularity:
            "docker://rkibioinf/bcftools:1.11--19c96f3"
        shell:
            r"""
                ( bcftools consensus \
                    -H I \
                    -o {output} \
                    -f {input.fasta} \
                    -m {input.mask} \
                    -p {params.sample} \
                    -i 'INFO/AF >= {params.varlow}' \
                    {input.vcf} ) &> {log}
            """
