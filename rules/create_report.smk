rule createReport:
    input:
        trimmed_read_length = os.path.join(DATAFOLDER["trimmed"], "merged_length.tsv"),
        coverage = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.coverage.tsv"), sample = SAMPLES.keys()),
        frag_size  = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.fragsize.tsv"), sample = SAMPLES.keys()),
        mapping_statistics = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.txt"), sample = SAMPLES.keys()),
        mapping_statistics_forR = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt"), sample = SAMPLES.keys()),
        bam_stats_forR = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt"), sample = SAMPLES.keys()),
        version = os.path.join(PROJFOLDER, "pipeline.version"),
        dedupe_metric = expand(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.dedup.txt"), sample = SAMPLES.keys()),
        consensus= expand(os.path.join(IUPAC_CNS_FOLDER, "{sample}.iupac_consensus.fasta"), sample = SAMPLES.keys())
    output:
        report = os.path.join(PROJFOLDER, "qc_report.html"),
        nstats = (os.path.join(DATAFOLDER["reporting"], "N_content_and_Ambigiuos_calls.csv"))
    params:
        p_folder = PROJFOLDER,
        l_folder = DATAFOLDER["reporting"],
        run_id = REPORT_RUNID,
        min_cov = CNS_MIN_COV,
        template = srcdir("../flupipe.Rmd")
    log:
        os.path.join(DATAFOLDER["logs"], "reporting", "reporting.log")
    conda:
        "../envs/r.yaml"
    threads:
        1
    shell:
        # maybe need to replace shell by r call
        r"""
            # create report
            echo "####### compiling report" >> {log}
            VERSION=$(cat {input.version})
            Rscript -e "rmarkdown::render('{params.template}',
                                            params=list(proj_folder='{params.p_folder}', list_folder='{params.l_folder}', run_name='{params.run_id}', min_cov='{params.min_cov}', version='$VERSION'),
                                            output_file=file.path('{output.report}'))" &> {log}
        """
