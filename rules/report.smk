"""
This will build the report

The inputs are:
    - SAN number
    - Sample ID
    - Host
    - H type
    - Cleavage site sequence
    - HA gene aa sequence in fasta format
    - Phylogenetic tree in image format
    - Location of the cleavage site
    - Which datasets were used for comparison purposes (e.g. NAIWB, Genback)
    - 
The outputs are:
    - a report in Microsoft word format
"""

configfile: "config.yaml"

rule make_report:
    message:
        """
        ** report **
        Making the report in R
        """
    input:
      # SAN = "21-00000",
      # SampleID = "21-000-ID",
      # Host = "wild bird",
      # H_type = "H7",
      # CleavageSite = "PEIPGKR*GLF",
       tree = "/flush5/sco308/aiv_pipeline/introduction/raw_data/H7N7_20-02853_tree.png",
       HA_gene = "/flush5/sco308/aiv_pipeline/introduction/raw_data/HA_alignment.align.fasta",
       rmarkdown = config["program_dir"] + "/reporting/rmarkdown_test.Rmd"
    output:
        report = "/flush5/sco308/aiv_pipeline/testing/reportTEST_01.docx"
    shell:
        """
        Rscript {config[program_dir]}/reporting/run_rmarkdown.R \
            --SAN 21-00344 \
            --SampleID 34 \
            --Host Duck \
            --H_type H7 \
            --CleavageSite PEKQTR*GLF \
            --tree {input.tree} \
            --HA_gene {input.HA_gene} \
            --rmarkdown {input.rmarkdown} \
            --output {output.report}
        """