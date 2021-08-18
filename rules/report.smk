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

rule make_report:
    message:
        """
        ** report **
        Making the report in R
        """
    input:
        SAN = "21-00000",
        SampleID = "21-000-ID",
        Host = "wild bird",
        H_type = "H7",
        CleavageSite = "PEIPGKR*GLF",
        tree = "/RFolder/Resources/H7N7_20-02853_tree.png",
        HA_gene = "/RFolder/Resources/HA_alignment.align.fasta",
        datasets = "NAIWB sequence datasets and recent Victorian H7 (chicken and emu) outbreaks"
    output:
        report = "../reporting/reportTEST_01.docx"
    shell:
        """
        Rscript {config[program_dir]}/reporting/ReportTest_RMarkdown.Rmd \
            --SAN {input.SAN} \
            --SampleID {input.SampleID} \
            --Host {input.Host} \
            --H_type {input.H_type) \
            --CleavageSite {input.CleavageSite} \
            --tree {input.tree} \
            --HA_gene {input.HA_gene} \
            --datasets {input.datasets} \
            --output {output.report}
        """