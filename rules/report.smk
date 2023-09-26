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


import textwrap
import requests
import json
from pathlib import Path
import glob
import csv


rule all_report:
    input:
        expand("05_report/{sample}_blast.csv", sample=config["samples"])


rule draw_dag:
    input:
        "AIV_pipe.snakefile"
    output:
        "benchmarks/dag.png"
    shell:
        "snakemake -s {input} --rulegraph 2> /dev/null | dot -T png > {output}"


rule get_package_versions:
    input:
        config["program_dir"] + "qc/config/software_list.txt"
    output:
        "logs/software_versions.txt"
    params:
        # for some reason, I have to define the sed pattern here, then pass it in
        # otherwise it is weirdly expanded in the actual shell command
        sed_pat1 = r"s/\(.*\)/^\1\t/g",
        sed_pat2 = r"s/ \+/\t/g"
    shell:
        # lists all packages in the conda environment
        # replaces large whitespace with tabs
        # greps for specific packages with ^ and \t to ensure complete match
        # cuts just the first 2 columns of interest
        """
        conda list | \
            sed "{params.sed_pat2}" | \
            grep -f <(cat {input} | sed "{params.sed_pat1}") | \
            cut -f 1,2 > {output}
        """


rule san_lookup:
    message: textwrap.dedent(f"""
        Query Sample Manager for SAN job/sample details
        via the SanQuery service by tunnelling through MySequences.
        This convoluted approach maintains security when accessing 
        Sample Manager.
        
        The resulting JSON file will contain the Job details along 
        with a list of all samples for that job.
        
        The list can be filtered for the desired sample
    """).rstrip()


    output: 
        #san_record = f"san_{list(config['samples'].keys())[0][0:8]}_dump.csv"
        san_record = "03_annotation/san_{sample}_dump.csv"

    params:
        san_endpoint = config["san_query_endpoint"],
        api_key = config["api_key"]
        
    run:
        # Query is performed for a Job
        # All associated sample details are returned
        #job_san = list(config['samples'].keys())[0][0:8]
        #job_san = {wildcards.sample}[0:8]
        job_san = wildcards.sample[0:8]
        
        headers = {"Authorization": f"Bearer {params.api_key}"}
        response = requests.get(f"{params.san_endpoint}/{job_san}", headers=headers, verify=False)
        
        # Check that we got a 200 response
        if response.status_code == requests.codes.ok:  # 200
        
            # All samples are returned. Filter for the requested sample
            raw_doc = response.json()
            print(f"Number of samples in job: {len(raw_doc['samples'])}")
          
            sample = [x for x in raw_doc['samples'] if x['id'] in wildcards.sample]
            
            raw_doc['samples']=[sample]
            
            # Write the result to the expected output file
            with open(output.san_record, "w") as fh:
                # header row
                for s in raw_doc["samples"][0]:
                    s_str = [str(item).strip() for item in s.values()]
                    fh.write(",".join(s_str) + "\n")                
        else:
            print(
                f"Failed to obtain details for SAN: {wildcards.sample}",
                f"({response.status_code}: {response.reason})"
            )
            

rule combine_annotations:
    # this rule combines annotations (e.g. SAN info, cleavage site, etc)
    # from a single sample
    input:
        run_annotation = "03_annotation/{sample}/ANNOTATION_COMPLETE",
        san_record = "03_annotation/san_{sample}_dump.csv",
    params:
        annotation_dir = "03_annotation/{sample}/",
    output:
        table1 = "05_report/table1_{sample}.csv"
    run:
        # combine SAN and annotation info into a table for the report
        # get SAN table info 
        with open(input.san_record) as infile:
            # SAN file looks like this
            # ['21-02023-0001', '21-02391-055', 'WB', 'SWAB', '780676', '0\n']
            # grab SAN, submitter ID, host and sample type
            SAN_info = infile.read().split(",")[0:4]
        # get cleavage site
        cleavage_table = glob.glob(params.annotation_dir + "*cs.txt")
        # assuming that there is only one file detected and opening
        # the first item in the list
        with open(cleavage_table[0]) as infile:
            # cleavage site is first item in the table
            SAN_info.append(infile.read().split(" ")[0])
        # get subtype from the annotation dir
        # this will throw an error if no HA file is produced
        H_type = glob.glob(params.annotation_dir + "A_HA_*")[0].split(".")[0].split("_")[3]
        N_type = glob.glob(params.annotation_dir + "A_NA_*")[0].split(".")[0].split("_")[3]
        SAN_info.append(H_type + N_type)
        
        # write out info into a table for the report
        with open(output.table1, "w") as f:
            f.write(",".join(SAN_info) + "\n")
       
       
rule combine_samples:
    # this rule takes individual sample annotations and cat's them 
    # together for the report
    input:
        expand("05_report/table1_{sample}.csv", sample=config["samples"])
    output:
        "05_report/table1.csv"
    shell:
        """
        cat {input} > {output}
        """


rule combine_blast_tables:
    # this rule will combine blast results for each segment for a single sample
    input:
        run_annotation = "03_annotation/{sample}/ANNOTATION_COMPLETE",
    params:
        annotation_dir = "03_annotation/{sample}/",
    output:
        blast_table = "05_report/{sample}_blast.csv"
    run:        
        # blast tables in the annotation directory
        # sorting the results to get some sort of order into the table
        blast_tables = sorted(glob.glob(params.annotation_dir + "*blast.csv"))
        # should return 8 tables if each segment is in there
        table_rows = []
        for tb in blast_tables:
            with open(tb) as infile:
                # blast tables look like this:
                # ,qseqid,sseqid,qlen,slen,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
                #0,A_PB1,gb|CY015025.1|,2341,2341,99.487,780,4,0,2,2341,1,2340,0.0,1854
                header = next(infile)
                for line in infile:
                    # just grab the query id (col2), subject id (col3), pident (col6) and length (col7)
                    cols = [line.split(",")[i] for i in [1,2,5,6]]
                    table_rows.append(cols)          
        # write out info into a table for the report
        with open(output.blast_table, "w") as f:
            writer = csv.writer(f)
            writer.writerows(table_rows)


rule make_report:
    message:
        """
        ** report **
        Making the report in R
        """
    input:
        #SAN_table = f"san_{list(config['samples'].keys())[0][0:8]}_dump.csv",
        table1 = "05_report/table1.csv",
        run_annotation = expand("03_annotation/{sample}/ANNOTATION_COMPLETE", sample=config["samples"]),
        run_tree = expand("04_phylogenetics/{sample}_tree_finished.txt", sample=config["samples"]),       
    params:
        #subtype_table = expand("02_irma_assembly/{sample}/irma_output/tables/READ_COUNTS.txt", sample=config["samples"]),
        annotation_dir = expand("03_annotation/{sample}/", sample=config["samples"]),
        rmarkdown = config["program_dir"] + "reporting/scripts/build_report.Rmd",
        word_template = config["program_dir"] + "reporting/scripts/custom-reference.docx",

    output:
        report = "report.docx"
    shell:
        """
        Rscript {config[program_dir]}/reporting/scripts/run_rmarkdown.R \
            --table1 {input.table1} \
            --annotation_dir {params.annotation_dir} \
            --rmarkdown {params.rmarkdown} \
            --word_template {params.word_template} \
            --output {output.report} \
            --output_dir .
        """
            

rule make_report_testing:
    message:
        """
        ** report **
        Making the report in R
        """
    input:
       SAN_table = f"san_{list(config['samples'].keys())[0][0:8]}_dump.csv",
       run_IRMA = expand("02_irma_assembly/{sample}/IRMA_COMPLETE", sample=config["samples"]),
       run_tree = expand("04_phylogenetics/{sample}_tree_finished.txt", sample=config["samples"]),
       HA_gene = "../introduction/raw_data/HA_alignment.align.fasta",
      # SampleID = "21-000-ID",
      # Host = "wild bird",
      # H_type = "H7",
      # CleavageSite = "PEIPGKR*GLF",
      # run_tree = "/flush5/sco308/aiv_pipeline/introduction/raw_data/H7N7_20-02853_tree.png"
       rmarkdown = config["program_dir"] + "/reporting/rmarkdown_test.Rmd",
       BLAST_table = "../introduction/raw_data/BLAST examples/blastresults.HA.csv"
    params:
       subtype_table = expand("02_irma_assembly/{sample}/irma_output/tables/READ_COUNTS.txt", sample=config["samples"])
    output:
        report = "reportTEST_01.docx"
    shell:
        """
        Rscript {config[program_dir]}/reporting/run_rmarkdown.R \
            --SAN_table 21-0000 \
            --subtype_table {params.subtype_table} \
            --CleavageSite PEKQTR*GLF \
            --tree {input.run_tree} \
            --HA_gene {input.HA_gene} \
            --BLAST_table {input.BLAST_table}.\
            --rmarkdown {input.rmarkdown} \
            --output {output.report} \
            --output_dir .
        """