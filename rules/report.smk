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
        san_record = f"san_{list(config['samples'].keys())[0][0:8]}_dump.csv"
    
    params:
        san_endpoint = config["san_query_endpoint"],
        api_key = config["api_key"]
        
    run:
        # Query is performed for a Job
        # All associated sample details are returned
        job_san = list(config['samples'].keys())[0][0:8]
        
        headers = {"Authorization": f"Bearer {params.api_key}"}
        response = requests.get(f"{params.san_endpoint}/{job_san}", headers=headers, verify=False)
        
        # Check that we got a 200 response
        if response.status_code == requests.codes.ok:  # 200
        
            # All samples are returned. Filter for the requested sample
            raw_doc = response.json()
            print(f"Number of samples in job: {len(raw_doc['samples'])}")
            
            #sample = [x for x in raw_doc['samples'] if x['id'] == san_config['san']][0]           
            sample = [x for x in raw_doc['samples'] if x['id'] in config['samples'].keys()]
            
            raw_doc['samples']=[sample]
            
            # Write the result to the expected output file
            with Path(output.san_record).open("w") as fh:
                #json.dump(raw_doc, fh, indent=4)
                # header row
                for s in raw_doc["samples"][0]:
                    s_str = [str(item).strip() for item in s.values()]
                    fh.write(",".join(s_str) + "\n")                
        else:
            print(
                f"Failed to obtain deails for SAN: {config['samples'].keys()}",
                f"({response.status_code}: {response.reason})"
            )
            

rule make_report:
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
       rmarkdown = config["program_dir"] + "reporting/rmarkdown_test.Rmd"
    params:
       subtype_table = expand("02_irma_assembly/{sample}/irma_output/tables/READ_COUNTS.txt", sample=config["samples"])
    output:
       report = "report01.docx"
    shell:
        """
        Rscript {config[program_dir]}/reporting/run_rmarkdown.R \
            --SAN_table {input.SAN_table} \
            --subtype_table {params.subtype_table} \
            --CleavageSite PEKQTR*GLF \
            --tree {input.run_tree} \
            --HA_gene {input.HA_gene} \
            --rmarkdown {input.rmarkdown} \
            --output {output.report}
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
       rmarkdown = config["program_dir"] + "/reporting/rmarkdown_test.Rmd"
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
            --rmarkdown {input.rmarkdown} \
            --output {output.report} \
            --output_dir {config["output_dir"]}
        """