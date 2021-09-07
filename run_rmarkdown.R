
library(rmarkdown)
library(optparse)
library(officedown)

# use optparse to grab command line arguments for the report

option_list <- list(
  # required args
  make_option(c("--SAN"), type="character", default=NULL,
              help="ACDP SAN number of job", metavar="character"),
  make_option(c("--SampleID"), type="character", default=NULL,
              help="ACDP Sample IDs analysed", metavar="character"),
  make_option(c("--Host"), type="character", default=NULL,
              help="Type of animal sample derived from", metavar="character"),
  make_option(c("--H_type"), type="character", default=NULL,
              help="H type identified", metavar="character"),
  make_option(c("--CleavageSite"), type="character", default=NULL,
              help="Cleavage Site sequence", metavar="character"),
  make_option(c("--tree"), type="character", default=NULL,
              help="Phylogenetic tree in image format", metavar="character"),
  make_option(c("--HA_gene"), type="character", default=NULL,
              help="HA sequence in fasta format", metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="name and file path for word doc", metavar="character"),
  make_option(c("--rmarkdown"), type="character", default=NULL,
              help="location of rmarkdown file", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# now read data for making the report

render(opt$rmarkdown,
    rdocx_document(reference_docx = "/flush5/sco308/aiv_pipeline/reporting/custom-reference.docx"),
    params = list(
    SAN = opt$SAN,
    Sample_ID = opt$SampleID,
    Host = opt$Host,
    H_type = opt$H_type,
    Cleavage_site = opt$CleavageSite,
    phylo_tree_Src = opt$tree,
    aa_seq_src = opt$HA_gene
    ),
    output_file = opt$output)

