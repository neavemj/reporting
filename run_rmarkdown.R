
library(rmarkdown)
library(optparse)
library(officedown)

# use optparse to grab command line arguments for the report

option_list <- list(
  # required args
  make_option(c("--SAN_table"), type="character", default=NULL,
              help="SAN table from the san_lookup rule", metavar="character"),

  make_option(c("--subtype_table"), type="character", default=NULL,
              help="subtype table READ_COUNTS.txt from irma", metavar="character"),
              
  make_option(c("--CleavageSite"), type="character", default=NULL,
              help="Cleavage Site sequence", metavar="character"),
              
  make_option(c("--tree"), type="character", default=NULL,
              help="Phylogenetic tree in png format from the draw_ggtree rule", metavar="character"),
              
  make_option(c("--HA_gene"), type="character", default=NULL,
              help="HA sequence in fasta format", metavar="character"),
              
  make_option(c("--output"), type="character", default=NULL,
              help="name and file path for word doc", metavar="character"),

  make_option(c("--output_dir"), type="character", default=NULL,
              help="directory to put word doc. Otherwise rmarkdown puts
              the output in the same directory as the markdown document", metavar="character"),
             
  make_option(c("--rmarkdown"), type="character", default=NULL,
              help="location of rmarkdown file", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# need to grab the working project directory here
# otherwise rmarkdown uses the directory where the markdown document
# is located as the working directory. 

working_dir <- getwd()

# now read data for making the report

render(opt$rmarkdown,
    rdocx_document(reference_docx = "../reporting/custom-reference.docx"),
    params = list(
    SAN_table = opt$SAN_table,
    subtype_table = opt$subtype_table,
    Cleavage_site = opt$CleavageSite,
    tree = opt$tree,
    aa_seq_src = opt$HA_gene
    ),
    output_file = opt$output,
    output_dir = opt$output_dir)

