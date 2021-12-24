
library(rmarkdown)
library(optparse)
library(officedown)

# use optparse to grab command line arguments for the report

option_list <- list(
  # required args
  make_option(c("--table1"), type="character", default=NULL,
              help="Table with SAN details and annotation results", metavar="character"),

  make_option(c("--annotation_dir"), type="character", default=NULL,
              help="directory containing the annotation module results", metavar="character"),

  make_option(c("--tree"), type="character", default=NULL,
              help="Phylogenetic tree in png format from the draw_ggtree rule", metavar="character"),
              
  make_option(c("--output"), type="character", default=NULL,
              help="name and file path for word doc", metavar="character"),

  make_option(c("--output_dir"), type="character", default=NULL,
              help="directory to put word doc. Otherwise rmarkdown puts
              the output in the same directory as the markdown document", metavar="character"),
             
  make_option(c("--rmarkdown"), type="character", default=NULL,
              help="location of rmarkdown file", metavar="character"),

  make_option(c("--word_template"), type="character", default=NULL,
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
    rdocx_document(reference_docx = opt$word_template),
    params = list(
    table1 = opt$table1,
	annotation_dir = opt$annotation_dir,
    tree = opt$tree,
    ),
    output_file = opt$output,
    output_dir = opt$output_dir)

