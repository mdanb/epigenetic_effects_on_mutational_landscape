library(optparse)

option_list <- list( 
  make_option("--cancer_type", type="character")
)

args = parse_args(OptionParser(option_list=option_list))

cancer_type = args$cancer_type
