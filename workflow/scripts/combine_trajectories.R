# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  snakemake-phylo-beast2
#        V\ Y /V    Combine BDMM-Prime trajectory files
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(data.table)
library(optparse)

source("workflow/scripts/utils.R")

combine_trajectories_dt <- function(traj_files, burnin_percentage = 10, subsample_n = NULL) {
  dt_l <- lapply(1:length(traj_files), function(i) {
    traj_file <- traj_files[[i]]
    dt <- fread(
      traj_file,
      fill = T,
      select = c("Sample", "type", "type2", "variable", "value", "age"),
      showProgress = TRUE)
    
    burnin_from <- max(dt$Sample) * (burnin_percentage/100)
    dt <- dt[Sample >= burnin_from]
    # make Sample unique across files
    dt[, Sample := as.character(paste0(i, "_", Sample))]
    
    dt[!is.na(variable)]
  })
  
  dt_combined <- rbindlist(dt_l, use.names = TRUE)
  
  if (!is.null(subsample_n)) {
    available_samples <- unique(dt_combined$Sample)
    n_available <- length(available_samples)
    
    if (subsample_n > n_available) {
      message(
        "Requested subsample_n = ", subsample_n,
        " but only ", n_available,
        " unique trajectories are available after burnin. Keeping all trajectories."
      )
      subsample_n <- n_available
    }
    
    s <- sample(available_samples, subsample_n)
    dt_combined <- dt_combined[Sample %in% s]
  }
  
  return(dt_combined)
}


option_list <- list(
  make_option(
      c("--input"),
      type = "character",
      action = "store",
      help = "Input .traj files as a single space-separated string",
      metavar = "files",
      dest = "input"
  ),
  make_option(
    c("--output"),
    type = "character",
    action = "store",
    help = "Output combined .traj file",
    metavar = "file",
    dest = "output"
  ),
  make_option(
    c("--burnin"),
    type = "double",
    default = 0.1,
    help = "Burn-in fraction [default %default]",
    dest = "burnin"
  ),
  make_option(
    c("--subsample_n"),
    type = "integer",
    default = NULL,
    help = "Optional number of posterior samples to keep",
    dest = "subsample_n"
  )
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input) || is.null(opt$output)) {
  print_help(parser)
  stop("Both --input and --output must be provided.", call. = FALSE)
}

traj_files <- strsplit(opt$input, ",")[[1]]

message("Combining ", length(traj_files), " trajectory files.")
message("Burn-in percentage: ", opt$burnin)

dt_combined <- combine_trajectories_dt(
  traj_files = traj_files,
  burnin_percentage = opt$burnin,
  subsample_n = opt$subsample_n
)

fwrite(dt_combined, file = opt$output, sep = "\t")
message("Wrote combined trajectories to: ", opt$output)

