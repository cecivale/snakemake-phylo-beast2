# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Subsample sequences script
#    (\   / - \     29 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


subsample <- function(ids_file, metadata_file, 
                      groups_file, grouping_var, weights_file,
                      include_file, exclude_file, 
                      n, method, weights, seed, 
                      output_file) {
  
  set.seed(seed)
  ids <- read_tsv(ids_file)
  df <- read_tsv(metadata_file) %>% 
    filter(sample_id %in% ids$sample_id | seq_id %in% ids$seq_id | strain %in% ids$strain)
    
  # Include and exclude specific sequences by id
  if (!is.null(include_file)) {
    include <- read_lines(include_file)
    include_df <- df %>% filter(sample_id %in% include | seq_id %in% include | strain %in% include)
  } else include_df <- tibble()
  
  if (!is.null(exclude_file)) {
    exclude <- read_lines(exclude_file)

    df <- df %>%  
      filter(!sample_id %in% exclude, !seq_id %in% exclude, !strain %in% exclude)
  }
  
  # Subsample
  
  if (method == "random") {
    # 1. Random in the full period
    subsample <- df %>%
      sample_n(size = min(n(), n), replace = F)
  
    }  else if (method == "groups") {
    # 2. Sample a number of samples per each specified group in a file
    groups <- read_tsv(groups_file)
    subsample <- df %>%
      left_join(groups) %>%
      group_by_at(grouping_var) %>%
      sample_n(size = min(n(), count), replace = F)
    
    if (n != -1) {
      subsample <- subsample %>% ungroup %>%
        sample_n(size = min(n(), n), replace = F)
    }
    
  } else if (method == "weights") {
    # 3. Sample according to given weights
    
    if (n >= nrow(df)) subsample <- df
    else {
      weights_specs <- read_tsv(weights_file)
      
      # detect join variable
      join_var <- intersect(names(df), names(weights_specs))
      if (length(join_var) != 1)
        stop("Need exactly one common join column")

      subsample <- df %>%
        left_join(weights_specs, by = join_var) %>%
        replace_na(list(w = 0)) %>%
        group_by_at(join_var) %>%
        mutate(n_var =  round(n * w),
               n_seqs = n(),
               n_take = min(n_var, n_seqs)) %>%
        group_modify(~ .x %>%
                       slice_sample(n = first(.x$n_take))) %>%
        ungroup()
      }
    
  } else if (n == -1) subsample <- df 
    else if (n == 0) subsample <- df %>% slice(0)
    
  subsample <- bind_rows(subsample, include_df) %>%
    ungroup %>%
    distinct()
  
  # Plot
  # gp <- ggplot(subsample %>% count(date, week)) +
  #   geom_bar(aes(date, n, fill = factor(week)), stat = "identity") +
  #   labs(title = paste("subsampled sequences, method:", method))
  
  # print(names(ids))
  # print(names(subsample))
  subsample_to_save <- left_join(subsample, ids) %>%
    select(names(ids))
  
  # Save plot and subsample
  # ggsave(gsub(".tsv", "_plot.pdf", output_file), gp, device = "pdf", width = 5)
  write_tsv(subsample_to_save, file = output_file)
  
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)

# Subsampling ------------------------------------------------------------------
subsample_output <- subsample(ids_file = snakemake@input[["ids"]], 
                              metadata_file = snakemake@input[["metadata"]], 
                              include_file = snakemake@params[["include"]],
                              exclude_file = snakemake@params[["exclude"]],
                              groups_file = snakemake@input[["groups"]],
                              grouping_var =  snakemake@params[["grouping_var"]],
                              weights_file = snakemake@input[["weights"]],
                              n = snakemake@params[["n"]],
                              method = snakemake@params[["method"]],
                              weights = snakemake@params[["weights"]],
                              seed = snakemake@wildcards[["seed"]],
                              output_file = snakemake@output[["ids"]]) 

#TODO set color for histogram


# if (method == "cases") {
#   # TODO implement proportional subsampling with weights instead of cases explicit
#   # 2. Proportional to cases in each week
#   cases <- read.csv(cases_file)
#   cases_filtered <- cases %>%
#     mutate(week = floor_date(ymd(date), "week", week_start = 1)) %>%
#     filter(ymd(date) > min(ymd(df$date)), ymd(date) < max(ymd(df$date)), 
#            country %in% unique(df$country))
#   p_week <- cases_filtered %>%
#     group_by(week) %>%
#     summarise(cases_week = sum(cases), .groups = "drop") %>%
#     mutate(p_cases = cases_week/sum(cases_week)) %>% ungroup
#   df_prop <- left_join(df, p_week, by = c("week")) %>%
#     replace_na(list(p_cases = 0))
#   subsample <- df_prop %>%
#     sample_n(size = min(n(), n), replace = F, weight = p_cases)
#   
# } 

