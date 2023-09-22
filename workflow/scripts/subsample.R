# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Subsample sequences script
#    (\   / - \     29 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


subsample <- function(ids_file, metadata_file, include_file, exclude_file, 
                      n, method, weights, seed, 
                      output_file) {
  
  set.seed(seed)
  ids <- read_tsv(ids_file)
  df <- read_tsv(metadata_file) %>% 
    filter(sample_id %in% ids$sample_id)

  if (n != -1) {
    
    # Include and exclude specific sequences by id
    if (!is.null(include_file)) {
      include <- read_lines(include_file)
      include_df <- df %>% filter(sample_id %in% include)
    } else include_df <- tibble()
    
    if (!is.null(exclude_file)) {
      exclude <- read_lines(exclude_file)
      df <- df %>%  
        filter(!sample_id %in% exclude)
    }
    
    # Subsample
    if (method == "random") {
      # 1. Random in the full period
      subsample <- df %>%
        sample_n(size = min(n(), n), replace = F)
    
    } else if (method == "cases") {
      # TODO implement proportional subsampling with weights instead of cases explicit
      # 2. Proportional to cases in each week
      cases <- read.csv(cases_file)
      cases_filtered <- cases %>%
        mutate(week = floor_date(ymd(date), "week", week_start = 1)) %>%
        filter(ymd(date) > min(ymd(df$date)), ymd(date) < max(ymd(df$date)), 
               country %in% unique(df$country))
      p_week <- cases_filtered %>%
        group_by(week) %>%
        summarise(cases_week = sum(cases), .groups = "drop") %>%
        mutate(p_cases = cases_week/sum(cases_week)) %>% ungroup
      df_prop <- left_join(df, p_week, by = c("week")) %>%
        replace_na(list(p_cases = 0))
      subsample <- df_prop %>%
        sample_n(size = min(n(), n), replace = F, weight = p_cases)
      
      
    } else if (method == "uniform") {
        
      # 3. Uniform in each week
      df <- df %>% mutate(week = floor_date(date, "week", week_start = 1))
      n_weeks <- nrow(df %>% count(week = floor_date(date, "week", week_start = 1)))
      n_seq_week <- ceiling(n/n_weeks)
      df_seq_week <- df %>% count(week) %>%
        rowwise() %>%
        mutate(n_subsample = min(n, n_seq_week)) %>%
        right_join(df, by = "week")
      subsample <- df_seq_week %>%
        group_by(week) %>%
        sample_n(size = min(n(), n_subsample), replace = F) %>%
        ungroup()
    }
    
    subsample <- bind_rows(subsample, include_df) %>%
      distinct()
  }
  else subsample <- df
  
  # Plot
  # gp <- ggplot(subsample %>% count(date, week)) +
  #   geom_bar(aes(date, n, fill = factor(week)), stat = "identity") +
  #   labs(title = paste("subsampled sequences, method:", method))
  
  print(names(ids))
  print(names(subsample))
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
                              snakemake@params[["n"]],
                              snakemake@params[["method"]],
                              snakemake@params[["weights"]],
                              snakemake@wildcards[["seed"]],
                              snakemake@output[["ids"]]) 

#TODO set color for histogram


