# This script intents to translate dorothea's data table to the format 
# required by viper format
library(tidyverse)
source("src/dorothea_analysis.R")

dest_path = "data/dorothea_benchmark/regulons/regulons_in_viper_format"
confidence_level_combi = c("A", "B", "C", "D", "E", 
                           "AB", "ABC", "ABCD", "ABCDE")
# human
h = read_csv("data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")

map(confidence_level_combi, function(c) {
  # dissect confidence level combinations
  splitted_c = c %>% str_split(pattern = "") %>% pluck(1)
  
  # transform data table to viper format
  viper_format = h %>%
    filter(confidence %in% splitted_c) %>%
    df2regulon()
  
  # construct save name path
  save_name = str_c("dorothea_human_", c, "_viper_format.rds")
  save_path = file.path(dest_path, save_name)
  saveRDS(viper_format, save_path)
})


# mouse
m = read_csv("data/dorothea_benchmark/regulons/dorothea_regulon_mouse_v1.csv")

map(confidence_level_combi, function(c) {
  # dissect confidence level combinations
  splitted_c = c %>% str_split(pattern = "") %>% pluck(1)
  
  # transform data table to viper format
  viper_format = m %>%
    filter(confidence %in% splitted_c) %>%
    df2regulon()
  
  # construct save name path
  save_name = str_c("dorothea_mouse_", c, "_viper_format.rds")
  save_path = file.path(dest_path, save_name)
  saveRDS(viper_format, save_path)
})
