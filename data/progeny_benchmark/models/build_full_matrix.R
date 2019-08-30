# load original full matrix of PROGENy (11 pathways)
o = get(load("data/progeny_benchmark/models/sub/model_matrix.RData"))$assocs %>%
  as_tibble()

# load extended full matrix with the 3 new pathways and filter for them
e = read_csv("data/progeny_benchmark/models/sub/full_model_extended.csv") %>%
  select(-X1) %>%
  filter(pathway %in% c("Androgen", "Estrogen", "WNT"))

# combine both matrices to a the full matrix containing all 14 pathways
full_matrix = bind_rows(o, e) %>%
  drop_na()

write_csv(full_matrix, 
          "data/progeny_benchmark/models/progeny_matrix_human_full_v1.csv")