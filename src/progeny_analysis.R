#' This function calculates pathway activity scores using PROGENy for a set of 
#' experiments. Pathway scores are scaled pathway wise
#' 
#' @param df tidy data frame containing gene expression values for
#'   a given set of experiments and their corresponding meta information
#'   (source, organism, accession etc...). 
#' @param M model matrix which is multiplied with the expression values
#' @return tidy data frame of scaled PROGENy results for a given set of 
#' experiments and their corresponding meta information 
run_progeny = function(df, M, ...) {
  meta_df = df %>%
    select(one_of("id", "accession", "pathway", "platform", "info", "treatment",
                  "effect", "source", "sign", "from", "disease",
                  "disease_name", "do_id")) %>%
    distinct()
  
  expr = df %>%
    select(gene, id, expression) %>%
    spread(id, expression, fill=0) %>%
    na.omit() %>%
    data.frame(stringsAsFactors = F, check.names = F)
  rownames(expr) = expr$gene
  expr$gene = NULL
  
  common_genes = intersect(rownames(expr), rownames(M))
  expr_matched = expr[common_genes,,drop=FALSE] %>%
    t()
  M_matched = M[common_genes,, drop=FALSE] %>%
    data.matrix()
  
  stopifnot(names(expr_matched) == rownames(M_matched))
  
  scores = expr_matched %*% M_matched %>%
    data.frame(stringsAsFactors = F, check.names = F) %>%
    mutate(id = rownames(.)) %>%
    gather(key=progeny_pathway, value="score", -id) %>%
    as_tibble() %>%
    inner_join(., meta_df, by="id") %>%
    group_by(progeny_pathway) %>%
    mutate(score = scale(score)) %>%
    ungroup()
}

#' This function ensures that only pathways which are perturbed in mouse and 
#' human data set are considered for the benchmark. If a pathway is only 
#' perturbed in one of the organism the corresponding experiments are discarded
#' 
#' @param df nested data frame containing pathway activity scores in column
#'   "activity"
#' @return same data frame but contains only experiment where the targeted 
#'   pathway was perturbed in mouse and human
filter_common_pws = function(df) {
  df_tmp = df %>%
    select(-activity)
  
  common_pws = df %>%
    unnest(activity) %>%
    distinct(organism, pathway) %>%
    group_by(pathway) %>%
    count() %>%
    filter(n == 2) %>%
    pull(pathway) 
  
  df %>%
    unnest(activity) %>%
    filter(pathway %in% common_pws) %>%
    nest(-c(organism), .key = activity) %>%
    inner_join(df_tmp, by="organism")
}

#' Taking run_progeny() output @seealso [run_progeny()] and prepare input for
#' roc-curve analysis. Column with true response values (0,1) and column with
#' predictor values, which are the PROGENy scores are generated.
#' For pathway inhibiting experiments the PROGENy scores are multiplied with
#' -1. 
#' 
#' @param df run_progeny() output @seealso [run_progeny()]
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for roc curve analysis
prepare_progeny_for_roc = function(df, filter_tn = F) {
  res = df %>%
    mutate(response = case_when(progeny_pathway == pathway ~ 1,
                                progeny_pathway != pathway ~ 0),
           predictor = score * sign)
  if (filter_tn == TRUE) {
    # Only TF which are perturbed and predicted are considered
    z = intersect(res$progeny_pathway, res$pathway)
    res = res %>%
      filter(progeny_pathway %in% z, pathway %in% z)
  }
  res %>%
    select(-c(pathway, score, accession, platform, info, treatment, sign, 
              from)) %>%
    select(pathway = progeny_pathway, everything())
}