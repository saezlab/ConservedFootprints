#' Function to read yaml files
#'
#' @param fname string containing file name 
#' @param n 
#'
#' @return parsed yaml files
read_yaml = function(fname, n=Inf) {
  flines = readLines(fname)
  splits = cumsum(flines == "---")
  
  mask = splits <= n
  flines = flines[mask]
  splits = splits[mask]
  
  docs = unname(base::split(flines, splits))
  parsed = lapply(docs, function(d) yaml::yaml.load(paste(d, collapse="\n")))
  return(parsed)
}

#' Function to download experiments and process the raw data and perform 
#' differential expression analysis
#'
#' @param accession string containg accession id of experiment
#' @param data data frame containing information from CREEDs data base
#' @param save_path string containing path where vector of differential 
#'   expressed genes should be saved
#' @return NULL
download_and_process = function(accession, data, save_path, ...) {
  file_name = paste0(accession,".RData")
  path = file.path(save_path, file_name)
  if (file.exists(path)) {
    message(paste(accession, "alreay exists"))
    return(as_tibble(NULL))
  }
  message(accession)
  geo_id = accession
  num = gsub(".*GSE([0-9]*).*", "\\1", geo_id)
  ae_id = paste0("E-GEOD-", num)
  eset_raw = try(ArrayExpress(ae_id))
  
  if (!(is.null(eset_raw) | class(eset_raw) == "try-error")) {
    eset = try(eset_raw %>%
                 ma$qc() %>%
                 ma$normalize() %>%
                 ma$annotate(summarize="hgnc_symbol"))
    
    if (!(is.null(eset) | class(eset) == "try-error")) {
      x = data %>% 
        mutate(accession = accession) %>%
        nest(-id) %>%
        mutate(eset = list(eset),
               load_data = F) %>%
        mutate(expr = pmap(., .f=match_samples)) %>%
        dplyr::select(id, expr) %>%
        mutate(limma = expr %>% map(calc_contrast)) %>%
        unnest(limma)
      
      if (nrow(x) != 0) {
        file_name = paste0(geo_id,".RData")
        save(x, file = file.path(save_path,file_name))
      }
    }
  }
}

#' This function matches samples from expression set to samples from CREEDs
#'
#' @param data data frame containing information from CREEDs data base
#' @param eset expression set with normalized expression values
#' @param load_data logical flag indicating if expression set is already 
#'   downloaded so it can be loaded from local disc
#' @return expression set with samples available in CREEDs data base
match_samples = function(data, eset=NULL, load_data = T, ...) {
  if (load_data == T) {
    print("load data")
    path = data %>% distinct(expr_path) %>% pull()
    platform = data %>% distinct(platform) %>% pull()
    eset = get(load(path))
    if (class(eset) == "list") {
      eset = eset[[platform]]
    }
  }
  print("hi")
  colnames(eset) = sub("(GSM[[:digit:]]*).*","\\1" ,colnames(eset))
  
  expr = eset %>%
    exprs() %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>%
    gather(sample, expression, -gene)
  
  merged_data = inner_join(data, expr, by="sample") %>%
    dplyr::select(-one_of("expr_path"))
  
  if (nrow(merged_data) == 0 | nrow(distinct(merged_data, group)) == 1) {
    return(as_tibble(NULL))
  }
  
  tmp = merged_data %>%
    distinct(group, sample) %>%
    group_by(group) %>%
    count() %>%
    spread(group, n) %>%
    mutate(keep = case_when(control >= 2 & perturbed >= 1 ~ "yes",
                            TRUE ~ "no"))
  
  if (tmp$keep == "no") {
    return(as_tibble(NULL))
  } else {
    return(merged_data)
  }
}

#' This function performs differential expression analysis
#'
#' @param expression set
#' @return vector with differential expressed genes
calc_contrast = function(expr) {
  if (nrow(expr) == 0) {
    return(as_tibble(NULL))
  }
  message(unique(expr$accession))
  
  meta = expr %>%
    dplyr::select(-c("gene", "expression", "sample", "group")) %>%
    distinct() %>%
    mutate(join_var = 1)
  
  
  # create target
  targets = expr %>%
    distinct(sample, group) 
  
  # create design matrix
  groups = targets %>%
    pull(group) %>%
    as_factor()
  
  emat = expr %>%
    dplyr::select(sample, gene, expression) %>%
    spread(sample, expression) %>%
    data.frame(row.names=1, check.names = F, stringsAsFactors = F)
  emat = emat[, targets$sample]
  
  design = model.matrix(~0+groups)
  colnames(design) = levels(groups)
  
  fit = lmFit(emat, design)
  
  # define contrasts
  contrasts = makeContrasts(perturbed-control, levels=design)
  
  fit2 = contrasts.fit(fit, contrasts)
  limma_result = eBayes(fit2) %>% 
    topTable(n=Inf) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    dplyr::select(p.val = P.Value, adj.p.val = adj.P.Val, everything()) %>%
    arrange(gene)
  
  
  mylogFC = expr %>%
    group_by(gene, group) %>%
    summarise(m = mean(expression)) %>%
    spread(group, m) %>%
    mutate(logFC = perturbed-control) %>%
    arrange(gene) %>%
    dplyr::select(gene, logFC)
  
  # check if the correct contrast is computed
  stopifnot(round(mylogFC$logFC,6) == round(limma_result$logFC, 6))
  
  # z-scores
  control_samples = targets %>%
    filter(group == "control") %>%
    pull(sample)
  
  perturbed_samples = targets %>%
    filter(group == "perturbed") %>%
    pull(sample)
  
  control = emat[, control_samples, drop=F]
  perturbed = emat[, perturbed_samples, drop=F]
  mean_control= apply(control, 1, mean)
  sd_control = apply(control, 1, sd)
  mean_perturbed = apply(perturbed, 1, mean)
  logFC = mean_perturbed - mean_control
  model = loess(sd_control ~ mean_control)
  
  z = logFC / predict(model, mean_perturbed)
  
  z_df = tibble(gene = names(z), z = unname(z)) 
  tmp = full_join(limma_result, z_df, by="gene") %>%
    mutate(join_var = 1) %>%
    inner_join(meta, by="join_var")%>%
    dplyr::select(-c(join_var, AveExpr, B))
}
