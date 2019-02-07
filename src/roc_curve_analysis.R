#' This function calculates precision recall curves
#' 
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column response
#' @return tidy data frame containing recall, precision, auc, tp, tn and 
#'   coverage
calc_pr_curve = function(df) {
  feature_coverage = df %>% 
    select(one_of("pathway", "tf")) %>% 
    distinct() %>% 
    nrow()
  
  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)
  
  r = pr.curve(scores.class0 = df$predictor,
               weights.class0 = df$response,
               curve=T)
  res = r$curve %>%
    as_tibble() %>%
    setNames(., c("recall", "precision", "th")) %>%
    mutate(auc = r$auc.davis.goadrich,
           type = r$type,
           n = sum(df$response),
           tp = nrow(tp),
           tn = nrow(tn),
           coverage = feature_coverage) %>%
    arrange(recall, desc(precision))
}

#' This function calculates receiver operating characteristic
#'
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column respons 
#' @param downsampling logical flag indicating if the number of TN should be 
#'   downsampled to the number of TP
#' @param times integer showing the number of downsampling
#' @return tidy data frame with tpr, fpr, auc, n, tp, tn and coverage
calc_roc_curve = function(df, downsampling=F, times = 1000) {
  if (sum(df$response) == 0) {
    return(as_tibble(NULL))
  } 
  
  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)
  
  feature_coverage = df %>% 
    select(one_of("pathway", "tf")) %>% 
    distinct() %>% 
    nrow()
  
  if (downsampling == T) {
    # number of true positives
    num_tp = nrow(tp)
    
    res = map_df(seq(from=1, to=times, by=1), .f= function(i) {
      df_sub = sample_n(tn, num_tp, replace=TRUE) %>% 
        bind_rows(tp)
      
      r_sub = roc(response = df_sub$response, predictor = df_sub$predictor, direction = "<")
      
      res_sub = tibble(tpr = r_sub$sensitivities,
                   fpr = 1-r_sub$specificities,
                   th = r_sub$thresholds,
                   auc = r_sub$auc,
                   n = sum(df$response),
                   tp = nrow(tp),
                   tn = nrow(tn),
                   coverage = feature_coverage) %>%
        mutate_("run" = i)

    })
  } else {
    r = roc(response = df$response, predictor = df$predictor, direction = "<")
    
    res = tibble(tpr = r$sensitivities,
                 fpr = 1-r$specificities,
                 th = r$thresholds,
                 auc = r$auc,
                 n = sum(df$response),
                 tp = nrow(tp),
                 tn = nrow(tn),
                 coverage = feature_coverage) %>%
      arrange(fpr, tpr)
    
    ci = calc_ci_auc(r)
    res = res %>%
      mutate(sig_level = ci$sig_level,
             ci = ci$ci)
  }
  return(res)
}

#' This function calculates confidence intervalls for area under the curve (auc)
#' values for different intervalls (90, 95, 99, 99.9 and 99.99 %)
#' 
#' @param r output (object) from pROC::roc() function
#' @return tibble with confidence intervall, upper bound, lower bound  and
#'   significance level (symbolic)
calc_ci_auc = function(r) {
  sig_levels = tibble(pval = c(0.1, 0.05, 0.01, 0.001),
                      symbol = c("+", "*", "**", "***"))
  
  sig_test = map_df(c(0.9, 0.95, 0.99, 0.999), function(conf_level) {
    x = ci(r, of="auc", conf.level=conf_level)
    list(ci = conf_level, lb = x[1], ub = x[3])
  })
  
  if (r$auc >= 0.5) {
    sig_symbol = sig_test %>%
      filter(lb >= 0.5) %>%
      slice(which.min(lb)) %>%
      mutate(
        sig_level = with(sig_levels,
                         as.character(symbol)[match(round(1-ci,8), pval)])
      ) %>%
      select(sig_level) %>%
      mutate(ci = list(sig_test))
  } else if (r$auc < 0.5) {
    sig_symbol = sig_test %>%
      filter(ub < 0.5) %>%
      slice(which.max(ub)) %>%
      mutate(
        sig_level = with(sig_levels, 
                         as.character(symbol)[match(round(1-ci,8), pval)])
      ) %>%
      select(sig_level) %>%
      mutate(ci = list(sig_test))
  }
  if (nrow(sig_symbol) == 0) {
    sig_symbol = tibble(sig_level = "", ci = list(sig_test))
  } else {
    return(sig_symbol)
  }
}

#' Returns ROC object for further analysis
#'
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column respons 
#' @return ROC object
get_roc_object = function(df) {
  if (sum(df$response) == 0) {
    return(as_tibble(NULL))
  } 
  roc(response = df$response, predictor = df$predictor, direction = "<")
}

#' This function calculates mathhews correlations coefficient
#'
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column respons 
#' @return tidy data frame with mathhew's correlation coefficient
calc_mcc = function(df) {
  p = df %>%
    filter(response == 1) %>%
    nrow()
  
  n = df %>%
    filter(response == 0) %>%
    nrow()
  
  r = roc(response = df$response, predictor = df$predictor, direction = "<")
  res = tibble(sensitivity = r$sensitivities,
               specificity = r$specificities,
               p = p,
               n = n) %>%
    arrange(sensitivity, desc(specificity)) %>%
    mutate(tp = sensitivity * p,
           fn = p - tp,
           tn = specificity * n,
           fp = n - tn) %>%
    mutate(denom = as.double(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn),
           denom = case_when(denom == 0 ~ 1,
                            TRUE ~ denom),
           mcc = ((tp * tn) - (fp * fn))/sqrt(denom))
  return(res)
}

#' This function calculates balanced accurary
#'
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column respons 
#' @return tidy data frame with balanced accurary
calc_bac = function(df) {
  p = df %>%
    filter(response == 1) %>%
    nrow()
  
  n = df %>%
    filter(response == 0) %>%
    nrow()
  
  r = roc(response = df$response, predictor = df$predictor, direction = "<")
  res = tibble(sensitivity = r$sensitivities,
               specificity = r$specificities,
               p = p,
               n = n)  %>%
    arrange(sensitivity, desc(specificity)) %>%
    mutate(tp = sensitivity * p,
           tn = specificity * n) %>%
    mutate(bac = 0.5*(tp/p + tn/n))
    
  return(res)
  
}
