# Cross validate shapley importance values

























#' Construct a matrix showing the effect of a model for different cutoff thresholds
#'
#' @param dp data.frame containing all information of the samples
#' @param model and lrm or XGBoost model 
#' @param day integer specifying untill which day variables were collected. This
#' is also the day at which should be decided whether a bed should be emptied.
#' @param df altered data.frame used for training and testing the model which
#' contains all the sample used for testing but not all variables. Must be
#' supplied together with recipe
#' @param test_rec preprocessing recipe. Must be supplied together with df.
#'
#' @return
#' @export
#'
#' @examples
construct_decision_matrix <- function(dp, model, df = NULL, test_rec = NULL,
                                      day = 7){
  
  # If a df has been supplied, it is possible that dp contains different
  # samples than the data the model will be tested on, therefore we here
  # ensure that only samples present in the testing data are retained, and in 
  # an identical order
  if(!is.null(df)){
    dp <- dp[dp$Record.Id %in% row.names(df),]
    df <- df[row.names(df) %in% dp$Record.Id,]
    dp <- arrange(dp, Record.Id) # order samples
    if(!all(dp$Record.Id == row.names(df)))
      stop("samples in dp and df are not ordered identically")
  }
  
  ## Afkappunten bepalen
  if(class(model)[1] == "lrm"){
    kans <- predict(model, type = "fitted")
  } else {
    kans <- 1 - predict_classprob.model_fit(model_xg, bake(prep(test_rec),
                                                           new_data = df))[[1]]
  }

  hist(kans)

  
  afkap <- seq(25, 95, 2.5)
  sens  <- rep(NA, length(afkap))
  spec  <- rep(NA, length(afkap))
  ppv   <- rep(NA, length(afkap))
  npv   <- rep(NA, length(afkap))
  abst  <- rep(NA, length(afkap))
  misc  <- rep(NA, length(afkap))
  ligd  <- rep(NA, length(afkap))
  ligp  <- rep(NA, length(afkap))
  
  # Create data.frame containing true and predicted classed
  df_eval <- data.frame(factor(dp$event, levels = c(0,1), labels = c(0,1)), kans)
  df_eval$pred <- factor(ifelse(kans < 0.5, 0, 1), levels = c(0,1))
  names(df_eval) <- c("truth", "prob", "pred")
  
  for (i in 1:length(afkap)){
    
    df_eval$pred  <- factor(ifelse(kans*100 < afkap[i], 0, 1), levels = c(0,1))
    sens[i] <- round(yardstick::sensitivity(df_eval, "truth", "pred")$.estimate, 2)
    spec[i] <- round(yardstick::specificity(df_eval, "truth", "pred")$.estimate, 2)
    ppv[i]  <- round(yardstick::ppv(df_eval, "truth", "pred")$.estimate, 2)
    npv[i]  <- round(yardstick::npv(df_eval, "truth", "pred")$.estimate, 2)
    abst[i] <- sum(df_eval$pred  == 1, na.rm = TRUE) # How many people are predicted to die
    misc[i] <- table(df_eval$pred, df_eval$truth)[2, 1] # How many people are incorrectly predicted to die
    ligd[i] <- sum(ifelse(dp$ICU_LoS[df_eval$pred == 1] - day < 1, 0, # How many bed days have been saved
                          dp$ICU_LoS[df_eval$pred == 1] - day), na.rm = TRUE) # How much percentage of total bed days has been saved
    ligp[i] <- round(ligd[i]/sum(dp$ICU_LoS, na.rm = TRUE)*100, 1)
    
  }
  
  ## Afkap in (%), testkenmerken, aantal abstineren (let op: aantal al voor 7 dagen
  ## van IC af, maar dat weten we bij aanvang niet), aantal onjuist abstineren,
  ## gespaarde ligdagen, gespaarde ligdagen (%).
  decision_matrix <- data.frame(afkap, sens, spec, ppv, npv, abst, misc, ligd, ligp)

  return(decision_matrix)
}





















sofa_rec <- prep(test_rec)
df_low <- bake(sofa_rec, df[df$day1 < 8,])
df_high <- bake(sofa_rec, df[df$day1 > 12,])
df_mid <- bake(sofa_rec, df[df$day >= 8 & df$day1 <=12,])
df_high$y <- factor(df_high$y, levels = c(0,1), labels = c("Alive", "Death"))
df_mid$y <- factor(df_mid$y, levels = c(0,1), labels = c("Alive", "Death"))
df_low$y <- factor(df_low$y, levels = c(0,1), labels = c("Alive", "Death"))

shap_low <- plot_shaply_summary(df_low, mod = model_xg, save = F,
                                  n_features = 5)
shap_low$plt <- shap_low$plt + ggtitle("Shapley summary plot - SOFA < 8 on day 1")

shap_mid <- plot_shaply_summary(df_mid, mod = model_xg, save = F,
                                n_features = 5)
shap_mid$plt <- shap_mid$plt + ggtitle("Shapley summary plot - SOFA inbetween 8 & 12 on day 1")

shap_high <- plot_shaply_summary(df_high, mod = model_xg, save = F,
                                 n_features = 5)
shap_high$plt <- shap_high$plt + ggtitle("Shapley summary plot - SOFA > 12 on day 1")

shap_low

shap_mid

shap_high

