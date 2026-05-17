####################################################################################################
# File name: 5b_Rebuild_Models.R
# Purpose: Convert xgboost model objects to be compatible with the current xgboost version.
#          The $raw bytes in saved boosters use xgboost's "serialize" (full-state) format, not
#          the model format. This script uses XGBoosterUnserializeFromBuffer_R to load them,
#          re-saves each booster to a temp file in the new format, then patches the workflows
#          and re-saves each .RData package so script 6 can run without retraining.
# Usage (R 4.6, run from project root):
#   Rscript R/5b_Rebuild_Models.R
####################################################################################################

library(here)
setwd(here::here("data"))

metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")

rebuild_booster <- function(booster) {
  raw <- booster$raw
  if (is.null(raw) || length(raw) == 0) {
    stop("booster$raw is empty — cannot reconstruct. Re-running 5_Regression_Model.R is required.")
  }
  # Use the serialize/unserialize path (XGBoosterUnserializeFromBuffer_R), not the model-load
  # path (XGBoosterLoadModelFromRaw_R), because $raw stores the full training state.
  bst    <- .Call(xgboost:::XGBoosterCreate_R, list())
  handle <- xgboost:::xgb.get.handle(bst)
  suppressWarnings(.Call(xgboost:::XGBoosterUnserializeFromBuffer_R, handle, raw))
  # Save to a temp file in the current format, then reload for a clean xgb.Booster object
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))
  xgboost::xgb.save(bst, tmp)
  xgboost::xgb.load(tmp)
}

fix_workflow <- function(wf) {
  wf$fit$fit$fit <- rebuild_booster(wf$fit$fit$fit)
  wf
}

#### Regression model packages (all metals) ####
for (metal.code in metal.codes) {
  cat("Rebuilding regression model:", metal.code, "\n")
  load(paste0("R_Output/", metal.code, "_ModelPackage_2step.RData"))

  final_model      <- lapply(final_model,      fix_workflow)
  final_full_model <- lapply(final_full_model, fix_workflow)

  save(final_model, final_full_model, df_test, df_train,
       train_recipes, model_workflows,
       file = paste0("R_Output/", metal.code, "_ModelPackage_2step.RData"))
  cat("  Saved\n")
}

#### Detection model (Cd only) ####
cat("Rebuilding detection model: Cd\n")
load("R_Output/Cd_DetectionModel.RData")

detection_model <- lapply(detection_model, fix_workflow)

save(detection_model, detection_df_train, detection_df_test,
     detection_model_workflows,
     file = "R_Output/Cd_DetectionModel.RData")
cat("  Saved\n")

cat("\nAll models rebuilt. You can now run 6_Model_Eval_Reg.R.\n")
