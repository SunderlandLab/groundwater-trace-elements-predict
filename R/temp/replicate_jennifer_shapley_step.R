final_xg_list <- extract_fit_engine(final_model_copy)
# Initialize output lists
df_pred_list <- list()
factorcols_list <- list()

# Loop over imputed datasets

df <- rbind(df_train_copy, df_test_copy)

# Process predictors using recipe
df_pred <- bake(
  prep(train_recipes_copy), # bake with the reference recipe to ensure consistency
  has_role('predictor'),
  new_data = df,
  composition = 'matrix'
)

# Extract factor names
factors <- df %>%
  dplyr::select(-c('location.id', 'data.source')) %>%
  dplyr::select_if(is.factor) %>%
  colnames()

# Extract factor-related columns
factorcols <- list()
for (factor in factors) {
  cols <- as.data.frame(df_pred) %>%
    dplyr::select(starts_with(factor)) %>%
    colnames()
  factorcols[[factor]] <- cols
}

# Align df_pred with model features
remove_features <- c(setdiff(colnames(df_pred), final_xg_list$feature_names))
print(paste0("features to remove from df_pred are ", remove_features))
remove_index <- which(colnames(df_pred) %in% remove_features)
if (length(remove_index) > 0) {
  df_pred <- df_pred[, -remove_index]
}

# Save results
df_pred_list <- df_pred
factorcols_list <- factorcols

#### 2a. Calculate xgboost shapley values --------------------------------------------------------------------------------------------------------


df <- rbind(df_train_copy, df_test_copy)
df_pred <- df_pred_list
factorcols <- factorcols_list

shap_list <- shapviz(
  object = final_xg_list,
  X_pred = df_pred,
  X = df,
  collapse = factorcols
)


# 1. Extract SHAP matrices
shap_matrices <- shap_list$S
# 2. Pool SHAP values by averaging across imputations
# Find common column names
common_cols <- Reduce(intersect, lapply(shap_matrices, colnames))

# # Subset each matrix to only those common columns
# shap_matrices_common <- lapply(shap_matrices, function(mat) mat[, common_cols, drop = FALSE])
# 
# # Average them
# pooled_S <- Reduce(`+`, shap_matrices_common) / length(shap_matrices_common)
# 3. Pool baseline values
baseline_values <- shap_list$baseline
pooled_baseline <- mean(baseline_values)
# 4. Use the feature matrix from the first object
X_ref <- shap_list$X
# 5. Reconstruct a new shapviz object
shap_pooled <- shapviz(
  shap_matrices,
  X = X_ref,
  baseline_value = pooled_baseline,
  model_class = "xgboost"
)
varImp.plot <- sv_importance(shap_pooled, kind = "both", show_numbers = TRUE) +
  scale_y_discrete(limits=rev(calc_meanabs_shap(shap_pooled))) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size=12))
varImp.plot
