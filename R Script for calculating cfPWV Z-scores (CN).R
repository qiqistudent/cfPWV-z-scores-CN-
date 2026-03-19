# ============================================================
# R script for calculating cfPWV z-scores
# ============================================================

library(gamlss)

# 1. Load the .rds file containing the male and female models
# Make sure the file is in the working directory
models_list <- readRDS("cfPWV_Reference_Models.rds")


# 2. Define a function to calculate z-scores by sex
get_zscore_by_sex <- function(models, new_data) {
  # new_data must contain the following columns:
  # 'age', 'ht', 'cfPWV', and 'sex'
  # sex can be coded as 1/2 or "Male"/"Female"
  
  required_cols <- c("age", "ht", "cfPWV", "sex")
  missing_cols <- setdiff(required_cols, names(new_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Initialize output columns
  new_data$z_score <- NA_real_
  new_data$percentile <- NA_real_
  new_data$pred_median <- NA_real_
  
  # Standardize sex coding
  sex_chr <- tolower(as.character(new_data$sex))
  
  # --- Process male data ---
  idx_male <- which(sex_chr %in% c("1", "male", "m"))
  
  if (length(idx_male) > 0) {
    data_m <- new_data[idx_male, , drop = FALSE]
    male_model <- models$Male$model
    male_ref <- models$Male$ref_data
    
    # Warn if extrapolation is required
    if (any(data_m$age < min(male_ref$age) | data_m$age > max(male_ref$age) |
            data_m$ht  < min(male_ref$ht)  | data_m$ht  > max(male_ref$ht))) {
      warning("Some male observations are outside the reference age/height range; results are extrapolated.")
    }
    
    # Predict distribution parameters
    preds_m <- predictAll(
      male_model,
      newdata = data_m,
      type = "response",
      data = male_ref
    )
    
    # Calculate cumulative probabilities under the fitted BCCGo distribution
    probs_m <- pBCCGo(
      q = data_m$cfPWV,
      mu = preds_m$mu,
      sigma = preds_m$sigma,
      nu = preds_m$nu
    )
    
    # Avoid infinite z-scores
    probs_m <- pmin(pmax(probs_m, 1e-10), 1 - 1e-10)
    
    # Fill in results
    new_data$z_score[idx_male] <- qnorm(probs_m)
    new_data$percentile[idx_male] <- probs_m * 100
    new_data$pred_median[idx_male] <- qBCCGo(
      0.5,
      mu = preds_m$mu,
      sigma = preds_m$sigma,
      nu = preds_m$nu
    )
  }
  
  # --- Process female data ---
  idx_female <- which(sex_chr %in% c("2", "female", "f"))
  
  if (length(idx_female) > 0) {
    data_f <- new_data[idx_female, , drop = FALSE]
    female_model <- models$Female$model
    female_ref <- models$Female$ref_data
    
    # Warn if extrapolation is required
    if (any(data_f$age < min(female_ref$age) | data_f$age > max(female_ref$age) |
            data_f$ht  < min(female_ref$ht)  | data_f$ht  > max(female_ref$ht))) {
      warning("Some female observations are outside the reference age/height range; results are extrapolated.")
    }
    
    # Predict distribution parameters
    preds_f <- predictAll(
      female_model,
      newdata = data_f,
      type = "response",
      data = female_ref
    )
    
    # Calculate cumulative probabilities under the fitted BCCGo distribution
    probs_f <- pBCCGo(
      q = data_f$cfPWV,
      mu = preds_f$mu,
      sigma = preds_f$sigma,
      nu = preds_f$nu
    )
    
    # Avoid infinite z-scores
    probs_f <- pmin(pmax(probs_f, 1e-10), 1 - 1e-10)
    
    # Fill in results
    new_data$z_score[idx_female] <- qnorm(probs_f)
    new_data$percentile[idx_female] <- probs_f * 100
    new_data$pred_median[idx_female] <- qBCCGo(
      0.5,
      mu = preds_f$mu,
      sigma = preds_f$sigma,
      nu = preds_f$nu
    )
  }
  
  # Warn if some rows have unrecognized sex coding
  unknown_idx <- which(!(sex_chr %in% c("1", "2", "male", "female", "m", "f")))
  if (length(unknown_idx) > 0) {
    warning("Some rows have unrecognized sex coding and were left as NA.")
  }
  
  return(new_data)
}


# ============================================================
# 3. Example
# ============================================================

# Example input data
# Here, 1 = Male and 2 = Female
test_patients <- data.frame(
  id = 1:4,
  sex = c(1, 2, 1, 2),
  age = c(10, 10, 14, 14),
  ht = c(140, 140, 165, 158),
  cfPWV = c(5.5, 5.5, 6.0, 6.0)
)

# Run the calculation
results <- get_zscore_by_sex(models_list, test_patients)

# View the results
print(results)