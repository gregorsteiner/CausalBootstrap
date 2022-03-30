
######### Preliminaries ######### 

# source functions
source("00_Functions.R")

######### Load and prepare data ######### 

# read data
dat <- as.data.frame(haven::read_dta("Replication Data/3_mainexp_responses.dta"))

# delete observation 167 due to encoding problem (Invalid byte sequence)
dat <- dat[-167, ]


# data transformations
dat <- within(dat, {
  # Create class factor
  class <- factor(dplyr::case_when(
    male == 0 & maritalstatus == 0 ~ "Female & Single",
    male == 0 & maritalstatus != 0 ~ "Female & Non-single",
    male == 1 & maritalstatus == 0 ~ "Male & Single",
    male == 1 & maritalstatus != 0 ~ "Male & Non-single"
  ))
  
  # treatment as factor
  treatment <- factor(ifelse(treatment == "A", "Private", "Public"))
  
  # add numerical treatment
  treatnum <- as.numeric(treatment == "Public")
  
  # create continous compensation and hours worked variables
  compensation <- sapply(regmatches(desiredcompensation, gregexpr("?[0-9]+[.]?[0-9]*", desiredcompensation)),
                         function(x) mean(as.numeric(x)))
  hours <- sapply(regmatches(hourswork, gregexpr("?[0-9]+[.]?[0-9]*", hourswork)),
                  function(x) mean(as.numeric(x)))
  
  
  
})



######### Tables & Plots #########


# compute compensation and hours by class and estimate the ATE
Results <- do.call(cbind, lapply(dat[, c("compensation", "hours")], function(x){
  mat <- tapply(x, list(dat$class, dat$treatment), mean, na.rm = TRUE)
  mat <- cbind(mat, "Difference" = mat[, "Public"]- mat[, "Private"])
  mat
}))

# export as tex table
knitr::kable(Results, format = "latex", digits = 2)





######### Bootstrap #########



# focus on single females
dat.female.single <- dat[dat$class == "Female & Single", ]

# population size
N <- round(mean(dat$maritalstatus == 0, na.rm = TRUE) * 21000000) # population size

# bootstrap CIs for compensation and hours worked
set.seed(1) # set seed such that results are reproducible
Boot.CIs <- do.call(rbind, lapply(c("compensation", "hours"), function(out){
  # compute CIs
  boot <- causal_boot(dat.female.single, dep.var = out, treatment = "treatnum", N = N)
  
  # combine table with CIs and point estimate
  res <- as.numeric(c("Lower Bound" = boot[["Confidence Interval"]]["Lower"],
           "Point Estimate" = boot["ATE Estimator"],
           "Upper Bound" = boot[["Confidence Interval"]]["Upper"]))
  return(res)
}))

colnames(Boot.CIs) <- c("Lower Bound", "Point Estimate", "Upper Bound")
rownames(Boot.CIs) <- c("Compensation", "Hours")

# export as tex table
knitr::kable(as.matrix(Boot.CIs), format = "latex", digits = 2)





