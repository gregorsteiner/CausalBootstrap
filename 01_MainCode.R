
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
  compensation <- sapply(regmatches(desiredcompensation, 
                                    gregexpr("?[0-9]+[.]?[0-9]*",
                                             desiredcompensation)),
                         function(x) mean(as.numeric(x)))
  hours <- sapply(regmatches(hourswork,
                             gregexpr("?[0-9]+[.]?[0-9]*", hourswork)),
                  function(x) mean(as.numeric(x)))
  
  
  
})



######### Tables & Plots #########

# summary table
vtable::sumtable(dat[, c("treatment", "male", "single",
                         "class", "compensation", "hours")],
                 out = "latex", file = "Tables/SummaryTable.tex",
                 anchor = "SumStats", title = "Summary Statistics",
                 fit.page = NA)


# compute compensation and hours by class and estimate the ATE
Results <- do.call(cbind, lapply(dat[, c("compensation", "hours")],
                                 function(x){
  mat <- tapply(x, list(dat$class, dat$treatment), mean, na.rm = TRUE)
  mat <- cbind(mat, "Difference" = mat[, "Public"]- mat[, "Private"])
  mat
}))

# export as tex table
writeLines(knitr::kable(Results, format = "latex", digits = 2), "Tables/Means.tex")




######### Bootstrap #########


# population sizes
p.single <- mean(dat$maritalstatus == 0, na.rm = TRUE)
pop.sizes <- c(c(1-p.single, p.single) * 21000000,
               c(1-p.single, p.single) * 21500000)

# bootstrap CIs for compensation and hours worked
set.seed(1) # set seed such that results are reproducible

# loop over classes
Boot.CIs <- t(mapply(function(class, pop.size){
  # select class
  data <- dat[dat$class == class, ]
  
  sapply(c("compensation", "hours"), function(out){
    # compute CIs
    boot <- causal_boot(data, dep.var = out,
                        treatment = "treatnum", N = pop.size)
    
    # combine table with CIs and point estimate
    res <- as.numeric(c("Lower Bound" = boot[["Confidence Interval"]]["Lower"],
                        "Upper Bound" = boot[["Confidence Interval"]]["Upper"]))
    return(res)
  })
  
}, levels(dat$class), pop.sizes))

colnames(Boot.CIs) <- rep(c("Lower Bound", "Upper Bound"), 2)

# convert to tex table
boot.res <- knitr::kable(Boot.CIs, format = "latex", digits = 2)

# add multicolummns as string
boot.res <- sub("\\hline",
                "\\hline & \\multicolumn{2}{c|}{Compensation} & \\multicolumn{2}{|c}{Hours} \\\\ \\hline",
                boot.res, fixed = TRUE)


# save table
writeLines(boot.res, "Tables/BootResults.tex")





