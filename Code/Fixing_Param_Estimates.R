#######################
# Fixing Functions & Re-doing Analyses 
# for Updated Parameter Estimates
#
# by Carly Fagnant
#######################

# to be copy & pasted into or edited together with Testing_on_Real_Data.R


# Model 1 -----------------------------------------------------------------

# Function to get parameter estimate matrix from CAR model fits, given...
#   ln.scale.fit: the CAR fit object for ln.scale
#   shape.fit: the CAR fit object for shape
#   rate.fit: the CAR fit object for rate
get_par_car_fix <- function(ln.scale.fit, shape.fit, rate.fit){
  car_fit <- rbind(exp(ln.scale.fit$fit$coef) * (1 + (diag(ln.scale.fit$fit$imat * ln.scale.fit$fit$s2)/2)),   # updated with transformation/approximation
                   # exp(ln.scale.fit$fit$coef),  # old estimate
                   shape.fit$fit$coef,
                   rate.fit$fit$coef)
  rownames(car_fit) <- c("scale", "shape", "rate")
  return(car_fit)
}

# summary(car_test_ln.scale)$fit$coef
# car_test_ln.scale$fit$coef
# exp(car_test_ln.scale$fit$coef)
# diag(car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2)
# 
# exp(car_test_ln.scale$fit$coef) * (1 + (diag(car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2)/2))


# Don't use get_se function! Use make_varcov_mats_car and get_se_from_varcov (from Testing_on_Simulated_Data.R)

car_fit_fix <- get_par_car_fix(car_test_ln.scale, car_test_shape, car_test_rate) 
car_fit
car_fit <- get_par_car_fix(car_test_ln.scale, car_test_shape, car_test_rate) 
get_se_from_varcov(make_varcov_mats_car(car_test_ln.scale, car_test_shape))

par_22 <- get_par_car_fix(car_ln.scale_22, car_shape_22, car_rate_22)
get_se_from_varcov(make_varcov_mats_car(car_ln.scale_22, car_shape_22))

par_52 <- get_par_car_fix(car_ln.scale_52, car_shape_52, car_rate_52)
get_se_from_varcov(make_varcov_mats_car(car_ln.scale_52, car_shape_52))



# Model 2 -----------------------------------------------------------------

# Function to get parameter estimate matrix from kriging model fits, given...
#   krig_regs: the co-krige fit object for ln.scale and shape
#   krig_regs_rate: the krige fit object for rate
get_par_krig_fix <- function(krig_regs, krig_regs_rate){
  krig_fit <- rbind(exp(krig_regs$ln.scale.pred)*(1 + (krig_regs$ln.scale.var/2)),    # updated with transformation/approximation
                    # exp(krig_regs$ln.scale.pred),   # old estimate
                    krig_regs$shape.pred,
                    krig_regs_rate$var1.pred)
  rownames(krig_fit) <- c("scale", "shape", "rate")
  colnames(krig_fit) <- c("Reg1", "Reg2", "Reg3")
  return(krig_fit)
}

krig_fit_fix <- get_par_krig_fix(krig_regs, krig_regs_rate) 
krig_fit
krig_fit <- get_par_krig_fix(krig_regs, krig_regs_rate) 
get_se_from_varcov(make_varcov_mats_krig(krig_regs$ln.scale.var, krig_regs$ln.scale.pred, krig_regs$shape.var, krig_regs$cov.ln.scale.shape))

par_krig_22 <- get_par_krig_fix(krig_regs_22, krig_regs_rate_22)
get_se_from_varcov(make_varcov_mats_krig(krig_regs_22$ln.scale.var, krig_regs_22$ln.scale.pred, krig_regs_22$shape.var, krig_regs_22$cov.ln.scale.shape))

par_krig_52 <- get_par_krig_fix(krig_regs_52, krig_regs_rate_52)
get_se_from_varcov(make_varcov_mats_krig(krig_regs_52$ln.scale.var, krig_regs_52$ln.scale.pred, krig_regs_52$shape.var, krig_regs_52$cov.ln.scale.shape))


