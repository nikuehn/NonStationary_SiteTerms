# different priors and ref_d
remove(list=ls())

# load necessay packages
library(INLA)
library(inlabru)

# define function o calculate total scpo (not used in paper)
slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}

# set working directory
setwd("/Users/nico/GROUNDMOTION/PROJECTS/NONERGODIC/DeltaS_NS/Git/NonStationary_SiteTerms")
dir_data <- "DATA"

# load mesh
load(file.path(dir_data, 'mesh.Rdata'))
# load distances to basin boundary for each mesh node
mesh_basin <- read.csv(file.path(dir_data, 'mesh_loc_basin.csv'))
D_basin_mesh <- mesh_basin$D_basin

# define reference distance for the non-stationary models
ref_d <- 2

# prior distributions for the standard deviations
prior_prec_tau    <- list(prec = list(prior = 'pc.prec', param = c(0.5, 0.01)))
prior_prec_phiS2S <- list(prec = list(prior = 'pc.prec', param = c(0.5, 0.01)))
prior_prec_phiSS  <- list(prec = list(prior = "loggamma", param = c(2, 0.5)))

# Prior on the fixed effects
prior.fixed <- list(mean.intercept = 0, prec.intercept = 5)

per <- 0.01
print(paste0("T = ",per))
if(per == 0.01) {
  data <- read.csv(file.path(dir_data, 'data_sitecorrected_PGA_socal_basin_utmeq.csv'))
} else {
  data <- read.csv(file.path(dir_data, sprintf('data_sitecorrected_T%.3fs_socal_basin_utmeq.csv',per)))
}
print(dim(data))


#number of data
n_stat <- length(unique(data$STATID))
n_eq <- length(unique(data$eqid))
n_rec <- length(data$recID)
n_basin <- length(unique(data$Basin_ID))

#identify station and eq ids
statid <- data$STATID
eqid <- data$eqid
basinid <- data$Basin_ID

#SET STATION INDICES
stat  <- array(0,n_rec)
statid_uni <- unique(statid)
for (i in 1:n_stat){
  indi <- ifelse(statid %in% statid_uni[i],i,0)
  stat <- stat + indi
}

#SET EVENT INDICES
eq  <- array(0,n_rec)
eqid_uni <- unique(eqid)
for (i in 1:n_eq){
  indi <- ifelse(eqid %in% eqid_uni[i],i,0)
  eq <- eq + indi
}

# other data/predictors
resid <- data$lnY - data$Median

# set basin distances
D_b_mesh <- D_basin_mesh
D_b_mesh[D_b_mesh > ref_d] <- ref_d

D_b <- data$D_basin
D_b[D_b > ref_d] <- ref_d

# data for regression
data_reg <- data.frame(Y = resid,
                       eq=eq,
                       stat=stat,
                       basin = basinid,
                       D_b = D_b,
                       intercept = 1
)


co_eq <- as.matrix(cbind(data$X_eq, data$Y_eq))
co_stat <- as.matrix(cbind(data$X_stat, data$Y_stat))
names(co_eq) <- c("X","Y")
names(co_stat) <- c("X","Y")

# ergodic model (event term/station term)
fit_inla <- inla(Y ~ 1  + f(eq, model = "iid", hyper = prior_prec_tau) + 
                   f(stat, model = "iid",hyper = prior_prec_phiS2S), 
                 data = data_reg,
                 family="gaussian",
                 control.fixed = prior.fixed,
                 control.family = list(hyper = list(prec = prior_prec_phiSS)),
                 control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                 quantiles = c(0.05, 0.5, 0.95)
)


#-------------------------------------------------------------------------------
# Spatial Models
#-------------------------------------------------------------------------------
# spde priors
spde_stat <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 0.100) = 0.99
  prior.range = c(100, 0.99),
  # P(sigma > .5) = 0.01
  prior.sigma = c(.5, 0.01)) 

# define projection matrix
A_stat   <- inla.spde.make.A(mesh, loc = co_stat)
idx.stat <- inla.spde.make.index("idx.stat",spde_stat$n.spde)

form_spatial_stat <- y ~ 0 + intercept + 
  f(eq, model = "iid", hyper = prior_prec_tau) +
  f(stat, model = "iid",hyper = prior_prec_phiS2S) +
  f(idx.stat, model = spde_stat)

# build stack
stk_spatial_stat <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_stat, 1),
  effects = list(idx.stat = idx.stat,
                 data_reg),
  tag = 'model_stk_spatial_stat')

fit_inla_spatial_stat <- inla(form_spatial_stat, 
                              data = inla.stack.data(stk_spatial_stat),
                              family="gaussian",
                              control.fixed = prior.fixed,
                              control.family = list(hyper = list(prec = prior_prec_phiSS)),
                              control.predictor = list(A = inla.stack.A(stk_spatial_stat)),
                              control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                              #control.inla = list(int.strategy='eb', strategy="gaussian"),
                              verbose=FALSE,
                              quantiles = c(0.05, 0.5, 0.95)
)
# get mode values of estimated parameters in internal scale
# to use as starting values for the non-stationary models
theta_inla_stat <- fit_inla_spatial_stat$mode$theta

#######################################
# group model with stations
# separate correlation based on basin id

A_statgr     <- inla.spde.make.A(mesh, loc = co_stat,
                                 group = basinid, n.group = n_basin)
idx.statgr   <- inla.spde.make.index("idx.statgr",spde_stat$n.spde,
                                     n.group = n_basin)

stk_spatial_statgr <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_statgr, 1),
  effects = list(idx.statgr = idx.statgr,
                 data_reg),
  tag = 'model_stk_spatial_statgr')

form_spatial_statgr <- y ~ 0 + intercept + 
  f(eq, model = "iid", hyper = prior_prec_tau) +
  f(stat, model = "iid",hyper = prior_prec_phiS2S) +
  f(idx.statgr, model = spde_stat, group = idx.statgr.group,
    control.group = list(model = "iid"))

theta <- theta_inla_stat

fit_inla_spatial_statgr <- inla(form_spatial_statgr, 
                                data = inla.stack.data(stk_spatial_statgr),
                                family="gaussian",
                                control.fixed = prior.fixed,
                                control.family = list(hyper = list(prec = prior_prec_phiSS)),
                                control.predictor = list(A = inla.stack.A(stk_spatial_statgr)),
                                control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                                #control.inla = list(int.strategy='eb', strategy="gaussian"),
                                verbose=FALSE,
                                control.mode=list(theta = theta, restart=TRUE),
                                quantiles = c(0.05, 0.5, 0.95)
)


#######################################
#### non-stationary inla model
# correlation depends on distance to basin

nu <- 1 
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0

# SPDE model
spde_stat_ns <- inla.spde2.matern(mesh, 
                                  B.tau = cbind(logtau0, -1, nu, nu * D_b_mesh), 
                                  B.kappa = cbind(logkappa0, 0, -1, -1 * D_b_mesh),
                                  theta.prior.mean = rep(0, 3), 
                                  theta.prior.prec = rep(1, 3)) 


idx.stat_ns <- inla.spde.make.index("idx.stat_ns",spde_stat_ns$n.spde)

form_spatial_stat_ns <- y ~ 0 + intercept + 
  f(eq, model = "iid", hyper = prior_prec_tau) +
  f(stat, model = "iid",hyper = prior_prec_phiS2S) +
  f(idx.stat_ns, model = spde_stat_ns)

# build stack
stk_spatial_stat_ns <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_stat, 1),
  effects = list(idx.stat_ns = idx.stat_ns,
                 data_reg),
  tag = 'model_stk_spatial_stat_ns')

# starting values
# staing values corespods 5 times  increase in range from zero to ref_d
theta <- c(theta_inla_stat[c(1,2,3,5)], theta_inla_stat[4] - log(5), log(5)/ref_d)

fit_inla_spatial_stat_ns <- inla(form_spatial_stat_ns, 
                                 data = inla.stack.data(stk_spatial_stat_ns),
                                 family="gaussian",
                                 control.fixed = prior.fixed,
                                 control.family = list(hyper = list(prec = prior_prec_phiSS)),
                                 control.predictor = list(A = inla.stack.A(stk_spatial_stat_ns)),
                                 control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                                 #control.inla = list(int.strategy='eb', strategy="gaussian"),
                                 verbose=FALSE,
                                 control.mode=list(theta = theta, restart=TRUE),
                                 quantiles = c(0.05, 0.5, 0.95)
)

#######################################
#### non-stationary inla model
# correlation depends on distance to basin
# separate correlation based on basin id

idx.statgr_ns   <- inla.spde.make.index("idx.statgr_ns",spde_stat_ns$n.spde,
                                        n.group = n_basin)

stk_spatial_statgr_ns <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_statgr, 1),
  effects = list(idx.statgr_ns = idx.statgr_ns,
                 data_reg),
  tag = 'model_stk_spatial_statgr_ns')

form_spatial_statgr_ns <- y ~ 0 + intercept + 
  f(eq, model = "iid", hyper = prior_prec_tau) +
  f(stat, model = "iid",hyper = prior_prec_phiS2S) +
  f(idx.statgr_ns, model = spde_stat_ns, group = idx.statgr_ns.group,
    control.group = list(model = "iid"))

fit_inla_spatial_statgr_ns <- inla(form_spatial_statgr_ns, 
                                   data = inla.stack.data(stk_spatial_statgr_ns),
                                   family="gaussian",
                                   control.fixed = prior.fixed,
                                   control.family = list(hyper = list(prec = prior_prec_phiSS)),
                                   control.predictor = list(A = inla.stack.A(stk_spatial_statgr_ns)),
                                   control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                                   #control.inla = list(int.strategy='eb', strategy="gaussian"),
                                   verbose=FALSE,
                                   control.mode=list(theta = theta, restart=TRUE),
                                   quantiles = c(0.05, 0.5, 0.95)
)



###################
c(fit_inla$waic$waic,
  fit_inla_spatial_stat$waic$waic,
  fit_inla_spatial_statgr$waic$waic,
  fit_inla_spatial_stat_ns$waic$waic,
  fit_inla_spatial_statgr_ns$waic$waic
)

summary(fit_inla)
summary(fit_inla_spatial_stat)
summary(fit_inla_spatial_statgr)
summary(fit_inla_spatial_stat_ns)
summary(fit_inla_spatial_statgr_ns)



