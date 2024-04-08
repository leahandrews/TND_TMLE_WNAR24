
## data_func creates simulated TND sample dataset
## inputs
# a_betas: vector of beta (log OR) coefficients with model coefficient 
#   names to generate exposure/vaccination variable A
# a_formula: string of formula to generate exposure/vaccination 
#    variable A must include a_betas variables
#    e.g., "~AGEGR1+SEX +RACE6c+ETHNICc+RISKGR1+REGION+CALTIME"
# a_y_beta: integer of exposure of interest beta (log OR) coefficient
#.  for generating SARS-infection variable Y. Can be unnamed
# not_a_y_betas: vector of confounder beta (log OR) coefficients 
#    with model coefficient names to generate SARS-infection variable Y, 
#    excluding exposure of interest coefficient
# y_formula: string of formula to generate SARS infection
#    variable Y must include y_betas variables
#.   e.g. "~A + AGE+SEX +RACE6c+ETHNICc+RISKGR1+REGION+CALTIME+CALTIMEs2+CALTIMEs3
# a_w_beta: integer of exposure of interest beta (log OR) coefficient
#.  for generating non-SARS infection variable W. Integer can be unnamed
# not_a_w_betas: vector of confounder beta (log OR) coefficients 
#    with model coefficient names to generate non-SARS infection variable W
# w_formula: string of formula to generate non-SARS infection
#    variable W must include w_betas variables
#    e.g. "~A + AGE+SEX +RACE6c+ETHNICc+RISKGR1+REGION+CALTIME+CALTIMEs2+CALTIMEs4 
# a_c_beta: integer of exposure of interest beta (log OR) coefficient
#.  for generating clinical symptom variable C. Can be unnamed
# not_a_c_betas: vector of confounder beta (log OR) coefficients
#     with model coefficient names to generate clinical symptom variable C
# c_formula: string of formula to generate clinical symptom 
#    variable C must include c_betas variables
#     e.g., "~Y+W+ A+Y*AGEGR1+ Y*SEX +Y*RISKGR1 +W*AGEGR1 + W*SEX+ W*RISKGR1"
# missing_rate_vec: vector of desired MCAR missing rates
#.   2 allowed currently
# mar20_d_betas: named vector of beta values for MAR beta parameters for first missing rate
# mar50_d_betas: named vector of beta values for MAR beta parameters for second missing rate
# mar_d_formula: string of formula to generate MAR missing vaccination
#.    e.g., "~ RISKGR1"
# pop_size: integer of desired population sample size
# sample_size: integer of desired TND sample size
# rct.df: original RCT data without missing data 
#   to get covariate distribution from. 
#   Must have covariates mentioned previously in the dataset
data_func <- function(a_betas,a_formula,
                       a_y_beta, not_a_y_betas, y_formula,
                       a_w_beta, not_a_w_betas,w_formula,
                       a_c_beta, not_a_c_betas,c_formula,
                       missing_rate_vec=c(0.2,0.5),
                       mar20_d_betas, mar50_d_betas,
                       mar_d_formula, 
                       pop_size,sample_size,rct.df){
  
  # Resample full RCT data with replacement to get new population
  pop_ids<- sample(1:nrow(rct.df), pop_size, replace=TRUE)
  simpop.df<- rct.df[pop_ids,]
  
  # No Missingness in A in the Population
  simpop.df$D <- 0
  
  ## Create Vaccination A
  # Design matrix to construct A
  x_a_matrix<- model.matrix(eval(parse(text=a_formula)),simpop.df)
  reorder_a_betas <- a_betas[match(colnames(x_a_matrix),names(a_betas))]
  a_betas.df<- data.frame(variable=names(reorder_a_betas),
                          a_beta=reorder_a_betas)
  # Create Vaccination Variable
  simpop.df$A <- rbinom(pop_size ,1, plogis(x_a_matrix %*% reorder_a_betas))
  
  ## Create SARS Infection Y and Non-SARS Infection W
  # Design matrix for Y
  x_y_matrix<- model.matrix(eval(parse(text=y_formula)),simpop.df)
  # Order y_betas to match design matrix
  names(a_y_beta) <- "A"
  y_betas <- c(a_y_beta, not_a_y_betas)
  reorder_y_betas <- y_betas[match(colnames(x_y_matrix),names(y_betas))]
  y_betas.df<- data.frame(variable=names(reorder_y_betas),
                          y_beta=reorder_y_betas)
  # design matrix for Non-SARS Infection W
  x_w_matrix<- model.matrix(eval(parse(text=w_formula)),simpop.df) 
  # Order w_betas to match design matrix
  names(a_w_beta) <- "A"
  w_betas <- c(a_w_beta, not_a_w_betas)
  reorder_w_betas <- w_betas[match(colnames(x_w_matrix),names(w_betas))]
  w_betas.df<- data.frame(variable=names(reorder_w_betas),
                          w_beta=reorder_w_betas)
  # Create SARS and NonSARS Variables
  simpop.df$Y <- rbinom(pop_size,1,plogis(x_y_matrix %*% reorder_y_betas))
  simpop.df$W <- rbinom(pop_size,1,plogis(x_w_matrix %*% reorder_w_betas))
  
  ## Create Symptoms Variable C
  # design matrix to construct C
  x_c_matrix<- model.matrix(eval(parse(text=c_formula)),simpop.df)
  # Order c_betas to match design matrix
  names(a_c_beta) <- "A"
  c_betas <- c(a_c_beta, not_a_c_betas)
  reorder_c_betas <- c_betas[match(colnames(x_c_matrix),names(c_betas))]
  c_betas.df<- data.frame(variable=names(reorder_c_betas),
                          c_beta=reorder_c_betas) 
  
  # Create Symptom Variable C
  simpop.df$C <- rbinom(pop_size,1,plogis(x_c_matrix %*% reorder_c_betas))
  
  # Create TND Sample
  simpop.df$S <- 0
  
  if(sum(simpop.df$C==1)< sample_size){
    sample_size <- sum(simpop.df$C==1)
    warning(paste0("Fewer symptomatic individuals than desired sample size. ",
                   sample_size," individuals sampled instead.") )}
  
  samp_ids <- sample(which(simpop.df$C==1),sample_size, replace=FALSE) 
  simpop.df[samp_ids,"S"] <- 1
  sample.df <- simpop.df[samp_ids,]
  
  ## Create Vaccination Missingness
  
  # Vaccination Missing Completely at Random
  sample.df$DMCAR20 <- rbinom(sample_size,1,missing_rate_vec[1])
  sample.df$DMCAR50 <- rbinom(sample_size,1,missing_rate_vec[2])
  
  # Vaccination Missing at Random
  mar_x_d_matrix<- model.matrix(eval(parse(text=mar_d_formula)),sample.df) 
  reorder_mar20_d_betas <- mar20_d_betas[match(colnames(mar_x_d_matrix),names(mar20_d_betas))]
  reorder_mar50_d_betas <- mar50_d_betas[match(colnames(mar_x_d_matrix),names(mar50_d_betas))]
  mar20_d_betas.df<- data.frame(variable=names(reorder_mar20_d_betas),
                                mar20_d_beta=reorder_mar20_d_betas,
                                mar50_d_beta=reorder_mar50_d_betas) 
  
  sample.df$DMAR20 <- rbinom(sample_size,1,plogis(mar_x_d_matrix %*% mar20_d_betas))
  sample.df$DMAR50 <- rbinom(sample_size,1,plogis(mar_x_d_matrix %*% mar50_d_betas))  
  
  # Collect all Data-Generating Parameters
  
  param.df <- merge(a_betas.df, y_betas.df,  all.x = TRUE, all.y =TRUE)
  param.df <- merge(param.df, w_betas.df,  all.x = TRUE, all.y = TRUE)
  param.df <- merge(param.df, c_betas.df,  all.x = TRUE, all.y = TRUE)
  param.df <- merge(param.df, mar20_d_betas.df,  all.x = TRUE, all.y = TRUE)
  
  
  return(list(Population=simpop.df,
              Sample=sample.df,
              DataParams = param.df))
}


### ordlogit_func runs an ordinary logistic regression
# and outputs list of model values
## ordlogit_tbl inputs:
# poi_var: string of exposure of interest variable name
#  eg. "A"
# outcome_var: string of outcome variable name
#  eg. "Y"
# missing_var: string of missing exposure variable name
# eg "D"
# model_formula: string of desired regression formula
#  eg. "Y~A+RISKGR1+SEX"
# df: dataset with no missing values 
#.   that contains variables in model_formula
#. e.g., sim.df
## ordlogi_tbl outputs:
# list that contains parameters and full glm model
ordlogit_func <- function(poi_var, outcome_var,model_formula,
                          missing_var="D",
                          oracle_formula = model_formula,
                          oracle_flag= F,learner="glm", df){
  
  # Remove Individuals with Missing Vaccination
  complete.df <- df %>% filter(get(missing_var)==0)
  
  # Run model using the oracle formula
  if(oracle_flag==T){
    model_formula <- oracle_formula
  }
  
  ## Run Logistic Regression Model
  mod <- glm(eval(parse(text=model_formula)), data=complete.df, family="binomial")
  glm.mod <- summary(mod)
  
  ## Calculating relevant statistics
  coef <- glm.mod$coefficients[poi_var,"Estimate"]
  se <- glm.mod$coefficients[poi_var,"Std. Error"]
  ve <- (1-exp(coef))*100
  z <-qnorm(0.975)
  veCIlow <- (1-exp(coef +z*se))*100
  veCIhigh <- (1-exp(coef - z*se))*100
  pval = glm.mod$coefficients[poi_var,"Pr(>|z|)"]
  
  return(list(StatMethod = "Ordinary", 
               Beta = coef,
               SE = se,
               VE = ve,
               VELL = veCIlow,
               VEUL = veCIhigh,
               PValue = pval,
               Reject = as.integer(pval<0.05),
               nsamp = nrow(df),
               ncomplete = nrow(complete.df),
              nVac = with(complete.df, sum(get(outcome_var)==1)),
              nCase = with(complete.df, sum(get(poi_var)==1)),
              nW = with(complete.df, sum(W==1)),
              nC = with(complete.df, sum(C==1)),
               nVacCase = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1)),
               nUnvacCase = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1)),
              nVacCaseCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1&RISKGR1=="At Risk")),
              nUnvacCaseCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1&RISKGR1=="At Risk")),
              nVacCaseNoCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1&RISKGR1=="Not At Risk")),
              nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1&RISKGR1=="Not At Risk")),
              nVacNoncaseCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==0&RISKGR1=="At Risk")),
              nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==0&RISKGR1=="At Risk")),
              nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==0&RISKGR1=="Not At Risk")),
              nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==0&RISKGR1=="Not At Risk")),
               missing_var = missing_var, 
               Formula = model_formula,
               Model = mod)) }

### semilogit_func runs an ordinary logistic regression
#   and outputs list of model values
##  inputs:
# poi_var: string of exposure of interest variable name
#  eg. "A"
# outcome_var: string of outcome variable name
#  eg. "Y"
# model_formula: string of desired regression formula
#  eg. "Y~A+RISKGR1+SEX"
# missing_var: string of missing exposure variable name
# eg "D"
# df: dataset with no missing values 
#.   that contains variables in model_formula
#. e.g., sim.df
## outputs:
# list that contains parameters and full glm model
semilogit_func <- function(poi_var, outcome_var,model_formula,
                           missing_var="D",
                           oracle_formula = model_formula,oracle_flag=F,
                           learner= "HAL", df){
  
  # Remove Individuals with Missing Vaccination
  complete.df <- df %>% filter(get(missing_var)==0)
  
  # Run model using the oracle formula
  if(oracle_flag==T){
    model_formula <- oracle_formula
  }
  
  ## Run Partially Linear Logistic Regression Model
  # To get data in right format, create dataset with 
  # only numeric and binary variables with model matrix
  # including outcome
  model_formula <- gsub(paste0(outcome_var," *~"),
                                paste0("~",outcome_var,"+"),
                                model_formula)
  
  model.df <- as.data.frame(model.matrix(eval(parse(text = model_formula)), data=complete.df))
  #set.seed(seed_val)
  if(learner=="HAL"){
  mod <- spglm( ~1, model.df,
                W = setdiff(names(model.df), c("(Intercept)",outcome_var,poi_var)), 
                A = poi_var, Y = outcome_var,
                learning_method = learner,
                estimand = "OR",
                HAL_args_Y0W = list(smoothness_orders = 1, 
                                    max_degree = 2, 
                                    num_knots = c(10, 10)))
  } else{
    mod <- spglm( ~1, model.df,
                  W = setdiff(names(model.df), c("(Intercept)",outcome_var,poi_var)), 
                  A = poi_var, Y = outcome_var,
                  learning_method = learner,
                  estimand = "OR")
  }
  
  ## Calculating relevant statistics
  return(list(StatMethod = "Semiparametric", 
              Beta = coef(mod)$tmle_est,
              SE = coef(mod)$se,
              VE = (1-coef(mod)$psi_exp)*100,
              VELL = (1-coef(mod)$upper_exp)*100,
              VEUL = (1-coef(mod)$lower_exp)*100,
              PValue = coef(mod)$p_value,
              Reject = as.integer(coef(mod)$p_value<0.05),
              nsamp = nrow(df),
              ncomplete = nrow(complete.df),
              nVac = with(complete.df, sum(get(outcome_var)==1)),
              nCase = with(complete.df, sum(get(poi_var)==1)),
              nW = with(complete.df, sum(W==1)),
              nC = with(complete.df, sum(C==1)),
              nVacCase = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1)),
              nUnvacCase = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1)),
              nVacCaseCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1&RISKGR1=="At Risk")),
              nUnvacCaseCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1&RISKGR1=="At Risk")),
              nVacCaseNoCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==1&RISKGR1=="Not At Risk")),
              nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==1&RISKGR1=="Not At Risk")),
              nVacNoncaseCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==0&RISKGR1=="At Risk")),
              nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==0&RISKGR1=="At Risk")),
              nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var)==1&get(poi_var)==0&RISKGR1=="Not At Risk")),
              nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var)==0&get(poi_var)==0&RISKGR1=="Not At Risk")),
              missing_var = missing_var, 
              Formula = model_formula,
              Model = mod)) }



