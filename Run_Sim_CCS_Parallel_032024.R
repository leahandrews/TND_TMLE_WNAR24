
# Set library path
.libPaths(c("/home/users/landrew2/R_lib"))

# install.packages(
#   "pbapply",# xgboost, "ranger", #"earth", # "origami", #"Rdpack",# c("R.utils","R.oo","delayed","rstackdeque", "checkmate","BBmisc"),
#   lib = "/home/users/landrew2/R_lib",
#   repos = "http://cran.us.r-project.org"
# )
#library(dplyr)
library(SimEngine)







run_on_cluster(
  first = {
    
    file_path <- "/home/users/landrew2/Dissertation/Project1/"
    #file_path <- "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/"
    #source("~/Desktop/UW Stuff/IS TND/Code/Simulation Params 031824.R")
    source(paste0(file_path,"Simulation_Params_031824.R"))
    source(paste0(file_path,"Simulation_Methods_032524.R"))
    #source("~/Desktop/UW Stuff/IS TND/Code/Simulation Methods.R") 
    nsims <- 1000 #500
    npop <- 50000
    ncores <- 20 #36 #10
    nseed <-  2024 #2025 #2024 #2023
    stop_error_flag <- TRUE
    
    sim <- new_sim()

    sim %<>% set_levels(
      Estimator = list(
        "Oracle" = list(func = "ordlogit_func", oracle_flag = T,
                                 learner="glm",
                                 missing_var = "D",
                                 model_formula = "true distribution"),
        "OrdLogit0" = list(func = "ordlogit_func",oracle_flag = F,
                          learner="glm",
                          missing_var = "D",
                          model_formula = "A~Y+FEMALE+RISKGR1+CALTIME" ),
        "SemiLogit0" = list(func = "semilogit_func", oracle_flag=F,
                           learner = "HAL",
                           missing_var="D",
                           model_formula = "A~Y+FEMALE+RISKGR1+CALTIME")),
        # "OrdLogitMCAR20" = list(func = "ordlogit_func",oracle_flag = F,
        #                    learner="glm",
        #                    missing_var = "DMCAR20",
        #                    model_formula = "A~Y+FEMALE+RISKGR1+CALTIME" ),
        # "SemiLogitMCAR20" = list(func = "semilogit_func", oracle_flag=F,
        #                        learner = "HAL",
        #                        missing_var="DMCAR20",
        #                        model_formula = "A~Y+FEMALE+RISKGR1+CALTIME"),
        # "OrdLogitMCAR50" = list(func = "ordlogit_func",oracle_flag = F,
        #                         learner="glm",
        #                         missing_var = "DMCAR50",
        #                         model_formula = "A~Y+FEMALE+RISKGR1+CALTIME" ),
        # "SemiLogitMCAR50" = list(func = "semilogit_func", oracle_flag=F,
        #                          learner = "HAL",
        #                          missing_var="DMCAR50",
        #                          model_formula = "A~Y+FEMALE+RISKGR1+CALTIME"),
        # "OrdLogitMAR20" = list(func = "ordlogit_func",oracle_flag = F,
        #                         learner="glm",
        #                         missing_var = "DMAR20",
        #                         model_formula = "A~Y+FEMALE+RISKGR1+CALTIME" ),
        # "SemiLogitMAR20" = list(func = "semilogit_func", oracle_flag=F,
        #                          learner = "HAL",
        #                          missing_var="DMAR20",
        #                          model_formula = "A~Y+FEMALE+RISKGR1+CALTIME"),
        # "OrdLogitMAR50" = list(func = "ordlogit_func",oracle_flag = F,
        #                         learner="glm",
        #                         missing_var = "DMAR50",
        #                         model_formula = "A~Y+FEMALE+RISKGR1+CALTIME" ),
        # "SemiLogitMAR50" = list(func = "semilogit_func", oracle_flag=F,
        #                          learner = "HAL",
        #                          missing_var="DMAR50",
        #                          model_formula = "A~Y+FEMALE+RISKGR1+CALTIME")),
        # "SemiLogitGLM" = list(func = "semilogit_func", oracle_flag=F,
        #                     learner = "glm",
        #                     missing_var="D",
        #                     model_formula = "A~Y+FEMALE+RISKGR1+CALTIME")),
      trueVE = c(0,50,90),
      nseed = c(nseed), # 2023
      distribution = list("Scenario 3" = list(a_betas = a_betas3,
                                              y_betas = y_betas3, 
                                              w_betas = w_betas,
                                              c_betas = c_betas,
                                              missing_rate_vec = c(0.2, 0.5),
                                              mar20_d_betas = comorb_mar20_d_betas,
                                              mar50_d_betas = comorb_mar50_d_betas,
                                              oracle_formula = "A~ Y+ FEMALE + RISKGR1+ FEMALE:RISKGR1+ CALTIME+CALTIMEs2+CALTIMEs3"),
                          "Scenario 2" = list(a_betas = a_betas2,
                                              y_betas = y_betas2, 
                                              w_betas = w_betas,
                                              c_betas = c_betas, 
                                              missing_rate_vec = c(0.2, 0.5),
                                              mar20_d_betas = comorb_mar20_d_betas,
                                              mar50_d_betas = comorb_mar50_d_betas,
                                              oracle_formula = "A~ Y+ FEMALE + RISKGR1+ FEMALE:RISKGR1+ CALTIME"),#"Y~ A+ FEMALE + RISKGR1+ FEMALE:RISKGR1+ CALTIME"),
                          "Scenario 1" = list(a_betas = a_betas1,
                                              y_betas = y_betas1, 
                                              w_betas = w_betas,
                                              c_betas = c_betas,
                                              missing_rate_vec = c(0.2, 0.5),
                                              mar20_d_betas = comorb_mar20_d_betas,
                                              mar50_d_betas = comorb_mar50_d_betas,
                                              oracle_formula = "A~ Y+ FEMALE + RISKGR1+ CALTIME")),
      nsamp = c(2000, 500, 1000, 2500)) #c(2000, 500), c(2000, 40000)
    
    sim %<>% set_config(
      num_sim = nsims,
      #parallel = FALSE,
      n_cores = ncores,
      batch_levels = c("trueVE","distribution","nsamp"), 
      return_batch_id = T,
      seed = nseed, #2023,
      stop_at_error= stop_error_flag,
      packages = c("causalglm","dplyr","readr")
    )
    
    # set_script is like a for loop
    sim %<>% set_script(function() {
      # Apply both methods to same dataset
      batch({
        data <- data_func(a_betas = L$distribution$a_betas, a_formula = a_formula,
                          a_y_beta = log(1-L$trueVE/100), 
                          not_a_y_betas = L$distribution$y_betas,
                          y_formula = y_formula,
                          a_w_beta = 0, 
                          not_a_w_betas = L$distribution$w_betas,
                          w_formula = w_formula,
                          a_c_beta = 0, 
                          not_a_c_betas = L$distribution$c_betas,
                          c_formula = c_formula,
                          missing_rate_vec = L$distribution$missing_rate_vec, #c(0.2, 0.5), 
                          mar20_d_betas = L$distribution$mar20_d_betas,
                          mar50_d_betas = L$distribution$mar50_d_betas,
                          mar_d_formula = comorb_mar_d_formula,
                          pop_size = npop, sample_size = L$nsamp, rct.df = simtest)
      })
      estimates <- use_method(L$Estimator$func, list(poi_var="Y",outcome_var="A",
                                                     model_formula = L$Estimator$model_formula, #"A~Y+FEMALE+RISKGR1+CALTIME", #"Y~A+RISKGR1+FEMALE+CALTIME", 
                                                     missing_var=L$Estimator$missing_var,
                                                     oracle_formula = L$distribution$oracle_formula,
                                                     oracle_flag = L$Estimator$oracle_flag,
                                                     learner = L$Estimator$learner,
                                                     df=data$Sample))
      return(list(
        "truelogOR" = log(1-L$trueVE/100),
        "logOR" = estimates$Beta,
        "logORSE" = estimates$SE,
        "logORVar" = (estimates$SE)^2, # new
        "VE" = estimates$VE,
        "VELL" = estimates$VELL,
        "VEUL" = estimates$VEUL,
        "PValue"= estimates$PValue,
        "Reject"= estimates$Reject,
        "PopMargA" = mean(data$Population$A),
        "PopMargY" = mean(data$Population$Y),
        "PopMargW" = mean(data$Population$W),
        "PopMargC" = mean(data$Population$C),
        "PopMargFem" = mean(data$Population$FEMALE),
        "PopMargCom" = prop.table(table(data$Population$RISKGR1))[2],
        "SampMargA" = estimates$nVac/estimates$ncomplete, #mean(data$Sample$A),
        "SampMargY" = estimates$nCase/estimates$ncomplete, #mean(data$Sample$Y),
        "SampMargW" = estimates$nW/estimates$ncomplete, #mean(data$Sample$W),
        "SampMargC" = estimates$nC/estimates$ncomplete, #mean(data$Sample$C),
        "SampMargnA" = estimates$nVac, #mean(data$Sample$A),
        "SampMargnY" = estimates$nCase, #mean(data$Sample$Y),
        "SampMargnW" = estimates$nW, #mean(data$Sample$W),
        "SampMargnC" = estimates$nC, #mean(data$Sample$C),
        "SampMargFem" = mean(data$Sample$FEMALE),
        "SampMargCom" = prop.table(table(data$Sample$RISKGR1))[2],
        "SampVacCase" = estimates$nVacCase, #with(data$Sample, sum(A==1&Y==1)),
        "SampUnvacCase" = estimates$nUnvacCase, #with(data$Sample, sum(A==0&Y==1)),
        "SampDMCAR20" = sum(data$Sample$DMCAR20==0),
        "SampDMCAR50" = sum(data$Sample$DMCAR50==0),
        "SampDMAR20" = sum(data$Sample$DMAR20==0),
        "SampDMAR50" = sum(data$Sample$DMAR50==0),
        "ModelSamp" = estimates$ncomplete,
        "MissingVar"= estimates$missing_var,
        "check" = data$Sample$SUBJID[1],
       # "samplesize"=nrow(data$Sample), # new
        "model_formula" = estimates$Formula,
       "SampVacCaseCom" = estimates$nVacCaseCom ,
       "SampUnvacCaseCom"=estimates$nUnvacCaseCom ,
       "SampVacCaseNoCom" =estimates$nVacCaseNoCom ,
       "SampUnvacCaseNoCom"=estimates$nUnvacCaseNoCom ,
       "SampVacNoncaseCom"= estimates$nVacNoncaseCom ,
       "SampUnvacNoncaseCom"=estimates$nUnvacNoncaseCom ,
       "SampVacNoncaseNoCom"=estimates$nVacNoncaseNoCom ,
       "SampUnvacNoncaseNoCom"= estimates$nUnvacNoncaseNoCom ,
        # Store larger, complex objects 
        ".complex" = list(
          # "Model" = estimates$Model,
          "Table" = with(data$Sample, table( Y,A,FEMALE, RISKGR1)))
        #"DataParams"= data$DataParams)
        # "Population" = data$Population,
        #  "Sample" = data$Sample)
      ))
    })
  },
  main = {
    sim %<>% run()
  },
  last = {
    
    result.tbl<- sim %>% SimEngine::summarize(
      list(stat="mean", x="logOR", name="logORMean"),
      list(stat="bias", truth="truelogOR", estimate="logOR", name="logORBias", mc_se=TRUE),
      list(stat="mean", x="truelogOR",name="truelogOR"),
      list(stat="mean", x="VE", name="VEMean"),
      list(stat="bias", truth="trueVE", estimate="VE", name="VEBias", mc_se=TRUE),
      list(stat="sd", x="logOR", name="logORMCSD"),
      list(stat="var", x="logOR", name="logORMCVar"),
      list(stat="mean", x="logORSE", name="logORSEMean"),
      list(stat="mean", x="logORVar", name="logORVarMean"),
      list(stat="coverage", lower="VELL", upper="VEUL", truth="trueVE", name="VECoverage", mc_se=TRUE),
      list(stat="mean", name="Power", x="Reject"),
      list(stat="mean", x="PopMargA", name="PopMargAMean"),
      list(stat="mean", x="PopMargY", name="PopMargYMean"),
      list(stat="mean", x="PopMargW", name="PopMargWMean"),
      list(stat="mean", x="PopMargC", name="PopMargCMean"),
      list(stat="mean", x="SampMargA", name="SampMargAMean"),
      list(stat="mean", x="SampMargY", name="SampMargYMean"),
      list(stat="mean", x="SampMargW", name="SampMargWMean"),
      list(stat="mean", x="SampMargC", name="SampMargCMean"))
    
    saveRDS(result.tbl, 
            file= paste0(file_path,"Output/","Sim_Summary_", npop/1000,"K_",nsims, "_sims_",Sys.Date(),".rds"))
    
    runtime<- sim %>% SimEngine::vars("total_runtime")
    print(paste0("Total run time is ",round(runtime/60), " minutes ", "(",round(runtime/60/60,2), " hours)"))
    
  },
  
  cluster_config = list(js="slurm")
)

saveRDS(sim$results, 
        file = paste0(file_path,"Output/","Sim_Raw_", npop/1000,"K_",nsims, "_sims_",Sys.Date(),".rds"))
saveRDS(sim$levels, 
        file = paste0(file_path,"Output/","Sim_Betas_", npop/1000,"K_",nsims, "_sims_",Sys.Date(),".rds"))

sim$warnings
sim$errors

runtime<- sim %>% SimEngine::vars("total_runtime")

print(paste0("This core's run time is ",round(runtime/60), " minutes ", "(",round(runtime/60/60,2), " hours)"))
