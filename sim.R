################################################################################
## The simulation code for manuscript "DIncorporating Sequential and Concurrent 
## Treatment Regimens in a Small n Sequential, Multiple Assignment, Randomized
## Trial (snSMART) with Application to the MISTIC Trial" 
################################################################################
# The code was written by Yuwei Cheng from the Department of Biostatistics at 
# the University of California San Diego. The code is modified and based on the  
# script written by Yan-Cheng Chao for his paper "Dynamic Treatment Regimens in 
# small n, Sequential, Multiple Assignment, Randomized Trials: an application in 
# focal segmental glomerulosclerosis". 
# Last updated on Jan 1st 2023. 
################################################################################


######## Setup ########
library(tidyverse)
library(rjags)

n.sim = 1000
NUM_ARMS = 3
NUM_COMBO = 13
MCMC_SAMPLE = 6000
BURN.IN = 1000
n_MCMC_chain = 1
nA = nB = nC = 60

# prior for first stage response rate 
# pi ~ uniform(0,1), beta0 ~ uniform(0,2), beta1 ~ uniform(0,2)
pi_prior.a = c(0,0,0)
pi_prior.b = c(1,1,1)  
beta0_prior.a = 0
beta0_prior.b = 2      
beta1_prior.r = 0         
beta1_prior.mu = 2 

# Since we assume only concurrent treatment involving A are affected (AB,AC,BA,BC), 
# Thus the priors for all delta_1,delta_seq, delta_con["BC"] and delta_con["CB"] 
# follow a uniform (0.9, 1.1) which allow them not to be impacted by the type parameter
# The deltas for the rest delta_con follows 
delta_1_prior.a=c(0.9,0,0.9,0,0.9,0,0.9,0.9,0.9,0,0.9,0.9,0.9)
delta_1_prior.b=c(1.1,3.5,1.1,3.5,1.1,3.5,1.1,1.1,1.1,3.5,1.1,1.1,1.1)

# set 1
pi_1A = 0.4
pi_1B = 0.2        
pi_1C = 0.2

# set 2 
# pi_1A = 0.45
# pi_1B = 0.3        
# pi_1C = 0.20

# set 3
# pi_1A = 0.45
# pi_1B = 0.45        
# pi_1C = 0.20

# set 4 : MISTIC
# pi_1A = 0.75
# pi_1B = 0.5        
# pi_1C = 0.5

# set 5 
# pi_1A = 0.6
# pi_1B = 0.4      
# pi_1C = 0.3

betas = list(
  "AA" = 1,
  "BB" = 1,
  "CC" = 1,
  "AB" = 0.8,
  "AC" = 0.8,
  "BA" = 0.6,
  "BC" = 0.6,
  "CA" = 0.4,
  "CB" = 0.4
)

# deltas for sequential DTRs are 1 for all cases
deltas_seq = list(
  "AB" = 1,
  "AC" = 1,
  "BA" = 1,
  "BC" = 1,
  "CA" = 1,
  "CB" = 1
)

# delta for first stage responders (AA BB and CC)
delta_1 = 1

# deltas for concurrent DTRs are differed by cases
# Marginal scenario
deltas_con = list(
  "AB" = (pi_1A * betas$AA)/(pi_1B * betas$AB),
  "AC" = (pi_1A * betas$AA)/(pi_1C * betas$AC),
  "BA" = (pi_1A * betas$AA)/(pi_1A * betas$BA),
  "BC" = 1,
  "CA" = (pi_1A * betas$AA)/(pi_1A * betas$CA),
  "CB" = 1
)

# Plateau scenario
# plateau_rate = 0.8
# deltas_con = list(
#   "AB" = (pi_1A * betas$AA)/(pi_1B * betas$AB)*plateau_rate,
#   "AC" = (pi_1A * betas$AA)/(pi_1C * betas$AC)*plateau_rate,
#   "BA" = (pi_1A * betas$AA)/(pi_1A * betas$BA)*plateau_rate,
#   "BC" = 1,
#   "CA" = (pi_1A * betas$AA)/(pi_1A * betas$CA)*plateau_rate,
#   "CB" = 1
# )

# Boost scenario
# boost_rate = 1.2
# deltas_con = list(
#   "AB" = (pi_1A * betas$AA)/(pi_1B * betas$AB)*boost_rate,
#   "AC" = (pi_1A * betas$AA)/(pi_1C * betas$AC)*boost_rate,
#   "BA" = (pi_1A * betas$AA)/(pi_1A * betas$BA)*boost_rate,
#   "BC" = 1,
#   "CA" = (pi_1A * betas$AA)/(pi_1A * betas$CA)*boost_rate,
#   "CB" = 1
# )


######## Data Generating function ########

data_generation = function(pi_1A , pi_1B, pi_1C,
                          betas, deltas_con,deltas_seq, delta_1,
                          nA, nB, nC) {
  
  # Generate second stage response rates based on our simulation:
  # B, C marginal: A + B = A, A + C = A, vice-versa
  
  # Second stage is the same as first stage
  pi_AA <<- pi_1A * betas[["AA"]] * delta_1
  pi_BB <<- pi_1B * betas[["BB"]] * delta_1
  pi_CC <<- pi_1C * betas[["CC"]] * delta_1
  
  # Concurrent treatment 2nd stage response rates
  pi_AB_con <<- pi_1B * betas[["AB"]] * deltas_con[["AB"]]
  pi_AC_con <<- pi_1C * betas[["AC"]] * deltas_con[["AC"]]
  pi_BA_con <<- pi_1A * betas[["BA"]] * deltas_con[["BA"]]
  pi_BC_con <<- pi_1C * betas[["BC"]] * deltas_con[["BC"]]
  pi_CA_con <<- pi_1A * betas[["CA"]] * deltas_con[["CA"]]
  pi_CB_con <<- pi_1B * betas[["CB"]] * deltas_con[["CB"]]
  
  # Sequential treatment 2nd stage response rates
  pi_AB_seq <<- pi_1B * betas[["AB"]] * deltas_seq[["AB"]] 
  pi_AC_seq <<- pi_1C * betas[["AC"]] * deltas_seq[["AC"]] 
  pi_BA_seq <<- pi_1A * betas[["BA"]] * deltas_seq[["BA"]] 
  pi_BC_seq <<- pi_1C * betas[["BC"]] * deltas_seq[["BC"]]
  pi_CA_seq <<- pi_1A * betas[["CA"]] * deltas_seq[["CA"]]
  pi_CB_seq <<- pi_1B * betas[["CB"]] * deltas_seq[["CB"]]
  
  # Initialize data 
  sim_data = tibble(
    stage1 = c(rep(as.integer(1), nA), rep(as.integer(2), nB), rep(as.integer(3), nC)),
    stage1_pi = c(rep(pi_1A , nA), rep(pi_1B, nB), rep(pi_1C, nC)),
    stage1_resp = map_int(stage1_pi, function(pi) {
      rbinom(1, size = 1, prob = pi)
    }),
    
    stage2 = map2_int(stage1, stage1_resp, function(s1, resp) {
      
      # If response in first stage is success, stay on treatment
      if (resp == 1) { return(s1) }
      
      # If not, randomize equally to remaining arms
      remaining_arms = setdiff(1:3, s1)
      arm = sample(remaining_arms, size = 1, replace = FALSE)
      
      as.integer(arm)
    }),
    
    stage2_type = map2_chr(stage1, stage1_resp, function(s1, resp) {
      
      # If response in first stage is success, 
      if (resp == 1) { return("s1_response") }
      
      # Decide type of treatment with equal probability
      t = sample(c("seq", "con"), size = 1, replace = FALSE)
      
      t
    }),
    
    stage2_pi = pmap_dbl(list(stage1, stage2, stage2_type), function(s1, s2, type) {
      
      pi = case_when(
        s1 == 1 & s2 == 1 ~ pi_AA,
        s1 == 2 & s2 == 2 ~ pi_BB,
        s1 == 3 & s2 == 3 ~ pi_CC,
        
        s1 == 1 & s2 == 2 & type == "con" ~ pi_AB_con,
        s1 == 1 & s2 == 3 & type == "con" ~ pi_AC_con,
        s1 == 2 & s2 == 1 & type == "con" ~ pi_BA_con,
        s1 == 2 & s2 == 3 & type == "con" ~ pi_BC_con,
        s1 == 3 & s2 == 1 & type == "con" ~ pi_CA_con,
        s1 == 3 & s2 == 2 & type == "con" ~ pi_CB_con,
        
        s1 == 1 & s2 == 2 & type == "seq" ~ pi_AB_seq,
        s1 == 1 & s2 == 3 & type == "seq" ~ pi_AC_seq,
        s1 == 2 & s2 == 1 & type == "seq" ~ pi_BA_seq,
        s1 == 2 & s2 == 3 & type == "seq" ~ pi_BC_seq,
        s1 == 3 & s2 == 1 & type == "seq" ~ pi_CA_seq,
        s1 == 3 & s2 == 2 & type == "seq" ~ pi_CB_seq
        
       
      )
      
      pi
      
    }),
    stage2_resp = map_int(stage2_pi, function(pi) {
      rbinom(1, size = 1, prob = pi)
    }),
    
   id =  pmap_int(list(stage1,stage1_resp),function(t1,r1){
    id=if_else(r1==0,2*t1-1,2*t1)
    as.integer(id)
   }),

    idType = pmap_int(list(stage1, stage2, stage2_type), function(s1, s2, type) {
      
      idtype = case_when(
        s1 == 1 & s2 == 1 ~ 1,
        s1 == 2 & s2 == 2 ~ 1,
        s1 == 3 & s2 == 3 ~ 1,
        
        s1 == 1 & s2 == 2 & type == "con" ~ 2,
        s1 == 1 & s2 == 3 & type == "con" ~ 4,
        s1 == 2 & s2 == 1 & type == "con" ~ 6,
        s1 == 2 & s2 == 3 & type == "con" ~ 8,
        s1 == 3 & s2 == 1 & type == "con" ~ 10,
        s1 == 3 & s2 == 2 & type == "con" ~ 12,
        
        s1 == 1 & s2 == 2 & type == "seq" ~ 3,
        s1 == 1 & s2 == 3 & type == "seq" ~ 5,
        s1 == 2 & s2 == 1 & type == "seq" ~ 7,
        s1 == 2 & s2 == 3 & type == "seq" ~ 9,
        s1 == 3 & s2 == 1 & type == "seq" ~ 11,
        s1 == 3 & s2 == 2 & type == "seq" ~ 13 
        
        
      )
      
      as.integer(idtype)
    })
  )
  
  sim_data %>% select(stage1, stage1_resp, stage2, stage2_resp,stage2_type, id, idType)

}

######## Simulations ########

# initiate vars
error_count = 0
pi_hat_new = pi_DTR_est_new = pi_hat_OG = pi_DTR_est_OG = c()
# simulation 
for (i in 1:n.sim){
  
  # generate dataset
  set.seed(i+10000)
  mydata = data_generation(pi_1A , pi_1B , pi_1C, betas, deltas_con, deltas_seq, delta_1, nA, nB, nC)
  
  # model path

  jags_path = "C:/Users/ychen/OneDrive/Documents/winter rotation/MIS-C"
  jag.model.name.new = "BJSM_12dtr_con_seq.bug"   
  jag.model.name.OG = "BJSM_OG_model.bug"   

  # new model
  tryCatch({
    jag_new = jags.model(file.path(jags_path,jag.model.name.new),
                      data=list(n = nrow(mydata),
                                num_arms = NUM_ARMS,
                                num_combo=NUM_COMBO,
                                Y1 = mydata$stage1_resp,
                                Y2 = mydata$stage2_resp, 
                                treatment_stageI = mydata$stage1,
                                treatment_stageII = mydata$stage2,
                                response_stageI_id = mydata$id,
                                id_type=mydata$idType,
                                #prior
                                pi_prior.a = pi_prior.a,
                                pi_prior.b = pi_prior.b, # pi
                                beta0_prior.a = beta0_prior.a,
                                beta0_prior.b = beta0_prior.b,
                                beta1_prior.a = beta1_prior.r,   
                                beta1_prior.c = beta1_prior.mu, # linkage
                                delta_1_prior.a=delta_1_prior.a,
                                delta_1_prior.b=delta_1_prior.b
                                 # concurrent
                      ),n.chains=n_MCMC_chain,n.adapt = BURN.IN)   
       posterior_sample_new = coda.samples(jag_new,c('pi','beta','delta'),MCMC_SAMPLE)
  },
  warning = function(war){
    warning_count = warning_count + 1
    err_war_message = rbind(paste("The warning of new model", warning_count, " is: ", war))
  },
  error = function(err){
    error_count = error_count + 1
    err_war_message = rbind(paste("The error of new model", error_count, " is: ", err))
    error_ind = 1
  },
  finally = {
    print(i)     # show the number of iterations run 
  }
  )
  
  # posterior 
  posterior_new = posterior_sample_new[[1]]
  
  # estimates pi
  pi_hat_new = rbind(pi_hat_new,apply(posterior_new[,20:22],2,mean))
  
  # calculate estimated DTRs: pi_jjj'= pi_j(pi_j*beta_1j*delta_1)+(1-pi_j)(pi_j'*beta_0j*delta_k)
  pi_AB_tt_con = posterior_new[,20]^2*posterior_new[,2]*posterior_new[,7]+(1-posterior_new[,20])*posterior_new[,21]*posterior_new[,1]*posterior_new[,8] 
  pi_AB_tt_seq = posterior_new[,20]^2*posterior_new[,2]*posterior_new[,7]+(1-posterior_new[,20])*posterior_new[,21]*posterior_new[,1]*posterior_new[,9] 

  pi_AC_tt_con = posterior_new[,20]^2*posterior_new[,2]*posterior_new[,7]+(1-posterior_new[,20])*posterior_new[,22]*posterior_new[,1]*posterior_new[,10]
  pi_AC_tt_seq = posterior_new[,20]^2*posterior_new[,2]*posterior_new[,7]+(1-posterior_new[,20])*posterior_new[,22]*posterior_new[,1]*posterior_new[,11]

  pi_BA_tt_con = posterior_new[,21]^2*posterior_new[,4]*posterior_new[,7]+(1-posterior_new[,21])*posterior_new[,20]*posterior_new[,3]*posterior_new[,12]
  pi_BA_tt_seq = posterior_new[,21]^2*posterior_new[,4]*posterior_new[,7]+(1-posterior_new[,21])*posterior_new[,20]*posterior_new[,3]*posterior_new[,13]
  
  pi_BC_tt_con = posterior_new[,21]^2*posterior_new[,4]*posterior_new[,7]+(1-posterior_new[,21])*posterior_new[,22]*posterior_new[,3]*posterior_new[,14]
  pi_BC_tt_seq = posterior_new[,21]^2*posterior_new[,4]*posterior_new[,7]+(1-posterior_new[,21])*posterior_new[,22]*posterior_new[,3]*posterior_new[,15]
  
  pi_CA_tt_con = posterior_new[,22]^2*posterior_new[,6]*posterior_new[,7]+(1-posterior_new[,22])*posterior_new[,20]*posterior_new[,5]*posterior_new[,16]
  pi_CA_tt_seq = posterior_new[,22]^2*posterior_new[,6]*posterior_new[,7]+(1-posterior_new[,22])*posterior_new[,20]*posterior_new[,5]*posterior_new[,17]
  
  pi_CB_tt_con = posterior_new[,22]^2*posterior_new[,6]*posterior_new[,7]+(1-posterior_new[,22])*posterior_new[,21]*posterior_new[,5]*posterior_new[,18]
  pi_CB_tt_seq = posterior_new[,22]^2*posterior_new[,6]*posterior_new[,7]+(1-posterior_new[,22])*posterior_new[,21]*posterior_new[,5]*posterior_new[,19]

  pi_DTR_est_new = rbind(pi_DTR_est_new,c(mean(pi_AB_tt_con),mean(pi_AB_tt_seq),mean(pi_AC_tt_con),mean(pi_AC_tt_seq),mean(pi_BA_tt_con),mean(pi_BA_tt_seq),mean(pi_BC_tt_con),mean(pi_BC_tt_seq),mean(pi_CA_tt_con),mean(pi_CA_tt_seq),mean(pi_CB_tt_con),mean(pi_CB_tt_seq)))
  
  # OG's model
  tryCatch({
    jag_OG = jags.model(file.path(jags_path,jag.model.name.OG),
                      data=list(n = nrow(mydata),
                                num_arms = NUM_ARMS,
                                Y1 = mydata$stage1_resp,
                                Y2 = mydata$stage2_resp,
                                treatment_stageI = mydata$stage1,
                                treatment_stageII = mydata$stage2,
                                response_stageI_id = mydata$id,
                                #prior
                                pi_prior.a = pi_prior.a,
                                pi_prior.b = pi_prior.b,
                                beta0_prior.a = beta0_prior.a,
                                beta0_prior.b = beta0_prior.b,
                                beta1_prior.a = beta1_prior.r,   
                                beta1_prior.c = beta1_prior.mu 
                              
                      ),
                      n.chains=n_MCMC_chain,n.adapt = BURN.IN)   
    posterior_sample_OG = coda.samples(jag_OG,
                                     c('pi','beta'),
                                     MCMC_SAMPLE)
  },
  warning = function(war){
    warning_count = warning_count + 1
    err_war_message = rbind(paste("The warning of OG's model ", warning_count, " is: ", war))
  },
  error = function(err){
    error_count = error_count + 1
    err_war_message = rbind(paste("The error of OG's model", error_count, " is: ", err))
    error_ind = 1
  },
  finally = {
    print(i)     # show the number of iterations run 
  }
  )
  
  # posterior
  posterior_OG = posterior_sample_OG[[1]]
  
  # estimates pi
  pi_hat_OG = rbind(pi_hat_OG,apply(posterior_OG[,7:9],2,mean))
  
  # calculate estimated DTRs: pi_jjj'= pi_j(pi_j*beta_1j)+(1-pi_j)(pi_j'*beta_0j)
  pi_AB_tt = posterior_OG[,7]^2*posterior_OG[,2]+(1-posterior_OG[,7])*posterior_OG[,8]*posterior_OG[,1]  
  pi_AC_tt = posterior_OG[,7]^2*posterior_OG[,2]+(1-posterior_OG[,7])*posterior_OG[,9]*posterior_OG[,1]
  pi_BA_tt = posterior_OG[,8]^2*posterior_OG[,4]+(1-posterior_OG[,8])*posterior_OG[,7]*posterior_OG[,3]
  pi_BC_tt = posterior_OG[,8]^2*posterior_OG[,4]+(1-posterior_OG[,8])*posterior_OG[,9]*posterior_OG[,3]
  pi_CA_tt = posterior_OG[,9]^2*posterior_OG[,6]+(1-posterior_OG[,9])*posterior_OG[,7]*posterior_OG[,5]
  pi_CB_tt = posterior_OG[,9]^2*posterior_OG[,6]+(1-posterior_OG[,9])*posterior_OG[,8]*posterior_OG[,5]
  pi_DTR_est_OG = rbind(pi_DTR_est_OG,c(mean(pi_AB_tt),mean(pi_AC_tt),mean(pi_BA_tt),mean(pi_BC_tt),mean(pi_CA_tt),mean(pi_CB_tt)))
}


######## Output Estimates for Both Models ########
# Treatment estimates from the modified model 
trt_output_new = data.frame(Treatment=c("A","B","C"),
                                true_pi_new = c(pi_1A,pi_1B,pi_1C),
                             pi_new = apply(pi_hat_new,2,mean),
                             sd_new = apply(pi_hat_new,2,sd)
                             )
trt_output_new = trt_output_new%>%
  mutate(
    bias_new = pi_new - true_pi_new,
    rMSE_new = sqrt(sd_new^2 + bias_new^2),
    abs_bias_new = abs(bias_new),
    mean_abs_bias = mean(abs_bias_new),
    mean_rMSE = mean(rMSE_new)
)
  
trt_output_new$bias_new = trt_output_new$pi_new - trt_output_new$true_pi_new
trt_output_new$rMSE_new = sqrt(trt_output_new$sd_new^2 + trt_output_new$bias_new^2)
trt_output_new$abs_bias_new=abs(trt_output_new$bias_new)
trt_output_new$mean_abs_bias=mean(trt_output_new$abs_bias_new)
trt_output_new$mean_rMSE=mean(trt_output_new$rMSE_new)

# expected DTR response rates
expected_new = c()  
# AB.con
expected_new[1] = pi_1A * pi_AA + (1 - pi_1A) * pi_AB_con
# AB.seq
expected_new[2] = pi_1A * pi_AA + (1 - pi_1A) * pi_AB_seq
# AC.con
expected_new[3] = pi_1A * pi_AA + (1 - pi_1A) * pi_AC_con
# AC.seq
expected_new[4] = pi_1A * pi_AA + (1 - pi_1A) * pi_AC_seq
# BA.con
expected_new[5] = pi_1B * pi_BB + (1 - pi_1B) * pi_BA_con
# BA.seq
expected_new[6] = pi_1B * pi_BB + (1 - pi_1B) * pi_BA_seq
# BC.con
expected_new[7] = pi_1B * pi_BB + (1 - pi_1B) * pi_BC_con
# BC.seq
expected_new[8] = pi_1B * pi_BB + (1 - pi_1B) * pi_BC_seq
# CA.con
expected_new[9] = pi_1C * pi_CC + (1 - pi_1C) * pi_CA_con
# CA.seq
expected_new[10] = pi_1C * pi_CC + (1 - pi_1C) * pi_CA_seq
# CB.con
expected_new[11] = pi_1C * pi_CC + (1 - pi_1C) * pi_CB_con
# CB.seq
expected_new[12] = pi_1C * pi_CC + (1 - pi_1C) * pi_CB_seq

# DTR estimates from the modified model 
dtr_output_new = data.frame(DTR=c("AB_con","AB_seq","AC_con","AC_seq","BA_con","BA_seq","BC_con","BC_seq","CA_con","CA_seq","CB_con","CB_seq"),
                                true_pi_dtr_new = expected_new,
                                 pi_DTR_new = apply(pi_DTR_est_new,2,function(x) mean(x,na.rm=T)),
                                 sd_DTR_new = apply(pi_DTR_est_new,2,sd)
                                  )
dtr_output_new = dtr_output_new%>%
  mutate(
    bias_new = pi_DTR_new - true_pi_dtr_new,
    rMSE_new = sqrt(sd_DTR_new^2 + bias_new^2),
    abs_bias_new = abs(bias_new),
    mean_abs_bias = mean(abs_bias_new),
    mean_rMSE = mean(rMSE_new)
  )

dtr_output_new$bias_new = dtr_output_new$pi_DTR_new - dtr_output_new$true_pi_dtr_new
dtr_output_new$rMSE_new = sqrt(dtr_output_new$sd_DTR_new^2 + dtr_output_new$bias_new^2)
dtr_output_new$abs_bias_new=abs(dtr_output_new$bias_new)
dtr_output_new$mean_abs_bias=mean(dtr_output_new$abs_bias_new)
dtr_output_new$mean_rMSE=mean(dtr_output_new$rMSE_new)

# write.csv(dtr_output_new,"simulations/0.4_0.2_0.2/N180 0.4_0.2_0.2/marginal_dtr_12_new_180.csv")
# write.csv(trt_output_new,"simulations/0.4_0.2_0.2/N180 0.4_0.2_0.2/marginal_trt_12_new_180.csv")

# Treatment estimates from OG's model 
trt_output_OG = data.frame(Treatment=c("A","B","C"), 
                                     true_pi_OG = c(pi_1A,pi_1B,pi_1C),
                             pi_OG = apply(pi_hat_OG,2,mean),
                             sd_OG = apply(pi_hat_OG,2,sd)
                             )

trt_output_OG = trt_output_OG%>%
  mutate(
    bias_OG = pi_OG - true_pi_OG,
    rMSE_OG = sqrt(sd_OG^2 + bias_OG^2),
    abs_bias_OG = abs(bias_OG),
    mean_abs_bias = mean(abs_bias_OG),
    mean_rMSE = mean(rMSE_OG)
)

expected_OG = c()  # expected DTR response rates
expected_OG[1] = pi_1A * pi_AA + (1 - pi_1A) * (pi_AB_seq+pi_AB_con)/2
expected_OG[2] = pi_1A * pi_AA + (1 - pi_1A) * (pi_AC_seq+pi_AC_con)/2
expected_OG[3] = pi_1B * pi_BB + (1 - pi_1B) * (pi_BA_seq+pi_BA_con)/2
expected_OG[4] = pi_1B * pi_BB + (1 - pi_1B) * (pi_BC_seq+pi_BC_con)/2
expected_OG[5] = pi_1C * pi_CC + (1 - pi_1C) * (pi_CA_seq+pi_CA_con)/2
expected_OG[6] = pi_1C * pi_CC + (1 - pi_1C) * (pi_CB_seq+pi_CB_con)/2

# DTR estimates from the modified model 
dtr_output_OG = data.frame(DTR=c("AB","AC","BA","BC","CA","CB"),
                                     true_pi_dtr_OG = expected_OG,
                                 pi_DTR_OG = apply(pi_DTR_est_OG,2,function(x) mean(x,na.rm=T)),
                                 sd_DTR_OG = apply(pi_DTR_est_OG,2,sd)
                                  )

dtr_output_OG = dtr_output_OG%>%
  mutate(
    bias_OG = pi_DTR_OG - true_pi_dtr_OG,
    rMSE_OG = sqrt(sd_DTR_OG^2 + bias_OG^2),
    abs_bias_OG = abs(bias_OG),
    mean_abs_bias = mean(abs_bias_OG),
    mean_rMSE = mean(rMSE_OG)
  )

# write.csv(dtr_output_OG,"simulations/0.4_0.2_0.2/N180 0.4_0.2_0.2/marginal_dtr_OG_180.csv")
# write.csv(trt_output_OG,"simulations/0.4_0.2_0.2/N180 0.4_0.2_0.2/marginal_trt_OG_180.csv")







































