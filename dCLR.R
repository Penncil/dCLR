# This is the code written by Jessie for dCLR


###############################
#### Step One: Compile the C code
###############################
# 1. open terminal
# 2. go the folder with the C code
# 3. R CMD SHLIB <xxx.c>
# (you will get .dll files in the current folder if you are using Windows)

# ++++++++++ possible errors: 
# “R” is not recognized as an internal or external command,
# solution: add path in the command, "E:/Apps/R/R-3.5.1/bin/R" R CMD SHLIB <file.c> (the R version and path migh not be the same)

###############################
#### Step Two: Load .so files and install library
###############################
dyn.load("mylik_odal.so")  
dyn.load("mylik_gradient.so")
dyn.load("cov_indep_odal.so")
library(meta)
library(logistf)
library(matlib)


###############################
#### Step Three: clean the real data
###############################
# K: number of hospitals
# n: a list number of patients in each hospital 
# y: outcome
#   1. each row represents a hospital
#   2. each column represents a patient
# x_all: covariates
#   1. each row represents a hospital
#   2. each column represents a patient
#   3. all the covariate are column-combined together
# length_par: number of covariates
# Please refer the sample data i shared: outcome.csv and variables.csv for reference

# please change the following number to the correct ones
K <- 6 # number of hospitals
n <- c(111, 2222, 3333, 4444, 555, 666) # sample size for six sites
length_par <- 4 # number of covairate included in the analysis


###############################
#### Step Forur: main functions to run
###############################

####################################################
##### pw.odal is the function to run the proposed method ##### 
####################################################
pw.odal <- function(K, n, y, x_all, length_par){
  
  
  ####################################################
  ################### meta method ######################
  ####################################################
  beta_meta_list = matrix(0,nrow = dim(y)[1],ncol = length_par+1)
  se_meta_list = matrix(0,nrow = dim(y)[1],ncol = length_par+1)
  # Jessie changed the code here: 04.07.2020
  for (i in c(1:dim(y)[1])){
    X_each = matrix(unlist(x_all[i,])[!is.na(unlist(x_all[i,]))], ncol = length_par)
    Y_each = unlist(y[i,])[!is.na(unlist(y[i,]))]
    fit_each= logistf(Y_each ~ X_each)
    beta_meta_list[i,] = fit_each$coefficients
    se_meta_list[i,] = sqrt(diag(fit_each$var))
  }
  
  
  #estimate from meta-analysis
  beta_meta_fix = c()
  beta_meta_random = c()
  beta_meta_fix_lower = c()
  beta_meta_fix_upper = c()
  beta_meta_random_lower = c()
  beta_meta_random_upper = c()
  for (i in 1:dim(beta_meta_list)[2]){
    ################### fixed effect meta method ######################
    tmp = metagen(beta_meta_list[,i], se_meta_list[,i],
                  comb.fixed = TRUE,comb.random = TRUE,sm="OR")
    beta_meta_fix[i] = tmp$TE.fixed
    beta_meta_fix_lower[i] = tmp$lower.fixed
    beta_meta_fix_upper[i] = tmp$upper.fixed
    ################### random effect meta method ######################
    beta_meta_random[i] = tmp$TE.random
    beta_meta_random_lower[i] = tmp$lower.random
    beta_meta_random_upper[i] = tmp$upper.random
  }
  
  ######################### Method 1 #######################
  ####################################################
  # ######## MLE (gold standard) ####### !!!!!!!!!!!!!!!!!!~!!!!!!!!!!!!!!!!!! here is hardcoding 
  ####################################################
  # Jessie changed the code here (04/07/2020)
  n_max = max(n)
  
  y_input = c(t(y))[!is.na(c(t(y)))]
  x1_tmp = x_all[,1:n_max]; x1_input =  c(t(x1_tmp))[!is.na(c(t(x1_tmp)))]
  x2_tmp = x_all[,((max(n)+1):(2*max(n)))]; x2_input =  c(t(x2_tmp))[!is.na(c(t(x2_tmp)))]
  x3_tmp = x_all[,((2*max(n)+1):(3*max(n)))]; x3_input =  c(t(x3_tmp))[!is.na(c(t(x3_tmp)))]
  x4_tmp = x_all[,((3*max(n)+1):(4*max(n)))]; x4_input =  c(t(x4_tmp))[!is.na(c(t(x4_tmp)))]
  dat_group =  unlist(c(mapply(rep, 1:K, n))) # site number

  dat = data.frame(y = y_input, 
                   x1 = x1_input, 
                   x2 = x2_input, 
                   x3 = x3_input,
                   x4 = x4_input,
                   group=dat_group)
  
  fit2 <- glmmboot(y ~ x1 + x2 + x3 + x4, cluster = group, data = dat)
  espar_esm_gold_1 =  c(fit2$coefficients)
  espar_esm_gold_var_1 = c(diag(fit2$variance))
  
  ####################################################
  # ######## stratified method (gold standard) #######
  ####################################################
  ######################### Method 2 #######################
  lik_gold=function(par){
    lik=rep(0,K)
    for(k in 1:K){
      temp<- .C("mylik_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par),result=double(1))
      lik[k]=temp[["result"]]
    }
    return(sum(lik))
  }
  
  # get the result for local site
  tryCatch(
    {
      espar_esm_gold_2=NA
      op_gold=optim(beta_meta_random[-1],
                    lik_gold,
                    control = list(fnscale=-1,maxit=1000),
                    method = "Nelder-Mead")
      espar_esm_gold_2=op_gold$par
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
    })
  espar_esm_gold_var_2 = var_func(espar_esm_gold_2)

  
  ######################### Method 3 #######################
  lik_gold_3=function(par){
    lik=rep(0,K)
    N_sum = 0
    logL = 0
    for(k in 1:K){
      temp<- .C("mylik_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n = n[k]*(n[k]-1)/2
      N_sum = N_sum + each_n
      lik[k]=temp[["result"]]/each_n
      logL <- logL + each_n*lik[k]
    }
    return(logL/N_sum)
  }
  
  
  # get the result for local site
  tryCatch(
    {
      espar_esm_gold_3=NA
      op_gold=optim(beta_meta_random[-1],
                    lik_gold_3,
                    control = list(fnscale=-1,maxit=1000),
                    method = "Nelder-Mead")
      espar_esm_gold_3=op_gold$par
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
    })
  espar_esm_gold_var_3 = var_func(espar_esm_gold_3)
  
  ####################################################
  ################## ODAL + pairwise ####################
  ####################################################
  ############ step 2: #########################
  ########### Initiation #######################
  # estimate the local within each site & return the variance to get the weighting initial value
  lik_local_list = matrix(NA, nrow = K, ncol = length_par)
  local_var_list = matrix(NA, nrow = K, ncol = length_par^2)
  for (local_num in 1:K){
    # local likelihood
    lik_local=function(par){
      temp<- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),
                as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n = n[local_num]*(n[local_num]-1)/2
      lik=temp[["result"]]/each_n
      return(lik)
    }
    # get the result for local site
    tryCatch(
      {
        espar_esm_local=NA
        op_local=optim(espar_esm_gold_1,
                       lik_local, 
                       control = list(fnscale=-1,maxit=1000),
                       method = "Nelder-Mead")
        espar_esm_local=op_local$par
        lik_local_list[local_num,] = espar_esm_local
        local_var_list[local_num,] = var_func(espar_esm_local)
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
    
  }
  
  #### weighting of the broadcast value
  # weighted average ### first gradient
  est_tmp = rep(0, length_par)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = local_var_list[index,]
    tmp2 = lik_local_list[index,]
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
      tryCatch(
        {
          est_tmp = est_tmp + t(Ginv(matrix(tmp1, length_par, length_par))%*% tmp2)
          est_tmp_2 = est_tmp_2 + as.vector(Ginv(matrix(tmp1, length_par, length_par)))
        },error=function(e){
          cat("ERROR :",conditionMessage(e), "\n")
        })
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  esm_init_bc = c(t(Ginv(matrix(est_tmp_2,length_par,length_par)) %*% t(est_tmp)))
  
  ############ step 3: #########################
  ## with the beta-bar obtain the gradients from all the sites with the initial values
  # first order 
  ourlik_all=function(par){
    grad=matrix(0,nrow = K,ncol=length_par)
    for(k in 1:K){
      temp<- .C("mylik_gradient",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par),result=double(length_par))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  grads = ourlik_all(esm_init_bc)
  mean_grads = apply(grads, 2, mean)
  
  
  # second order
  ourlik_all_second=function(par)
  {
    grad=matrix(0,nrow = K,ncol=length_par^2)
    for(k in 1:K){
      temp<- .C("mylik_gradient_second",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par), result=double(length_par^2))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  second_grads = ourlik_all_second(esm_init_bc)
  mean_2_grads = matrix(c(apply(second_grads, 2, mean)), length_par, length_par)
  
  
  
  # every site its own data and the 
  est_matrix_second = matrix(NA, nrow = K, ncol = length_par)
  first_grd_var = second_grd_var = matrix(NA, nrow = K, ncol = length_par^2)
  for (local_num in 1:K){
    print(local_num)
    # local likelihood
    lik_local=function(par){
      temp<- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n = n[local_num]*(n[local_num]-1)/2
      lik=temp[["result"]]/each_n
      return(lik)
    }
    
    ############ step 4: #########################
    ## construct surrogate likelihood function
    # l̃1(β)=l1(β)+{∇l(β¯)−∇l1(β¯)}(β−β¯)+(β−β¯)T{∇2l(β¯)−∇2l1(β¯)}(β−β¯)/2,
    surr_lik_second = function(par){
      lik_local(par) +
        (mean_grads - grads[1,])%*%(par - esm_init_bc) +
        (0.5*(t(par - esm_init_bc)%*%(mean_2_grads - matrix(second_grads[1,],length_par,length_par)) %*%(par - esm_init_bc)))
    }

    
    tryCatch(
      {
        esm_final_tmp_second = NA
        op_final_tmp_second=optim(esm_init_bc,
                                  surr_lik_second, 
                                  control = list(fnscale=-1,maxit=1000),
                                  method = "Nelder-Mead")
        esm_final_tmp_second = op_final_tmp_second$par
        
        # assign the final answer
        est_matrix_second[local_num, ] = esm_final_tmp_second
        second_grd_var[local_num,] = c(var_func(esm_final_tmp_second))
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
    
    
    
  }
  
  #### weighting of the final value
  # weighted average ### first gradient
  #### weighting of the broadcast value
  # weighted average ### second gradient
  est_tmp = rep(0, length_par)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = second_grd_var[index,]
    tmp2 = est_matrix_second[index,]
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
      est_tmp = est_tmp + t(solve(matrix(tmp1, length_par, length_par))%*% est_matrix_second[index,])
      est_tmp_2 = est_tmp_2 + as.vector(solve(matrix(tmp1, length_par, length_par)))
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  esm_final_second = c(t(solve(matrix(est_tmp_2,length_par,length_par)) %*% t(est_tmp)))
  esm_final_second_var = var_func(esm_final_second)
  
  
  
  ######## return final answer #########
  return(list(beta_meta_fix = beta_meta_fix,
              beta_meta_random = beta_meta_random,
              beta_meta_fix_lower = beta_meta_fix_lower,
              beta_meta_fix_upper = beta_meta_fix_upper,
              beta_meta_random_lower = beta_meta_random_lower,
              beta_meta_random_upper = beta_meta_random_upper,
              MLE_gold_standard = espar_esm_gold_1,
              MLE_espar_esm_gold_var = espar_esm_gold_var_1,
              stratified_gold_standard_old = espar_esm_gold_2,
              stratified_espar_esm_gold_var_old = espar_esm_gold_var_2,
              stratified_gold_standard_new = espar_esm_gold_3,
              stratified_espar_esm_gold_var_new = espar_esm_gold_var_3,
              esm_init_bc = esm_init_bc,
              second_gradient_matrix = second_grads,
              est_matrix_second = est_matrix_second,
              second_gradient_ODAL_pw = esm_final_second,
              beta_meta_list  = beta_meta_list,
              se_meta_list = se_meta_list,
              second_grd_var = second_grd_var,
              esm_final_second_var = esm_final_second_var))
  
}

###################################
## variance function (sandwich) ##### 
var_func <- function(par){
  #######
  # A
  #######
  dev_score_tmp=array(dim = c(K,length_par^2))
  for(k in 1:K){
    temp1<- .C("cal_dev_score_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
               as.double(par),as.double(length_par),result=double(length_par^2))
    dev_score_tmp[k,] = temp1[["result"]]
  }
  
  # diagnal of the A matrix
  dev_score_not_full=array(dim = c(K,length_par^2))
  for (i in 1:length_par){
    dev_score_not_full[,1+(length_par+1)*(i-1)] = dev_score_tmp[,i]
  }
  
  # funtion to get a symmatric matrix
  makeSymm <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  
  # rearrange the matrix
  dev_score=array(dim = c(K,length_par^2))
  for (row in 1:K){
    tmp_matrix <- matrix(dev_score_not_full[row,], nrow = length_par, byrow = TRUE)
    tmp_matrix[upper.tri(tmp_matrix, diag=FALSE)] <- dev_score_tmp[row, ((length_par+1):((((length_par^2)+length_par))/2))]
    tmp_matrix2 = makeSymm(tmp_matrix)
    dev_score[row,] = as.vector(tmp_matrix2)
  }
  
  # A in sandwich
  A.tmp = apply(dev_score, 2, sum)
  A = matrix(A.tmp, nrow = length_par, byrow = TRUE)
  
  #######
  # B
  #######
  score=array(dim = c(K,length_par))
  for(k in 1:K){
    temp2<- .C("cal_score_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
               as.double(par),as.double(length_par),result=double(length_par))
    score[k,] = temp2[["result"]]
  }
  B = t(score) %*% score 
  
  # get the result
  tmp = solve(A) %*% B %*% solve(A)
  esvar=  tmp
  
  return(esvar)
}

###################################
##### expit function ##### 
expit <- function(x){
  exp(x)/(1+exp(x))
}

####################################################
####### Run the function
####################################################
results = pw.odal(K, n, y, x_all, length_par)



