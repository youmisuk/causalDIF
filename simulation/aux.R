##### data generation for Design 1

DGP_D1 <- function(N=1000, focal.frac = 1/4, J=10, DIFmag=0, itemImpact=FALSE) {
  
  # N = number of test-taker
  # J = number of DIF items 
  group <- c(rep(1, N * focal.frac), rep(0, N * (1- focal.frac))) # focal dummy
  
  theta <- rnorm(N, 0, 1) # true theta for all the examinees
  
  if (itemImpact) {
    group_mean <- group*0.5
    theta = rnorm(N, group_mean)
  }
  
  tJ = J+J/2
  
  # item parameters
  set.seed(0)
  b <- round(runif(tJ), 3) 
  set.seed(0)
  a <- rep(1, tJ) 
  
  set.seed(NULL)
  
  probs0 <- probs1 <- probs <- counter_probs <-  matrix(NA, nrow = N, ncol = tJ)
  Y0 <- Y1 <- Y <- matrix(NA, nrow = N, ncol = tJ)
  
  
  for(j in 1:tJ){
    probs_l <- theta - b[j]  # Item generation under the control, reference group
    probs0[,j] <- 1 / (1 + exp(-probs_l))
    Y0[,j] <- rbinom(N, size = 1, prob =  probs0[,j])
  }
  
  for(j in 1:J){
    probs_l <- theta - b[j]  # non-DIF items under the treatment, focal group
    probs1[, j] <- 1 / (1 + exp(-probs_l))
    Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
  }
  
  if (DIFmag==0) {
    
    for(j in (J+1):tJ){ # there are non-DIF items under the treatment group
      probs_l <-   theta- b[j]  
      probs1[,j] <- 1 / (1 + exp(-probs_l))
      Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
    }
    
  } else {
    
    for(j in (J+1):tJ){ # DIF items under the treatment group
      probs_l <-   theta- b[j]-DIFmag  
      probs1[,j] <- 1 / (1 + exp(-probs_l))
      Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
    }
    
  }
  
  for(j in 1:tJ){
    Y[,j] <- ifelse(group == 1, Y1[,j], Y0[, j])
    probs[,j] <- ifelse(group == 1, probs1[,j], probs0[, j])  
    counter_probs[,j] <- ifelse(group == 1, probs0[,j], probs1[, j])  
  }
  
  return(list(Y=Y, Y1=Y1, Y0=Y0, probs=probs, counter_probs=counter_probs, 
              probs1=probs1, probs0=probs0, theta=theta, group=group))
}


##### data generation for Design 2

DGP_D2 <- function(N=1000, focal.frac = 1/4, J=10, DIFmag=0.5, itemImpact=FALSE) {
  
  # N = number of test-taker
  # J = number of DIF items 
  X <- runif(N, -1, 1)
  
  group <- rbinom(N, size = 1, prob = (1 / (1 + exp(-0.8*X)))) # focal dummy
  
  theta <- rnorm(N, 0, 1) # true theta for all the examinees
  if (itemImpact) {
    group_mean <- group*0.5
    theta = rnorm(N, group_mean)
  }
  
  tJ = J+J/2
  
  # item parameters
  set.seed(0)
  b <- round(runif(tJ), 3) 
  set.seed(0)
  a <- rep(1, tJ) 
  
  set.seed(NULL)
  
  probs0 <- probs1 <- probs <- counter_probs <-  matrix(NA, nrow = N, ncol = tJ)
  Y0 <- Y1 <- Y <- matrix(NA, nrow = N, ncol = tJ)
  
  
  for(j in 1:J){
    
    probs_l <- theta - b[j]   # Item generation under the control, reference group
    probs0[,j] <- 1 / (1 + exp(-probs_l))
    Y0[,j] <- rbinom(N, size = 1, prob =  probs0[,j])
    
    probs_l <- theta - b[j]   # non-DIF items under the treatment, focal group
    probs1[, j] <- 1 / (1 + exp(-probs_l))
    Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
  }
  
  for(j in (J+1):tJ){  # DIF items
    probs_l <- theta - b[j] + 2*X  # under the control, reference group
    probs0[,j] <- 1 / (1 + exp(-probs_l))
    Y0[,j] <- rbinom(N, size = 1, prob =  probs0[,j])
    
    probs_l <-   theta - b[j]- DIFmag  +  2*X # under the treatment, focal group
    probs1[,j] <- 1 / (1 + exp(-probs_l))
    Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
  }
  
  
  for(j in 1:tJ){
    Y[,j] <- ifelse(group == 1, Y1[,j], Y0[, j])
    probs[,j] <- ifelse(group == 1, probs1[,j], probs0[, j])  
    counter_probs[,j] <- ifelse(group == 1, probs0[,j], probs1[, j])  
  }
  
  return(list(Y=Y, Y1=Y1, Y0=Y0, probs=probs, counter_probs=counter_probs, 
              probs1=probs1, probs0=probs0, theta=theta, group=group, X=X))
}


##### data generation for Design 3 - two confounders

DGP_D3 <- function(N=1000, focal.frac = 1/4, J=10, DIFmag=0.5, itemImpact=FALSE) {
  
  # N = number of test-taker
  # J = number of DIF items 
  
  X <- runif(N, -1, 1) # confounder between mediator and outcome
  
  X2 <- runif(N, -1, 1) # confounder between protected variable and outcome
 
  group <- rbinom(N, size = 1, prob = (1 / (1 + exp(-0.8*X2)))) # focal dummy
  
  theta <- rnorm(N, 1*X, 1) # true theta for all the examinees
  
  if (itemImpact) {
    group_mean <- group*0.5
    theta = rnorm(N, 1*X-group_mean)
  }
  
  tJ = J+J/2
  
  # item parameters
  set.seed(0)
  b <- round(runif(tJ), 3) 
  set.seed(0)
  a <- rep(1, tJ) 
  
  set.seed(NULL)
  
  probs0 <- probs1 <- probs <- counter_probs <-  matrix(NA, nrow = N, ncol = tJ)
  Y0 <- Y1 <- Y <- matrix(NA, nrow = N, ncol = tJ)
  
  
  for(j in 1:J){
    
    probs_l <- theta - b[j]   # Item generation under the control, reference group
    probs0[,j] <- 1 / (1 + exp(-probs_l))
    Y0[,j] <- rbinom(N, size = 1, prob =  probs0[,j])
    
    probs_l <- theta - b[j]   # non-DIF items under the treatment, focal group
    probs1[, j] <- 1 / (1 + exp(-probs_l))
    Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
  }
  
  for(j in (J+1):tJ){  # DIF items
    probs_l <- theta - b[j] + 1*X + 1*X2 # under the control, reference group
    probs0[,j] <- 1 / (1 + exp(-probs_l))
    Y0[,j] <- rbinom(N, size = 1, prob =  probs0[,j])
    
    probs_l <-   theta - b[j]- DIFmag  + 1*X + 1*X2 # under the treatment, focal group
    probs1[,j] <- 1 / (1 + exp(-probs_l))
    Y1[,j] <- rbinom(N, size = 1, prob = probs1[, j])
  }
  
  
  for(j in 1:tJ){
    Y[,j] <- ifelse(group == 1, Y1[,j], Y0[, j])
    probs[,j] <- ifelse(group == 1, probs1[,j], probs0[, j])  
    counter_probs[,j] <- ifelse(group == 1, probs0[,j], probs1[, j])  
  }
  
  return(list(Y=Y, Y1=Y1, Y0=Y0, probs=probs, counter_probs=counter_probs, 
              probs1=probs1, probs0=probs0, theta=theta, group=group, X=X, X2=X2))
}


##### DIF model

DIFmodel <- function(itemdata, group, theta, confounder=NULL, inter=FALSE) {
  
  J <- ncol(itemdata)
  
  if (is.null(confounder) & inter==FALSE) {
    
    rlst <- matrix(0, ncol=5, nrow=J)
    colnames(rlst) <- c("p_value", "(Intercept)", "theta", "group", "theta:group")
    
    for (i in 1:J) {
      
      temp.item <- itemdata[, i]
      
      M1 <- glm(temp.item ~ theta, family=binomial)
      M2 <- glm(temp.item ~ theta +  group + theta:group, family=binomial)
      
      temp_anova <- anova(M1, M2)
      
      temp.p.value <- pchisq(temp_anova$Deviance[2], df = temp_anova$Df[2], lower.tail=FALSE) 
      
      rlst[i, 1] <- temp.p.value
      
      if (temp.p.value > 0.05) {
        
        rlst[i, 2:3] <- M1$coefficients
        
      } else {
        
        rlst[i, 2:5] <- M2$coefficients
        
      }
      
    }
    
  } else if (!is.null(confounder)) {
    
    confounder <- as.matrix(confounder)
    nconf <- ncol(confounder)
    
    rlst <- matrix(0, ncol=5+nconf, nrow=J)
    colnames(rlst) <- c("p_value", "(Intercept)", "theta", paste0("confounder_", 1:nconf), "group", "theta:group")
    
    for (i in 1:J) {
      
      temp.item <- itemdata[, i]
      
      M1 <- glm(temp.item ~ theta + confounder, family=binomial)
      M2 <- glm(temp.item ~ theta + confounder +  group + theta:group, family=binomial)
      
      temp_anova <- anova(M1, M2)
      
      temp.p.value <- pchisq(temp_anova$Deviance[2], df = temp_anova$Df[2], lower.tail=FALSE) # DIF
      
      rlst[i, 1] <- temp.p.value
      
      if (temp.p.value > 0.05) {
        
        rlst[i, 2:(3+nconf)] <- M1$coefficients
        
      } else {
        
        rlst[i, 2:(5+nconf)] <- M2$coefficients
        
      }
      
    }
    
  }
  
  return(rlst)
  
}





