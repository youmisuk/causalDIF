# :::: Paper: Rethinking Item Fairness Using Single World Intervention Graphs ::::
# :::: Youmi Suk & Weicong Lyu :::::::::::::::::::::::::::::::::::::::::::::::::::

library(mirt)
library(difR)
library(xtable)
library(doParallel) # parallel cores
library(mgcv)
source("aux.R")

# :::: Simulations: Design 3 :::: ####

J = 60 
reps = 500 # the number of replications
tJ=J+J/2

design_value <- list(c(0, FALSE), c(0.5, FALSE), c(0, TRUE), c(0.5, TRUE)) # the first value = DIF effect; the second = the presence of impact

cores <- detectCores() - 1
cluster <- makeCluster(cores)

registerDoParallel(cluster)

result <- list()

foreach(bb = 1:reps, .packages = c('mgcv','nnet','foreach','doParallel','difR', 'caret', 'ltm', "mirt"), .combine = 'rbind') %dopar% {
  
   temp <- c()
  
for (d in 1:4) {
    
    print(paste("Desgin", d))  
    
    repeat{  
      itemlist <- DGP_D3(N = 2000, J = J, focal.frac = 0.5, DIFmag = design_value[[d]][1], itemImpact=design_value[[d]][2])
      
      df <- data.frame(theta=itemlist$theta, group=itemlist$group)
      
      tJ <- ncol(itemlist$Y)
      
      # Divide data into anchor items and pilot items
      item_data_theta <- data.frame(itemlist$Y[, 1:J]) # the first J items for ability estimation
      item_data_pilot <- data.frame(itemlist$Y[, -(tJ+1)])
      
      
      ability_lvl <- cut(itemlist$theta, quantile(itemlist$theta, seq(0, 1, 0.1)))
      
      # conduct MH test with the true ability categories
      MH_rlst <- tryCatch({difMH(Data=item_data_pilot, group=itemlist$group, focal.name=1,  match =ability_lvl)},
                          error = function(e){NULL} )
      
      if (design_value[[d]][2] == FALSE) {
        
        model_2pl <- mirt(data = item_data_theta, itemtype = "2PL", verbose = F, technical = list(NCYCLES = 500), covdata = data.frame(X=itemlist$X), formula = ~ X)
        theta.est <- fscores(model_2pl, method="EAP")[,1]
        
      } else {
        
        model_mg <- mirt(data = item_data_theta, itemtype = "2PL", verbose = F, technical = list(NCYCLES = 500), covdata = data.frame(X=itemlist$X, group=itemlist$group), formula = ~ X + group)
        theta.est <- fscores(model_mg, method="EAP")[,1]
        
      } 
 
      sum_score <- apply(item_data_theta, 1, sum)
      
      
      MH_rlst_sum <-  tryCatch({difMH(Data=item_data_pilot, group=itemlist$group, focal.name=1,  match =sum_score)},
                               error = function(e){NULL})
      
      
      if (!is.null(MH_rlst) & !is.null(MH_rlst_sum)) # 
        break  
      
      
    } 
    
    
    # true ability
    Logistic_rlst <- difLogistic(Data=item_data_pilot, group=itemlist$group, focal.name=1,  match =itemlist$theta)
    Std_rlst <- difStd(Data=item_data_pilot, group=itemlist$group, focal.name=1, match=ability_lvl)
    CDIF_rlst <- DIFmodel(itemdata=item_data_pilot, group = itemlist$group, theta =  itemlist$theta, confounder=data.frame(itemlist$X, itemlist$X2))  ##  
    CDIF_rlst2 <- DIFmodel(itemdata=item_data_pilot, group = itemlist$group, theta =  itemlist$theta, confounder=itemlist$X)  ##  
    
    # estimated ability for CDIF or sum scores for traditional DIF
    Logistic_rlst_sum <- difLogistic(Data=item_data_pilot, group=itemlist$group, focal.name=1,  match =sum_score)
    Std_rlst_sum <- difStd(Data=item_data_pilot, group=itemlist$group, focal.name=1, match=sum_score)
    CDIF_rlst_est <- DIFmodel(itemdata=item_data_pilot, group = itemlist$group, theta =  theta.est, confounder=data.frame(itemlist$X, itemlist$X2))  ##  
    CDIF_rlst_est2 <- DIFmodel(itemdata=item_data_pilot, group = itemlist$group, theta =  theta.est, confounder=itemlist$X)  ##  
    
    temp.est <- as.numeric(c(CDIF=CDIF_rlst[(J+1):tJ,1], CDIF_p=CDIF_rlst2[(J+1):tJ,1], Std=Std_rlst$PDIF[(J+1):tJ], MH=MH_rlst$p.value[(J+1):tJ], LR=Logistic_rlst$p.value[(J+1):tJ],
                        CDIF_est=CDIF_rlst_est[(J+1):tJ,1], CDIF_p_est=CDIF_rlst_est2[(J+1):tJ,1], Std_sum=Std_rlst_sum$PDIF[(J+1):tJ], MH_sum=MH_rlst_sum$p.value[(J+1):tJ], LR_sum=Logistic_rlst_sum$p.value[(J+1):tJ]))
    
    temp <- c(temp, temp.est)
    
  }
  
  result[as.character(bb)] <- temp
  
} -> simu_results

stopCluster(cluster)

# :::: Summarize results :::: ####

# Scenario 1 with true ability
rlst_D1 <- cbind.data.frame(RegAdj=apply(simu_results[,1:30] < 0.05, 2, mean),
                             RegAdj_p=apply(simu_results[,31:60] < 0.05, 2, mean),
                             Std=apply(abs(simu_results[,61:90]) > 0.05, 2, mean),
                             MH=apply(simu_results[,91:120] < 0.05, 2, mean),
                             Logistic=apply(simu_results[,121:150] < 0.05, 2, mean))

# Scenario 1 with estimated ability
rlst_D1_sum <- cbind.data.frame(RegAdj=apply(simu_results[,151:180] < 0.05, 2, mean),
                                RegAdj_p=apply(simu_results[,181:210] < 0.05, 2, mean),
                                Std=apply(abs(simu_results[,211:240]) > 0.05, 2, mean),
                                MH=apply(simu_results[,241:270] < 0.05, 2, mean),
                                Logistic=apply(simu_results[,271:300] < 0.05, 2, mean))

# Scenario 2 with true ability
rlst_D2 <- cbind.data.frame(RegAdj=apply(simu_results[,301:330] < 0.05, 2, mean),
                             RegAdj_p=apply(simu_results[,331:360] < 0.05, 2, mean),
                             Std=apply(abs(simu_results[,361:390]) > 0.05, 2, mean),
                             MH=apply(simu_results[,391:420] < 0.05, 2, mean),
                             Logistic=apply(simu_results[,421:450] < 0.05, 2, mean))

# Scenario 2 with estimated ability
rlst_D2_sum <- cbind.data.frame(RegAdj=apply(simu_results[,451:480] < 0.05, 2, mean),
                                RegAdj_p=apply(simu_results[,481:510] < 0.05, 2, mean),
                                Std=apply(abs(simu_results[,511:540]) > 0.05, 2, mean),
                                MH=apply(simu_results[,541:570] < 0.05, 2, mean),
                                Logistic=apply(simu_results[,571:600] < 0.05, 2, mean))

# Scenario 3 with true ability
rlst_D3 <- cbind.data.frame(RegAdj=apply(simu_results[,601:630] < 0.05, 2, mean),
                             RegAdj_p=apply(simu_results[,631:660] < 0.05, 2, mean),
                             Std=apply(abs(simu_results[,661:690]) > 0.05, 2, mean),
                             MH=apply(simu_results[,691:720] < 0.05, 2, mean),
                             Logistic=apply(simu_results[,721:750] < 0.05, 2, mean))

# Scenario 3 with estimated ability
rlst_D3_sum <- cbind.data.frame(RegAdj=apply(simu_results[,751:780] < 0.05, 2, mean),
                                RegAdj_p=apply(simu_results[,781:810] < 0.05, 2, mean),
                                Std=apply(abs(simu_results[,811:840]) > 0.05, 2, mean),
                                MH=apply(simu_results[,841:870] < 0.05, 2, mean),
                                Logistic=apply(simu_results[,871:900] < 0.05, 2, mean))

# Scenario 4 with true ability
rlst_D4 <- cbind.data.frame(RegAdj=apply(simu_results[,901:930] < 0.05, 2, mean),
                             RegAdj_p=apply(simu_results[,931:960] < 0.05, 2, mean),
                             Std=apply(abs(simu_results[,961:990]) > 0.05, 2, mean),
                             MH=apply(simu_results[,991:1020] < 0.05, 2, mean),
                             Logistic=apply(simu_results[,1021:1050] < 0.05, 2, mean))

# Scenario 4 with estimated ability
rlst_D4_sum <- cbind.data.frame(RegAdj=apply(simu_results[,1051:1080] < 0.05, 2, mean),
                                RegAdj_p=apply(simu_results[,1081:1110] < 0.05, 2, mean),
                                Std=apply(abs(simu_results[,1111:1140]) > 0.05, 2, mean),
                                MH=apply(simu_results[,1141:1170] < 0.05, 2, mean),
                                Logistic=apply(simu_results[,1171:1200] < 0.05, 2, mean))


rlst_true <- cbind.data.frame(Scenario1=apply(rlst_D1, 2, mean),
                              Scenario2=apply(rlst_D2, 2, mean),
                              Scenario3=apply(rlst_D3, 2, mean),
                              Scenario4=apply(rlst_D4, 2, mean))

round(rlst_true, 3)

rlst_true <- t(rlst_true)
rownames(rlst_true) <- c("Neither CDIF nor impact", "CDIF only", "Impact only", "Both CDIF and impact")
xtable(rlst_true, digits=3)


rlst_sum <- cbind.data.frame(Scenario1=apply(rlst_D1_sum, 2, mean),
                              Scenario2=apply(rlst_D2_sum, 2, mean),
                              Scenario3=apply(rlst_D3_sum, 2, mean),
                              Scenario4=apply(rlst_D4_sum, 2, mean))

round(rlst_sum, 3) 

rlst_sum <- t(rlst_sum)
rownames(rlst_sum) <-  c("Neither CDIF nor impact", "CDIF only", "Impact only", "Both CDIF and impact")
xtable(rlst_sum, digits=3)
