#' Least Root Mean Squared Deviation for Molecule Conformations
#' @description 
#' LRMSD or RMSD calculation between two coordinate sets (3 by n matrices).
#' @usage  
#' rmsd_calc(pred, truth, n=NULL, m=NULL, atype='all', optimal=FALSE)
#' @param pred matrix containing predicted coordinates of atoms in protein conformation.
#' @param truth matrix containing true coordinates of atoms in protein conformation. 
#' @param n The Resno to start at for the rmsd calculation.
#' @param m The Resno to end at for the rmsd calculation.
#' @param atype Can be Null, CA, CaCNO, specifies the types of items to consider.
#' Null means consider all atoms.
#' @param optimal TRUE to apply optimal rotations as described in 
#' \url{https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures}
#' Otherwise calculates RMSD without rotations
#' @return
#' Returns the calculated LRMSD value and the rotation matrix used to achieve optimal rotation.
#'
#'
#' @examples
#' predicted <- bio3d::read.pdb("TR928prediction.pdb")
#' truthful <- bio3d::read.pdb("TR928truth.pdb")
#' rmsd_calc(predicted, truthful,n=6, m=8, 'all', T)
#' rmsd_calc(predicted, truthful,n=6, m=8, 'CA', T)
#' rmsd_calc(predicted, truthful,n=6, m=8, "CaCNO", T)
#' rmsd_calc(predicted, truthful,n=6, m=8, "CaCNO", F)
#' rmsd_calc(predicted, truthful, n=8, m=NULL, "all", T)
#' rmsd_calc(predicted, truthful, n=NULL, m=NULL, 'all', F)
#' rmsd_calc(predicted, truthful)
#'
#' @export


rmsd_calc <- function(pred, truth, n=NULL, m=NULL, atype='all', optimal=FALSE){
  
  ####Combine into one df ----------------------------------------------------------####
  prednew <- pred$atom[c("elety", "resid", "resno", "x", "y", "z")]
  truthnew <- truth$atom[c("elety", "resid", "resno", "x", "y", "z")]
  total <- merge(prednew,truthnew,by=c("elety", "resid", "resno")) # , all.x=TRUE)
  total <- as.data.frame.list(total)
  
  #####-------------------Select a range-----------------------------####  
  #Specific Resno/ Resno Range
  if(!is.null(n)){
    if(!is.null(m)){
      df<-dplyr::filter(total, resno >=n & resno <=m)
    }
    else{
      df<-dplyr::filter(total, resno >=n)
    }
  }
  else{
    #do entire range
    df<-total
  }
  
  #####----------------------atype selection and calculation of rmsd----------------------####
  #CA only 
  if('CA' %in% atype){
    df<- dplyr::filter(df, elety=="CA")
  }
  
  #CA, C, N, O
  if('CaCNO' %in% atype){
    df<- dplyr::filter(df, elety=="CA" | elety=="C" | elety=="N" | elety=="O")
  }
  ####----------------------------Chose the df (optimal/not optimal)----------------------####
  #If the user types the word optimal, then calculate the optimal matrix, else do nothing
  if(optimal){
    result <- LRMSD(t(as.matrix(df[c("x.x", "y.x", "z.x")])), t(as.matrix(df[c("x.y", "y.y", "z.y")])))
    rmsd <- result$RMSD
    U <- t(result$rotation)
  } 
  else{
    rmsd <- RMSD(t(as.matrix(df[c("x.x", "y.x", "z.x")])), t(as.matrix(df[c("x.y", "y.y", "z.y")])))
    U<-NULL
  } 
  resno_out<-list("Rmsd"=rmsd, "Optimal matrix"= U)
  return(resno_out)
}

rmsd_calc_old <- function(pred, truth, n=NULL, m=NULL, atype='all', optimal=FALSE){
  
  ####Combine into one df ---------------------------------------------------####
  prednew <- pred$atom[c("elety", "resid", "resno", "x", "y", "z")]
  truthnew <- truth$atom[c("elety", "resid", "resno", "x", "y", "z")]
  total <- merge(prednew,truthnew,by=c("elety", "resid", "resno")) # , all.x=TRUE)
  total <- as.data.frame.list(total)

  
  #####-----------------------------Select a range---------------------------####  
  #Specific Resno/ Resno Range
  if(!is.null(n)){
    if(!is.null(m)){
      df<-dplyr::filter(total, resno >=n & resno <=m)
    }
    else{
      df<-dplyr::filter(total, resno >=n)
    }
  }
  else{
    #do entire range
    df<-total
  }
  

  #####----------------atype selection and calculation of rmsd---------------####
  #all 
  if('all' %in% atype){
    NULL
  }
  
  #CA only 
  if('CA' %in% atype){
    df<-dplyr::filter(df, elety=="CA")
  }else{
    NULL
  }
  
  #CA, C, N, O
  if('CaCNO' %in% atype){
    df<-dplyr::filter(df, elety=="CA" | elety=="C" | elety=="N" | elety=="O")
  } else{
    NULL
  }
  
  ####----------------------------Chose the df (optimal/not optimal)----------------------####
  #If the user types the word optimal, then calculate the optimal matrix, else do nothing
  if(optimal){
    
    #remove na rows 
    df2<-df[complete.cases(df), ]; df2
    
    #Calculate the means of the predicted x, y, z coordinates
    meanpredx<-mean(df2$x.x, na.rm=T); meanpredx
    meanpredy<-mean(df2$y.x, na.rm=T); meanpredy
    meanpredz<-mean(df2$z.x, na.rm=T); meanpredz
    
    #Calculate the means and the truth x, y, z coordinates
    meantruthx<-mean(df2$x.y, na.rm=T); meantruthx
    meantruthy<-mean(df2$y.y, na.rm=T); meantruthy
    meantruthz<-mean(df2$z.y, na.rm=T); meantruthz
    
    #Subtract the means from xi, yi, zi (predicted)
    df2<-mutate(df2, predprimex=(x.x-meanpredx))
    df2<-mutate(df2, predprimey=(y.x-meanpredy))
    df2<-mutate(df2, predprimez=(z.x-meanpredz))
    
    #Subtract the means from xi, yi, zi (truth)
    df2<-mutate(df2, truthprimex=(x.y-meantruthx))
    df2<-mutate(df2, truthprimey=(y.y-meantruthy))
    df2<-mutate(df2, truthprimez=(z.y-meantruthz))
    
    #Create an indicator variable to be used a the row names in the calculation of the matrix
    # df2$indicator<-paste(df2$elety.x, df2$resid.x, df2$resno.x); df2
    df2$indicator<-paste(df2$elety, df2$resid, df2$resno); df2
    
    
    opt_data_P<-df2[c("indicator", "predprimex", "predprimey", "predprimez")]
    opt_data_P1<-opt_data_P[,-1]
    rownames(opt_data_P1) <- opt_data_P[,1]
    head(opt_data_P1)
    
    opt_data_Q<-df2[c("indicator","truthprimex", "truthprimey", "truthprimez")]
    opt_data_Q1<-opt_data_Q[,-1]
    rownames(opt_data_Q1) <- opt_data_Q[,1]
    head(opt_data_Q1)
    
    
    # information_df<-df2[c("elety.x", "resid.x", "resno.x", "elety.y", "resid.y", "resno.y")]
    information_df<-df2[c("elety", "resid", "resno")]
    
    #Part 2: Covariance Matrix A 
    #Pmat is the predicted matrix; Qmat is the truth matrix 
    Pmat<-data.matrix(opt_data_P1)
    Qmat<-data.matrix(opt_data_Q1)
    head(Pmat)
    head(Qmat)
    
    #Calculate the cross product of the matrices 
    A<-crossprod(Pmat, Qmat); A
    
    #Part 3: Optimal Rotation matrix 
    #SVD 
    SVDA<-svd(A); SVDA
    
    #Next, decide whether we need to correct our rotation matrix to ensure a right-handed coordinate system
    dd<-det(SVDA$u %*% t(SVDA$v)); dd
    sidd<-sign(dd)
    
    #Finally, calculate our optimal rotation matrix, U
    D<-matrix(c(1,0,0,0,1,0,0,0,sidd), nrow=3, ncol=3, byrow=TRUE)
    U<-(SVDA$u %*% D %*% t(SVDA$v)); U
    
    
    #multiply original by optimal matrix 
    opt_ready_rmsd_P<-Pmat%*%U; opt_ready_rmsd_P
    
    #Transform back into a data frame 
    opt_ready_rmsd_PQ<-cbind(opt_ready_rmsd_P, Qmat); opt_ready_rmsd_PQ
    opt_read_df<-data.frame(opt_ready_rmsd_PQ); head(opt_read_df)
    colnames(opt_read_df) <- c("x.x","y.x", "z.x", "x.y", "y.y", "z.y"); head(opt_read_df)
    rownames(opt_read_df) <- c()
    opt_read_df$rn<-NULL
    optimal_df<-cbind(opt_read_df, information_df)
    newlist <- list("Optimal matrix" = U, "Optimized Data Frame" = optimal_df)
    
    df<-optimal_df
    
  } 
  else{
    U<-NULL
  } 
  
  #####--------Calculating the rmsd------######
  
  sumsq<-(mutate(df, sumsq=((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2)))
  rmsd2<-mean(sumsq$sumsq, na.rm=TRUE)
  rmsd1<-sqrt(rmsd2)
  
  ####----------------------Returning information---------------------------####
  resno_out<-list("Rmsd"=rmsd1, "Optimal matrix"= U)
  return(resno_out)
}


####----------------------------Testing-------------------------------####
# predicted <- bio3d::read.pdb("TR928prediction.pdb")
# truthful <- bio3d::read.pdb("TR928truth.pdb")
# test11<-rmsd_calc(predicted, truthful,n=6, m=8, 'all', T); test11 
# test22<-rmsd_calc(predicted, truthful,n=6, m=8, 'CA', T); test22   
# test33<-rmsd_calc(predicted, truthful,n=6, m=8, "CaCNO", T); test33 
# test44<-rmsd_calc(predicted, truthful,n=6, m=8, "CaCNO", F); test44 
# test55<-rmsd_calc(predicted, truthful, n=8, m=NULL, "all", T); test55  #not working 
# test66<-rmsd_calc(predicted, truthful, n=NULL, m=NULL, 'all', F); test66
# test77<-rmsd_calc(predicted, truthful); test77

