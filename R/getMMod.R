#' getMMod function
#'
#' This function allows you to identify APA usage modalities.
#' The model used in the function is Gaussion Mixture model.
#'
#' @param APAData a PACdataset
#' @param type select "PAmax" or "PAmin"
#'
#' @keywords getMMod
#' @export
#'
getMMod <-
function(APAdata,type){
  tUTR.del <- as.matrix(log2(APAdata$fildata+1))
  dgmm2 <- apply(tUTR.del,1, function(x){
    # return(x)
    # x <- tUTR.del[8,]
    if(length(which(is.na(x)==TRUE))>0){
      dat <- x[-which(is.na(x))]
      if(length(which(dat == 0)) > 0){
        dat <- dat[-which(dat == 0)]
      }
      if(length(dat) >=3){
        ocgmm <- Optimal_Clusters_GMM(as.data.frame(dat), 3,criterion = "BIC",
                                      dist_mode = "maha_dist",plot_data = F)
      }else{return(0)}
    }
    else{
      if(length(which(x == 0)) > 0){
        x <- x[-which(x == 0)]
      }
      if(length(x) >= 3){
        ocgmm <- Optimal_Clusters_GMM(as.data.frame(x), 3,criterion = "BIC",
                                      dist_mode = "maha_dist",plot_data = F)
      }else{return(0)}

    }
    # return(as.numeric(which.min(ocgmm)))
    # return(as.numeric(ocgmm))
    list(ordgmm2 = as.numeric(which.min((ocgmm))),dgmm2 = as.numeric(ocgmm))

  })
  ordgmm2 <- unlist(lapply(dgmm2,function(d){
    if(type(d) == "list"){d[["ordgmm2"]]}
    else{return(0)}}))

  if(type == "PAmax"){
    dgmm2 <- as.data.frame(lapply(dgmm2,function(d){
      if(type(d) == "list"){d[["dgmm2"]]}}))}
  if(type == "PAmin"){
    dgmm2.tmp <- lapply(dgmm2,function(d){
      if(type(d) == "list"){return(d[["dgmm2"]])}
        else{return(0)}})
    # dgmm2.tmp <- dgmm2.tmp[which(dgmm2.tmp != "NULL")]
    dgmm2 <- as.data.frame(dgmm2.tmp)}
  table(ordgmm2)
  #==========================================#
  comp <- c()
  for (l in 1:ncol(dgmm2)) {
    # l <- 1
    m <- dgmm2[,l]
    # v <- vary[,l]
    if(length(which(m == 0)) < 1){
      c <- m[-which.max(m)]
      if((sign(c)[1]*sign(c)[2]) > 0){
        if(round(c[which.min(abs(c))]/c[which.max(abs(c))],2) >= 0.95){
          if(which.max(m)==1){comp <- c(comp,2)}
          if(which.max(m)==2){comp <- c(comp,1)}
          if(which.max(m)==3){comp <- c(comp,1)}
        }else{comp <- c(comp,ordgmm2[l])}
      }else{comp <- c(comp,ordgmm2[l])}
    }else{comp <- c(comp,0)}

  }
  table(comp)
  #==========================================#
  clust.lable2 <- lapply(c(1:nrow(tUTR.del)), function(x) {
    # x <- 8
    dat <- tUTR.del[x,]
    if(comp[x] != 0){
      if(length(which(is.na(dat)==TRUE))>0){
        mat <- dat[-which(is.na(dat))]
        if(length(which(mat == 0)) > 0){
          mat <- mat[-which(mat == 0)]
        }
        gmm <- GMM(as.data.frame(mat) , gaussian_comps = comp[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(mat), as.matrix(gmm$centroids),
                          gmm$covariance_matrices, gmm$weights)


      }
      else{
        if(length(which(dat == 0)) > 0){
          dat <- dat[-which(dat == 0)]
        }
        gmm <- GMM(as.data.frame(dat) , gaussian_comps = comp[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(dat), gmm$centroids,
                          gmm$covariance_matrices, gmm$weights)

      }
      return(list(gmm,pr))
      # return(gmm)
    }

  })
  ## calculate cells' number in every component
  comCELLs <- lapply(c(1:length(clust.lable2)), function(i){
    # i <- 1
    x <- tUTR.del[i,]
    if(!is.null(clust.lable2[[i]])){
      len <- length(clust.lable2[[i]][[1]][["centroids"]])
      id <- clust.lable2[[i]][[2]][["cluster_labels"]]
      if(len == 1){
        return(length(id))
      }else if(len == 2){
        return(c(length(which(id == 0)), length(which(id == 1))))
      }else{return(c(length(which(id == 0)), length(which(id == 1)), length(which(id == 2))))}
    }

  })
  newdgmm2 <- unlist(lapply(c(1:length(comCELLs)), function(x){
    if(!is.null(comCELLs[[x]])){
      if(sum(comCELLs[[x]]) <= 10){
        comp[x] <- 1
      }else if(length(comCELLs[[x]]) == 2 & length(which(comCELLs[[x]] < 10)) > 0){
        comp[x] <- 1}else if(length(comCELLs[[x]]) == 3 & length(which(comCELLs[[x]] < 10))> 0){
          comp[x] <- 2
        }else{comp[x]}
    }else{comp[x] <- 0}

  }))
  table(newdgmm2)
  # finaldgmm2 <-unlist(lapply(c(1:length(newdgmm2)), function(y){
  #   y <- 1
  #   la <- clust.lable2[[y]][[2]]$cluster_labels
  #   if (length(table(la)) != newdgmm2[y]) {
  #     newdgmm2[y] <- length(table(la))
  #   }
  #   return(newdgmm2[y])
  # }))
  # table(finaldgmm2)
  #==========================================#
  nq.vec <- newdgmm2
  pattern.type.all=lapply(c(1:length(nq.vec)), function(x){
    # x <- 32
    # mycen=as.data.table(clust.lable2[[x]][[1]][["centroids"]]) ;
    # 1 components
    if(nq.vec[x]==0){return(0)}
    if (nq.vec[x]==1){return ("Unimodal")}
    if (nq.vec[x]==2){
      # 2 components
      return("Bimodal");
    }else{
      # 3 components
      return ("Multimodal");
    }
  })
  pattern.type.all=unlist(pattern.type.all)
  table(pattern.type.all)
  return(list(modalities = pattern.type.all, results = clust.lable2))
}
