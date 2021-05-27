#' getMod function
#'
#' This function allows you to identify APA usage modalities.
#' The model used in the function is Gaussion Mixture model.
#'
#' @param PUIData a matrix which is the result of extrPairPA or exnon3UTRPA
#'
#' @keywords getMod
#' @export
#'
getMod <-
function(PUIData){

    tUTR.del <- as.matrix(PUIData)
    # print(dim(tUTR.del))
    #==========================================#
    dgmm2 <- apply(tUTR.del,1, function(x){
      # print(x)
      if(length(which(is.na(x)==TRUE))>0){
        dat <- x[-which(is.na(x))]
        if(length(dat) >=3){
          ocgmm <- Optimal_Clusters_GMM(as.data.frame(dat), 3,criterion = "BIC",
                                        dist_mode = "maha_dist",plot_data = F)
        }else{return(NULL)}



      }
      else{
        ocgmm <- Optimal_Clusters_GMM(as.data.frame(x), 3,criterion = "BIC",
                                      dist_mode = "maha_dist",plot_data = F)


      }
      # print(ocgmm)
      list(ordgmm2 = as.numeric(which.min((ocgmm))),dgmm2 = as.numeric(ocgmm))

    })
    ordgmm2 <- unlist(lapply(dgmm2,function(d){d[["ordgmm2"]]}))
    dgmm2 <- as.data.frame(lapply(dgmm2,function(d){d[["dgmm2"]]}))
    table(ordgmm2)

    comp <- c()
    for (l in 1:ncol(dgmm2)) {
      # l <- 1
      m <- dgmm2[,l]
      # v <- vary[,l]
      c <- m[-which.max(m)]
      # if(abs(v[1]-v[2]) <= 0.05 & abs(v[2]-v[3]) <= 0.05){
      #   comp <- c(comp,1)
      # }else
      if((sign(c)[1]*sign(c)[2]) > 0){
        if(round(c[which.min(abs(c))]/c[which.max(abs(c))],2) >= 0.95){
          if(which.max(m)==1){comp <- c(comp,2)}
          if(which.max(m)==2){comp <- c(comp,1)}
          if(which.max(m)==3){comp <- c(comp,1)}
        }else{comp <- c(comp,ordgmm2[l])}
      }else{comp <- c(comp,ordgmm2[l])}

    }
    table(comp)
    #==========================================#
    clust.lable2 <- lapply(c(1:nrow(tUTR.del)), function(x) {
      # x<-784
      if(length(which(is.na(tUTR.del[x,])==TRUE))>0){
        dat <- tUTR.del[x,][-which(is.na(tUTR.del[x,]))]
        gmm <- GMM(as.data.frame(dat) , gaussian_comps = comp[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(dat), gmm$centroids,
                          gmm$covariance_matrices, gmm$weights)


      }
      else{

        gmm <- GMM(as.data.frame(tUTR.del[x,]) , gaussian_comps = comp[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(tUTR.del[x,]), gmm$centroids,
                          gmm$covariance_matrices, gmm$weights)

      }

      return(list(gmm,pr))

    })
    ## correct component
    ## calculate cells' number in every component
    comCELLs <- lapply(c(1:length(clust.lable2)), function(i){
      x <- tUTR.del[i,]
      len <- length(clust.lable2[[i]][[1]][["centroids"]])
      id <- clust.lable2[[i]][[2]][["cluster_labels"]]
      if(len == 1){
        return(length(id))
      }else if(len == 2){
        return(c(length(which(id == 0)), length(which(id == 1))))
      }else{return(c(length(which(id == 0)), length(which(id == 1)), length(which(id == 2))))}

    })
    newdgmm2 <- unlist(lapply(c(1:length(comCELLs)), function(x){
      if(length(comCELLs[[x]]) == 2 & length(which(comCELLs[[x]] < 10)) > 0){
        comp[x] <- 1}else if(length(comCELLs[[x]]) == 3 & length(which(comCELLs[[x]] < 10))> 0){
          comp[x] <- 2
        }else{comp[x]}
    }))
    table(newdgmm2)
    #=================================================================#
    finaldgmm2 <- unlist(lapply(c(1:length(newdgmm2)), function(y){
      # y <- 71
      la <- clust.lable2[[y]][[2]]$cluster_labels
      if (length(table(la)) != newdgmm2[y]) {
        newdgmm2[y] <- length(table(la))
      }
      return(newdgmm2[y])
    }))
    table(finaldgmm2)
    clust.lable2 <- lapply(c(1:nrow(tUTR.del)), function(x) {
      # x <- 71
      if(length(which(is.na(tUTR.del[x,])==TRUE))>0){
        dat <- tUTR.del[x,][-which(is.na(tUTR.del[x,]))]
        gmm <- GMM(as.data.frame(dat) , gaussian_comps = finaldgmm2[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(dat), gmm$centroids,
                          gmm$covariance_matrices, gmm$weights)


      }
      else{
        gmm <- GMM(as.data.frame(tUTR.del[x,]) , gaussian_comps = finaldgmm2[x] , dist_mode = "maha_dist", seed_mode = "random_spread",
                   verbose = F )
        pr <- predict_GMM(as.data.frame(tUTR.del[x,]), gmm$centroids,
                          gmm$covariance_matrices, gmm$weights)

      }

      return(list(gmm,pr))

    })

    #=================================================================#
    # pattern
    # nq.vec <- comp
    nq.vec <- finaldgmm2
    pattern.type.all=lapply(c(1:length(nq.vec)), function(x){
      # x <- 32
      # 1 components
      if (nq.vec[x]==1){return ("Unimodal")}
      if (nq.vec[x]==2){
        # 2 components
        return("Bimodal")
      }else{
        # 3 components
        return ("Multimodal");
      }
    })
    pattern.type.all=unlist(pattern.type.all)
    table(pattern.type.all)
    return(list(modalities = pattern.type.all, results = clust.lable2))

}
