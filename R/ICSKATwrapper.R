ICSKATwrapper <- function(left_dmat, right_dmat, initValues, lt, rt, obs_ind, tpos_ind, gMat, PH=TRUE, nKnots=1, maxIter=3, returnNull = FALSE) {

	xMat <- left_dmat[, 1:(ncol(left_dmat) - nKnots - 2)]	
	counter <- 0
	pass <- FALSE
	while (counter < maxIter) {

		counter <- counter + 1
	
		# initial values 
		if (counter == 1) {
			init_beta <- initValues
		} else {
			init_beta <- runif(n=ncol(left_dmat), min = -1, max = 1)
		}

		# null fit
		if (PH) {
			nullFit <- ICSKAT_fit_null(init_beta=init_beta, lt=lt, rt=rt,
                              left_dmat=left_dmat, right_dmat=right_dmat,
                              obs_ind=obs_ind, tpos_ind=tpos_ind)
		}	else {
			nullFit <- ICSKAT_fit_null_PO(init_beta=init_beta,
                              lt=lt, rt=rt,
                              ZL=left_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)],
                              ZR=right_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)], xMat = xMat,
                              obs_ind=obs_ind, tpos_ind=tpos_ind)
		}

		# if null fit not good, go to next loop
		if ( nullFit$err == 1 ) {
			next
		}
	
		# test
		if (PH) {
			skatOutput <- ICskat(left_dmat=left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                        right_dmat=right_dmat, gMat=gMat, lt=lt, rt=rt,
                        null_beta=as.numeric(nullFit$beta_fit), Itt=nullFit$Itt)
		} else {
			skatOutput <- ICskatPO(ZL=left_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)],
                    ZR=right_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)], xMat = xMat,
                    lt=lt, rt=rt,
                    obs_ind=obs_ind, tpos_ind=tpos_ind, gMat=gMat, nullCoef=nullFit$nullCoef, Itt = nullFit$Itt)
		}

		# if worked, break
		if (skatOutput$err == 0 | skatOutput$err == 22) {
			pass <- TRUE	
			break
		}
	} # end while loop

	# did it work
	if (!pass) {
		if (nullFit$err == 1) {
      skatOutput <- list(p_SKAT=NA, p_burden=NA, complex=NA, err=1, errMsg="Failed null fit")
    } else {
      skatOutput <- list(p_SKAT=NA, p_burden=NA, complex=NA, err=1, errMsg="Failed testing")
    }	
	}

	# return
	if (returnNull) {
      return(list(skatOutput = skatOutput, nullFit = nullFit))
  } else {
      return(skatOutput)
  }

}

