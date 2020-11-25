#' Fit RO-logit model using Cox-PH
#'
#' @param design.mat The design matrix without intercept term.
#' @param y vector of numeric. The outcomes.
#' @param svar vector of numeric. The strata variable.
#' @param method string. Use Efron (\code{"efron"}) or Breslow
#'   (\code{"breslow"}) method for handling ties in the outcome. The default is
#'   \code{"efron"}. See \code{\link[survival]{coxph}} for details.
#'
#' @return Returns the model fitted using \code{coxph}.
#' @import survival
rologit.coxph<-function(design.mat,y,svar=NULL,method="efron"){
  method <- match.arg(tolower(method), choices = c("efron", "breslow"))
  n<-length(y)

  ## Fit a stratified model?
  if(is.null(svar)){#No. Do not fit a stratified model.
    # Create the simplest test data set
    ry<-rank(-1*y,ties.method = "average")
    test1 <- data.frame(y=ry,
                        status=rep(1,n),
                        design.mat)
    # R will automatically convert ":" to ".", but if we assign column names of
    # design.mat back we will have to deal with this ourselves, and this is
    # troublesome if there are more than 1 interaction terms.
    # colnames(test1)[-c(1:2)]<-colnames(design.mat)

    out<-coxph(as.formula(paste("Surv(y, status) ~", paste(colnames(test1)[-c(1:2)],collapse="+"))),data=test1, ties = method)

  }else{#Yes. Fit a stratified model.
    # Create the simplest test data set
    ry<-rep(NA,n)
    for(i.grp in unique(svar)){
      ii<-which(svar==i.grp)
      ry[ii]<-rank(-1*y[ii],ties.method = "average")
    }

    test1 <- data.frame(y=ry,
                        status=rep(1,n),
                        design.mat)
    # colnames(test1)[-c(1:2)]<-colnames(design.mat)

    out<-coxph(as.formula(paste("Surv(y, status) ~", paste(colnames(test1)[-c(1:2)],collapse="+"),"+strata(svar)")),data=test1, ties = method)
  }

  return(out)
}
#' The negative log likelihood function for obtaining heuristic residuals
#' @description Compute the negative log-likelihood for obtaining heuristic
#'   residuals.
#'
#' @param par vector of numeric. Contain the intercept, first entry, and
#'   coefficient of log(sigma), second entry.
#' @param y vector of numeric. Centered outcomes within each stratum.
#' @param x matrix. The design matrix including intercept term.
#' @param par1 vector of numeric. The estimated coefficients from RO-logit
#'   model.
loglikhresid<-function(par,y,x,par1){
  beta0 = par[1]
  sigma = exp(par[2])
  hat.y<-c(x%*%c(beta0,par1*sigma))
  res.y<-y-hat.y
  return(-sum(evd::dgumbel(res.y,loc=0,scale=sigma,log=TRUE)))
}
#' Obtain the heuristic residuals.
#' @description Obtain the heuristic residuals.
#'
#' @param y vector of numeric. The outcomes.
#' @param svar vector of numeric. The strata variable.
#' @param design.mat matrix. The design matrix including intercept term.
#' @param o matrix. Contains the estimated coefficients of the RO-logit and
#'   their SEs.
#' @param initial.res.par The initial values of the intercept and log(scale), to
#'   be passed to the \code{optim} function. The default values are set to
#'   \code{c(0, 0)}, yet users are recommended to try a few initial values to
#'   make sure global optimum is reached.
#' @param ... Other parameters to be passed to the \code{optim} function. See
#'   \code{\link[stats]{optim}} for details.
#' @return Returns a list containing the estimated intercept and
#'   \code{log(scale)}, the covariance matrix of these two parameters,
#'   convergence status from \code{optim}, and the heuristic residuals.
#' @importFrom stats as.formula model.matrix optim
hresid<-function(y,svar,design.mat,o,initial.res.par=c(0, 0), ...){
  design.mat.c<-design.mat
  if (is.null(svar)) {
    svar <- rep(1, length(y))
  }

  if(ncol(design.mat)==2){

    nexp<-ncol(design.mat)
    ##Get the mean exposure for each grouping
    mx<-tapply(design.mat[,nexp],svar,mean)
    mmx<-mx[match(svar,names(mx))]
    ##Center the outcome such that the matched effects are removed
    x.c<-design.mat[,nexp]-mmx

    design.mat.c[,nexp]<-x.c

  }
  if(ncol(design.mat)>2){

    for(nexp in 2:ncol(design.mat)){
      ##Get the mean exposure for each grouping
      mx<-tapply(design.mat[,nexp],svar,mean)
      mmx<-mx[match(svar,names(mx))]
      ##Center the outcome such that the matched effects are removed
      x.c<-design.mat[,nexp]-mmx

      design.mat.c[,nexp]<-x.c
    }
  }

  ##Get the mean outcome for each grouping
  my<-tapply(y,svar,mean)
  mmy<-my[match(svar,names(my))]
  ##Center the outcome such that the matched effects are removed
  y.c<-y-mmy

  sol<-optim(par = initial.res.par,fn = loglikhresid,hessian=T,y=y.c,x=design.mat.c,par1=o[,1], ...)
  resid<-y.c-c(design.mat.c%*%c(sol$par[1],o[,1]*exp(sol$par[2])))

  info<-list(est=sol$par,
             cov=solve(sol$hessian),
             convergence=sol$convergence,
             hresid=resid)

  return(info)
}
#' Make Q-Q plot for residual diagnostics.
#'
#' @param hresid vector of numeric. The heuristic residuals.
#' @param scale numeric. The scale parameter.
#'
#' @examples
#' \dontrun{
#' hresid <- evd::rgumbel(n = 100, loc = 0, scale = 3)
#' qqplot.EVT1(hresid = hresid, scale = 3)
#' }
#' @importFrom graphics abline plot
#' @export
qqplot.EVT1<-function(hresid,scale){
  n<-length(hresid)
  qq<-c(1:n)/(n+1)
  tq<-evd::qgumbel(qq,scale=scale)
  o<-order(hresid)
  plot(tq,hresid[o],xlab="Theoretical quantiles",ylab="Observed quantiles")
  abline(0,1,col="gray")
}

#' Fit RO-logit model and obtain heuristic residuals
#'
#' @param yvar string. Name of outcome variable.
#' @param evar string (or vector of strings). Name of exposure(s).
#' @param cfdr string (or vector of strings). Names of confounder(s). Default is
#'   \code{NULL}.
#' @param emod string (or vector of strings). Name of effect modifier(s).
#'   Default is \code{NULL}.
#' @param svar string. Name of stratum variable. Use \code{NULL} to fit model
#'   without stratification.
#' @param dat \code{data.frame}. Contains all the variables needed for the
#'   analysis.
#' @param method string. Use Efron (\code{"efron"}) or Breslow
#'   (\code{"breslow"}) method for handling ties in the outcome. The default is
#'   \code{"efron"}. See \code{\link{coxph}} for details.
#' @param initial.res.par The initial values of the intercept and log(scale), to
#'   be passed to the \code{optim} function. The default values are set to
#'   \code{c(0, 0)}, yet users are recommended to try a few initial values to
#'   make sure global optimum is reached.
#' @param plot logic. To plot the Q-Q plot for the heuristic residuals. Default
#'   is \code{TRUE}.
#' @param ... Other parameters to be passed to the \code{optim} function for the
#'   second stage analysis.
#' @return Returns a list containing \code{obj} (the RO-Logit model fitted using
#'   \code{coxph}), \code{hresid} (the vector of heuristic residuals),
#'   \code{logscale} (log of scale parameter of the heuristic residuals), and
#'   \code{coefficients} (a data.frame with estimated coefficients before and
#'   after scaling).
#'
#' @references
#' \itemize{
#'  \item{Allison PD, Christakis NA. Logit-models for sets of ranked items.
#'    Sociological Methodology 1994, Vol 24. 1994;24:199-228.}
#'  \item{Beggs S, Cardell S, Hausman J. Assessing the Potential Demand for
#'    Electric Cars. J Econometrics. 1981;17:1-19.}
#'  \item{Tan CS, Støer NC, Chen Y, Andersson M, Ning Y, Wee HL, Khoo EY,
#'    Tai ES, Kao SL, Reilly M. A stratification approach using logit-based
#'    models for confounder adjustment in the study of continuous outcomes.
#'    Statistical methods in medical research. 2017 Jan 1:0962280217747309.}
#'  \item{Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox
#'    Model: Springer New York; 2000.}
#' }
#'
#' @examples
#' # Fit an RO-logit model to determine whether the glycaemic control of
#' # patients differs between medical and surgical wards.
#' data(inpat_bg)
#' # Divide patients into strata based on age, gender, duration of monitoring
#' # episodes, and frequency of daily BG measurements.
#' inpat_bg$group <- paste(inpat_bg$age_group, inpat_bg$sex, inpat_bg$los_group,
#'                         inpat_bg$bg_freq_group, sep = "|")
#' # Fit an RO-logit model with mean BG reading as the outcome and ward as the
#' # exposure:
#' obj <- rologit(yvar = "bg_mean", evar = "ward", svar = "group",
#'                dat = inpat_bg, initial.res.par = c(2, 2))
#' summary(obj)
#' @export
rologit <- function(yvar, evar, cfdr = NULL, emod = NULL, svar, dat,
                    method = "efron", initial.res.par = c(0, 0), plot = TRUE,
                    ...) {
  if (!("data.frame" %in% class(dat))) {
    stop(simpleError("dat should be a data.frame."))
  }
  if (length(evar) > 1) {
    evar_str <- paste(evar, collapse = "+")
  } else {
    evar_str <- evar
  }
  if (!is.null(emod)) {
    emod_str <- paste(sapply(evar, function(e) paste(e, emod, sep = ":")),
                      collapse = "+")
  }
  if (is.null(cfdr)) {
    # Because paste("evar", NULL, sep = "+") gives "evar+"
    if(is.null(emod)){
      form.rologit<-as.formula(paste("~",evar_str,sep=""))
    }else{
      form.rologit<-as.formula(paste("~",paste(evar_str, emod_str,sep="+"),sep=""))
    }
  } else {
    # Paste confounders into a string first, because
    # paste("evar", paste("evar","emod",sep=":"), c("cfdr1","cfdr2"), sep = "+")
    # gives [1] "evar+evar:emod+cfdr1" "evar+evar:emod+cfdr2"
    #
    # paste(c("a"), collapse = "+") gives "a", so this works for single
    # confounder.
    cfdr_str <- paste(cfdr, collapse = "+")
    if(is.null(emod)){
      form.rologit<-as.formula(paste("~",paste(evar_str,cfdr_str,sep="+"),sep=""))
    }else{
      form.rologit<-as.formula(paste("~",paste(evar_str, emod_str,
                                               cfdr_str,sep="+"),sep=""))
    }
  }
  # Design matrix with intercept (for second stage):
  res.design.mat.c<-model.matrix(form.rologit, data=dat)
  # Design matrix without intercept (for fitting the model):
  design.mat<-matrix(res.design.mat.c[,-1],nrow = nrow(dat))
  colnames(design.mat) <- colnames(res.design.mat.c)[-1]
  # Fit RO-Logit model using coxph
  if (!is.null(svar)) {
    svar <- dat[, svar] # Because dat[, NULL] is not NULL
  }
  obj<-rologit.coxph(design.mat=design.mat,y = dat[, yvar], svar = svar,
                     method = method)
  # Get heuristic residuals
  coef_mat <- summary(obj)$coef # Coefficient matrix from coxph
  if (nrow(coef_mat) == 1) {
    o<-matrix(coef_mat[1,c(1,3)],nrow=1)
    colnames(o)<-colnames(coef_mat)[c(1,3)]
    rownames(o)<-rownames(coef_mat)
  } else {
    o<-coef_mat[,c(1,3)]
  }
  oresid.c<-hresid(y = c(dat[,yvar]),svar = svar,
                   design.mat = res.design.mat.c,o = o,
                   initial.res.par = initial.res.par, ...)
  est.logscale.c<-oresid.c$est[length(oresid.c$est)] # log(scale)
  if (plot) {
    # km<-ks.test(oresid.c$resid,"pgumbel",scale=exp(est.logscale.c))
    qqplot.EVT1(hresid = oresid.c$hresid,scale=exp(est.logscale.c))
  }
  ##Gather the matched design results together
  est.scale <- exp(est.logscale.c)
  if(nrow(coef_mat)==1){
    ms.res <- matrix(c(
      coef_mat[1], est.scale * coef_mat[1], # coef and scaled.coef
      coef_mat[2], # exp(coef)
      coef_mat[3], est.scale * coef_mat[3], # se(coef) and se(scaled.coef)
      coef_mat[4:5] # z and p-value
    ), nrow = 1)
    ms.res <- data.frame(variable = rownames(coef_mat), ms.res)
    colnames(ms.res)[-1]<-c("coef","scaled.coef", "exp(coef)","se(coef)",
                            "se(scaled.coef)","z","Pr(>|z|)")
    rownames(ms.res) <- NULL
  }else{
    ms.res <- cbind(
      coef_mat[, 1], est.scale * coef_mat[, 1], # coef and scaled.coef
      coef_mat[, 2], # exp(coef)
      coef_mat[, 3], est.scale * coef_mat[, 3], # se(coef) and se(scaled.coef)
      coef_mat[, 4:5] # z and p-value
    )
    ms.res <- data.frame(variable = rownames(coef_mat), ms.res)
    colnames(ms.res)[-1]<-c("coef","scaled.coef", "exp(coef)","se(coef)",
                            "se(scaled.coef)","z","Pr(>|z|)")
    rownames(ms.res) <- NULL
  }

  output <- list(obj = obj, hresid = oresid.c$hresid,
                 logscale = est.logscale.c, coefficients = ms.res)
  class(output) <- "rologit"
  return(output)
}
#' Summarise RO-Logit Model
#' @description Prints the estimated coefficients of an RO-Logit model.
#'
#' @param object Model object fitted using \code{\link{rologit}}.
#' @param ... Additional arguments affecting the summary produced.
#'
#' @export
summary.rologit <- function(object, ...) {
  stopifnot(inherits(object, "rologit"))
  coef_mat <- object$coefficients
  coef_mat[, -1] <- apply(coef_mat[, -1], 2, function(x) round(x, 3))
  print(coef_mat)
}

#' Inpatient blood glucose data for 2487 patients
#'
#' @description A simulated dataset containing inpatient point-of-care blood
#'   glucose (BG) measurements from 2487 non-critical care adult patients above
#'   40 years old. Data was simulated based on real data.
#'
#' @format A data frame with 2487 rows and 13 variables:
#' \describe{
#'   \item{bg_mean}{Mean BG readings within each episode, in mmol/L.}
#'   \item{bg_sd}{Standard deviation of BG readings within each episode.}
#'   \item{sex}{Gender of patients.}
#'   \item{ward}{Whether each patient is in the medical ward (\code{ward = 0})
#'     or surgical ward (\code{ward = 1}).}
#'   \item{age_group}{Patient group defined based on the median age.}
#'   \item{los_group}{Patient group defined based on the median of duration of
#'     monitoring episode.}
#'   \item{bg_freq_group}{Patient group defined based on the median of daily
#'     BG monitoring frequency.}
#'   \item{age}{Patients' age.}
#'   \item{los}{Patients' duration of monitoring episode.}
#'   \item{bg_freq}{Patients' daily BG monitoring frequency.}
#' }
#' @references Tan CS, Støer NC, Chen Y, Andersson M, Ning Y, Wee HL, Khoo EY,
#'  Tai ES, Kao SL, Reilly M. A stratification approach using logit-based models
#'  for confounder adjustment in the study of continuous outcomes. Statistical
#'  methods in medical research. 2017 Jan 1:0962280217747309.
"inpat_bg"