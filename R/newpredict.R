#' @noRd
#' @keywords internal
# @importFrom stats delete.response
# @importFrom stats terms
#' @import stats
get.pred.vars <- function(object) sapply(attr(stats::delete.response(
  stats::terms(object)), "variables"), as.character)[-1]

#' @noRd
#' @keywords internal
.get.miss.vars<-function(object, newdata) {
  pred.vars <- get.pred.vars(object)
  pred.vars[!(pred.vars %in% colnames(newdata))]
}

#' @noRd
#' @keywords internal
# @importFrom stats contrasts
#' @import stats
getvarn <- function(k, nk) 
  if (is.factor(k)) paste(nk, colnames(stats::contrasts(k)), sep='') else nk

#' @noRd
#' @keywords internal
# @importFrom stats delete.response
# @importFrom stats terms
# @importFrom stats model.frame
#' @import stats
getvarlv <- function(object) {
  mf <- stats::model.frame(object)
  if (!length(dim(mf))) mf <- stats::model.frame(
    stats::delete.response(stats::terms(object)), object$data)
  pred.vars <- get.pred.vars(object)
  pred.vars.lv <- lapply(pred.vars, function(k) 
    getvarn(mf[[which(colnames(mf)==k)]], k))
  names(pred.vars.lv) <- pred.vars
  interactions <- attr(stats::terms(object), 'term.labels')
  interactions <- interactions[grep(':', interactions, fixed=TRUE)]
  rgma <- regmatches(interactions,
                     regexpr(':', interactions, fixed=TRUE), invert=TRUE)
  interactions.lv <- lapply(rgma, function(k) 
    apply(expand.grid(pred.vars.lv[k]), 1, paste, sep='', collapse=':'))
  names(interactions.lv) <- interactions
  c(pred.vars.lv, interactions.lv)
}

#' @noRd
#' @keywords internal
deinteract<-function(n) {
  res <- lapply(n, function(k) 
  tmp <- regmatches(k, regexec(':', k, fixed=TRUE), invert=TRUE)[[1]])
  names(res) <- n
  res
}

#' @noRd
#' @keywords internal
renameIntlev<-function(Vnam,Lnam){
  if (grepl(':',Vnam,fixed = TRUE)) {
    t1<-unname(unlist(deinteract(Vnam)))
    t2<-deinteract(Lnam)
    apply(sapply(t2, function(k) paste(t1,k,sep='')),2,paste,collapse=':')
  } else {
    paste(Vnam,Lnam,sep='')
  }
}

#' @noRd
#' @keywords internal
getRE.lme4 <- function(object, RElist = ranef(object), frame = object@frame){
  tmpRE <- lapply(names(ngrps(object)), function(j) {
    REmean <- as.matrix(RElist[[j]])
    dat <- model.matrix(object, data = frame)
    datm <- dat[,colnames(dat) %in% colnames(REmean), drop=FALSE]
    tmpw <- which(!colnames(REmean) %in% colnames(datm))
    if (length(tmpw)) stop(msgList(24))
    datm <- datm[,colnames(REmean)]
    rownames(REmean) <- unname(renameIntlev(j,rownames(REmean)))
    reform <- as.formula(paste('~',j,'-1'))
    MM <- model.matrix(reform, data = frame)
    MM <- MM[ ,rownames(REmean)]
    if (length(colnames(MM)) != length(rownames(REmean))) stop (msgList(17))
    if (any(colnames(MM) != rownames(REmean))) stop (msgList(17))
    pred1 <- MM %*% REmean
    if (any(dim(datm) != dim(pred1))) stop (msgList(17))
    pred2 <- datm * pred1
    rowSums(pred2)
  })
  sumRE <- 0
  for (k in seq_along(tmpRE)) sumRE <- sumRE + tmpRE[[k]]
  sumRE
}

#' @noRd
#' @keywords internal
simRE.lme4<-function(n.sims, object, frame = object@frame){
  RElist = ranef(object, condVar = TRUE)
  REnames <- names(ngrps(object))
  # j = REnames[2]
  tmpREv <- lapply(REnames, function(j) {
    REcv <- RElist[[j]] 
    REvariance <- attr(REcv, which = "condVar")
    if (!length(REvariance)) REvariance <- attr(REcv, which = "postVar") else
      stop(msgList(26),call. = FALSE)
    REmean <- as.matrix(REcv)
    rm(REcv)
    
    if (ncol(REmean)>1) {
      newMeans <- lapply(seq_len(nrow(REmean)), 
                         function(k) as.matrix(mvtnorm::rmvnorm(
                           n = n.sims, 
                           mean = REmean[k,], 
                           sigma = REvariance[,,k], 
                           method = "chol")))
    } else {
      newMeans <- lapply(seq_len(nrow(REmean)), 
                         function(k) as.matrix(stats::rnorm(
                           n = n.sims,
                           mean = REmean[k,],
                           sd = sqrt(REvariance[,,k])
                         )))
    }
    
    newMeans <- array(as.numeric(unlist(newMeans)), dim=c(n.sims, ncol(REmean), nrow(REmean)))
    dim(newMeans)<-c(n.sims, ncol(REmean), nrow(REmean))
    dimnames(newMeans) <- list(seq_len(n.sims),colnames(REmean),rownames(REmean))
    newMeans
  })  
  
  REsim<-lapply(seq_len(n.sims), function(k) {
    g <- lapply(seq_along(REnames), function(j) {
      raw <- tmpREv[[j]][k,,,drop=FALSE]
      z <- tmpREv[[j]][k,,]
      dim(z) <- dim(raw)[-1]
      dimnames(z) <- dimnames(raw)[-1]
      t(z)
    })
    names(g)<-REnames
    g
  })
  sapply(REsim, function(j) as.matrix(getRE.lme4(RElist=j, object=object, frame=frame)))
}

#' All \%in\% function
#'
#' @param x,y numeric vectors
#' @usage x \%allin\% y
#' @author Maciej J. Danko
#' @keywords internal
'%allin%' <- function(x,y) all(x %in% y)

#' @noRd
#' @keywords internal
get.import.var.lv<-function(object, newdata){
  if (length(newdata)) {
    var.lv <- getvarlv(object)
    import <- unlist(var.lv[which(sapply(deinteract(names(var.lv)), 
                              function(k) k %allin% colnames(newdata)))]) 
    import <- unname(c('(Intercept)', import))
  } else import <- NULL
  import
}

#' @noRd
#' @keywords internal
#' @import stats
.get.avg.missing<-function(object, newdata, coefi){
  import <- get.import.var.lv(object, newdata)
  MM <- MM00 <- if (!isS4(object) && length(object$data)) 
    stats::model.matrix(object, 
                        data=object$data) else stats::model.matrix(object)
  contr <- contr0 <- attr(MM00, 'contrasts')
  contr <- contr[colnames(newdata)]
  attr(MM00,'contrasts')<-contr
  MM00[ , which((colnames(MM00)%in%import))] <- 0
  diff <- MM00 %*% coefi
  
  mf <- stats::model.frame(object)
  if (!length(dim(mf))) mf<-stats::model.frame(
    stats::delete.response(stats::terms(object)), object$data)
  mis.vars <- .get.miss.vars(object, newdata)
  newdata.ext <- cbind(newdata, matrix(0, nrow(newdata), length(mis.vars)))
  colnames(newdata.ext) <- c(colnames(newdata), mis.vars)
  for (k in mis.vars) if (is.factor(mf[[k]])) newdata.ext[[k]] <- 
    factor(levels(mf[[k]])[1],levels=levels(mf[[k]]))  
  for (k in colnames(newdata.ext)) attributes(newdata.ext[[k]]) <- 
    attributes(mf[[k]]) 
  MMn <- stats::model.matrix(stats::delete.response(stats::terms(object)), 
                       data=newdata.ext,
                       contrast.arg=contr)
  MMn[,which(!(colnames(MMn) %in% import))] <- 0
  
  list(diff=diff, model.matrix.incomplete=MMn)
} 

#' @noRd
#' @keywords internal
#' @import stats 
proc.nlme<-function(object){
  terms<-object$terms
  yname <- colnames(
    stats::model.frame(terms,object$data))[attr(terms,"response")]
  MM <- stats::model.matrix(stats::delete.response(terms), object$data)
  oformula <- eval(object$call$fixed)
  mf <- stats::model.frame(delete.response(terms), object$data)
  categorical <- names(attr(MM,'contrasts'))
  term.labels <- attr(terms,"term.labels")
  categorical <- categorical[categorical %in% term.labels]
  pos<-attr(object$terms, "offset")
  if (length(pos)) {
    offset <- mf[,pos]
  } else offset <- 0
  fixed <- object$coefficients$fixed
  fixed.nam <- names(fixed)
  RE <- stats::model.matrix(oformula,object$data)%*%fixed - 
    object$fitted[,object$dims$Q+1] #predict(object, type='link')
  
  RE.nam <- names(object$groups)
  pred.vars <- get.pred.vars(object)
  varcov <- as.matrix(object$varFix)
  linkinv <- identity
  linkfun <- identity
  family <- 'gaussian'
  link <-'identity'
  
  getMM <- function(object, newdata=NULL){
    MM <- stats::model.matrix(stats::delete.response(
      object$terms), object$data)
    if (length(newdata)) {
      MM <- stats::model.matrix(stats::delete.response(object$terms),
                                data=newdata,
                                contrast.arg=attr(MM,'contrasts'))
    } 
    MM
  }
  predFix <- function(object, betas=object$coefficients$fixed, 
                      newdata=NULL, MM=NULL, .offset=0) {
    opos <- attr(terms(object),'offset')
    if (length(opos)) stop(msgList(10),call. = FALSE)
    if (!length(MM)) stop(msgList(11),call.=FALSE)
    if(length(nrow(betas)) && (nrow(betas)>1)) {
      sapply(seq_len(nrow(betas)), 
             function(k) as.vector(MM%*%betas[k,] + .offset))
    } else MM%*%betas + .offset
  }
  lformula<-length(oformula)
  theta <-NULL
  simRE <- function(...) stop(msgList(25), call. = FALSE)
  get.avg.missing <- function(object, newdata, coef=object$coefficients$fixed) 
    .get.avg.missing(object, newdata, coef)
  get.miss.vars <- .get.miss.vars
  list(categorical = categorical,
       term.labels = term.labels,
       offset = offset,
       RE = RE,
       RE.nam = RE.nam,
       frame = mf,
       MM = MM,
       fixed = fixed,
       fixed.nam = fixed.nam,
       resp.nam = yname,
       varcov = varcov,
       family = family,
       link = link,
       linkinv = linkinv,
       linkfun = linkfun,
       lformula = lformula,
       pred.vars = pred.vars,
       predFix = predFix,
       simRE = simRE,
       get.avg.missing = get.avg.missing,
       get.miss.vars = get.miss.vars, 
       getMM = getMM)
}

#' @noRd
#' @keywords internal
#' @import stats mgcv 
proc.gamm<-function(object){
    tmp <- attr(stats::terms(object$gam),"dataClasses")
    categorical <- names(tmp)[which(tmp=='factor')]
    term.labels <- attr(stats::terms(object$gam),"term.labels")
    categorical <- categorical[categorical %in% term.labels]
    mf <- stats::model.frame(object$gam)
    pred.vars <- get.pred.vars(object$gam)
    offset <- stats::model.offset(mf)
    if (!length(offset)){
      Terms <- list(object$gam$pterms)
      offi <-attr(Terms, "offset")
      if (length(offi)) {
        offset <- mf[[names(attr(Terms[[1]], "dataClasses"))[offi]]]
      } else offset <- 0
    }
    #RE <- predict(object$lme,type='link') - mgcv::predict.gam(object$gam,type='link')
    RE <- object$lme$fitted[,object$lme$dims$Q+1] - mgcv::predict.gam(object$gam,type='link')
    # offset is included in gam predict!
    RE <- RE + offset
    RE.nam <- names(object$lme$groups)
    fixed <- stats::coef(object$gam)
    fixed.nam <- names(fixed)
    varcov <- as.matrix(mgcv::vcov.gam(object$gam))
    linkinv <- object$gam$family$linkinv
    linkfun <- object$gam$family$linkfun
    family <- object$gam$family$family
    link <- object$gam$family$link
    if(!length(family)){
      linkinv <- identity
      linkfun <- identity
      family <- 'gaussian'
      link <-'identity'
    }
    yname <- names(attr(object$gam$terms, 
                        "dataClasses"))[attr(object$gam$terms,"response")]
    MM <- mgcv::model.matrix.gam(object$gam, data = mf)
    predFix <- function(object, betas=object$gam$coefficients, 
                        newdata=NULL, MM=NULL, .offset=0) {
      obj <- object$gam
      opos<-attr(stats::terms(obj),'offset')
      if (length(opos)) 
        stop(msgList(1), call.=FALSE)
      if (length(MM)) warning(msgList(2), call.=FALSE)
      if(length(nrow(betas)) && (nrow(betas)>1)) {
        sapply(seq_len(nrow(betas)), function(k) {
          obj$coefficients<-betas[k,]
          mgcv::predict.gam(obj, newdata=newdata, 
                            type='link', se.fit=FALSE) + .offset
        })
      } else {
        obj$coefficients<-betas
        mgcv::predict.gam(obj, newdata=newdata, 
                          type='link', se.fit=FALSE) + .offset
      }
    }
    simRE <- function(...) stop(msgList(25), call. = FALSE)
    get.avg.missing <- function(...) stop(msgList(3),call. = FALSE)
    get.miss.vars <- function(object, newdata) 
      .get.miss.vars(object$gam, newdata)
    getMM <- function(...) stop(msgList(4),call. = FALSE)
    lformula<-length(mgcv::formula.gam(object$gam))
    list(categorical = categorical,
         term.labels = term.labels,
         offset = offset,
         RE = RE,
         RE.nam = RE.nam,
         frame = mf,
         MM = MM,
         fixed = fixed,
         fixed.nam = fixed.nam,
         resp.nam = yname,
         varcov = varcov,
         family = family,
         link = link,
         linkinv = linkinv,
         linkfun = linkfun,
         lformula = lformula,
         pred.vars = pred.vars,
         predFix = predFix,
         simRE = simRE,
         get.avg.missing = get.avg.missing,
         get.miss.vars = get.miss.vars, 
         getMM = getMM)
}

#' @noRd
#' @keywords internal
#' @import stats mgcv
proc.gam <- function(object){
  tmp <- attr(stats::terms(object),"dataClasses")
  categorical <- names(tmp)[which(tmp=='factor')]
  term.labels <- attr(terms(object),"term.labels")
  categorical <- categorical[categorical %in% term.labels]
  mf <- stats::model.frame(object)
  pred.vars <- get.pred.vars(object)
  offset <- stats::model.offset(mf)
  if (!length(offset)){
    Terms <- list(object$pterms)
    offi <-attr(Terms, "offset")
    if (length(offi)) {
      offset <- mf[[names(attr(Terms[[1]], "dataClasses"))[offi]]]
    } else offset <- 0
  }
  RE <- NULL
  RE.nam <- NULL
  fixed <- stats::coef(object)
  fixed.nam <- names(fixed)
  varcov <- as.matrix(mgcv::vcov.gam(object))
  linkinv <- object$family$linkinv
  linkfun <- object$family$linkfun
  family <- object$family$family
  link <- object$family$link
  if(!length(family)){
    linkinv <- identity
    linkfun <- identity
    family <- 'gaussian'
    link <-'identity'
  }
  yname <- names(attr(object$terms, 
                      "dataClasses"))[attr(object$terms,"response")]
  MM <- mgcv::model.matrix.gam(object, data = mf)
  predFix <- function(object, betas=object$coefficients, 
                      newdata=NULL, MM=NULL, .offset=0) {
    obj <- object
    opos<-attr(terms(obj),'offset')
    if (length(opos)) 
      stop(msgList(1), call.=FALSE)
    if (length(MM)) warning(msgList(2), call.=FALSE)
    if(length(nrow(betas)) && (nrow(betas)>1)) {
      sapply(seq_len(nrow(betas)), function(k) {
        obj$coefficients<-betas[k,]
        mgcv::predict.gam(obj, newdata=newdata, 
                          type='link', se.fit=FALSE) + .offset
      })
    } else {
      obj$coefficients<-betas
      mgcv::predict.gam(obj, newdata=newdata, 
                        type='link', se.fit=FALSE) + .offset
    }
  }
  simRE <- function(...) stop(msgList(25), call. = FALSE)
  get.avg.missing <- function(...) stop(msgList(3),call. = FALSE)
  get.miss.vars <- function(object, newdata) 
    .get.miss.vars(object, newdata)
  getMM <- function(...) stop(msgList(4),call. = FALSE)
  lformula<-length(mgcv::formula.gam(object))
  list(categorical = categorical,
       term.labels = term.labels,
       offset = offset,
       RE = RE,
       RE.nam = RE.nam,
       frame = mf,
       MM = MM,
       fixed = fixed,
       fixed.nam = fixed.nam,
       resp.nam = yname,
       varcov = varcov,
       family = family,
       link = link,
       linkinv = linkinv,
       linkfun = linkfun,
       lformula = lformula,
       pred.vars = pred.vars,
       predFix = predFix,
       simRE = simRE,
       get.avg.missing = get.avg.missing,
       get.miss.vars = get.miss.vars, 
       getMM = getMM)
}
  
#' @noRd
#' @keywords internal
#' @import stats gamm4 
proc.gamm4<-function(object){
  tmp <- attr(terms(object$gam),"dataClasses")
  categorical <- names(tmp)[which(tmp=='factor')]
  term.labels <- attr(terms(object$gam),"term.labels")
  categorical <- categorical[categorical %in% term.labels]
  offset<-offsettmp <- try(object$mer@frame[,'(offset)'],silent=TRUE)
  mf <- stats::model.frame(object$gam)
  pred.vars <- get.pred.vars(object$gam)
  if (class(offset)=="try-error") {
    offset <- model.offset(mf)
    if(!length(offset)){
      offc <- grep('offset(',colnames(mf),fixed=TRUE)
      if (length(offc)>1) stop(msgList(5),call.=FALSE)
      if (length(offc)) offset<-mf[,offc] else offset<-0
    }
    offsettmp <- 0
  }
   #RE <- predict(object$mer,type='link')- mgcv::predict.gam(object$gam,type='link')
  RE <- object$mer@resp$eta - mgcv::predict.gam(object$gam,type='link')
  #RE <- stats::fitted.values(object$mer) - mgcv::predict.gam(object$gam,type='link')
  # offset is included in mer predict
  RE <- RE - offsettmp
  RE.nam <- names(getME(object$mer,'flist'))
  fixed <- stats::coef(object$gam)
  fixed.nam <- names(fixed)
  varcov <- as.matrix(mgcv::vcov.gam(object$gam))
  linkinv <- object$gam$family$linkinv
  linkfun <- object$gam$family$linkfun
  family <- object$gam$family$family
  link <- object$gam$family$link
  if(!length(family)){
    linkinv <- identity
    linkfun <- identity
    family <- 'gaussian'
    link <-'identity'
  }
  
  yname <- names(attr(object$gam$terms, 
                      "dataClasses"))[attr(object$gam$terms,"response")]
  MM <- mgcv::model.matrix.gam(object$gam, data = mf)
  predFix <- function(object, betas=object$coefficients, 
                      newdata=NULL, MM=NULL, .offset=0) {
    obj <- object$gam
    opos<-attr(terms(obj),'offset')
    if (length(opos)) stop(msgList(6),call. = FALSE)
    if (length(MM)) warning(msgList(7), call.=FALSE)
    if(length(nrow(betas)) && (nrow(betas)>1)) {
      sapply(seq_len(nrow(betas)), function(k) {
        obj$coefficients<-betas[k,]
        mgcv::predict.gam(obj, newdata=newdata, 
                          type='link', se.fit=FALSE) + .offset
      })
    } else {
      obj$coefficients<-betas
      mgcv::predict.gam(obj, newdata=newdata, 
                        type='link', se.fit=FALSE) + .offset
    }
  }
  lformula<-length(mgcv::formula.gam(object$gam))
  theta<-object$mer@theta
  simRE <- function(...) stop(msgList(27), call. = FALSE)
  get.avg.missing <- function(...) stop(msgList(8),call. = FALSE)
  get.miss.vars <- function(object, newdata) 
    .get.miss.vars(object$gam, newdata)
  getMM <- function(...) stop(msgList(9),call. = FALSE)
  list(categorical = categorical,
       term.labels = term.labels,
       offset = offset,
       RE = RE,
       RE.nam = RE.nam,
       frame = mf,
       MM = MM,
       fixed = fixed,
       fixed.nam = fixed.nam,
       resp.nam = yname,
       varcov = varcov,
       family = family,
       link = link,
       linkinv = linkinv,
       linkfun = linkfun,
       lformula = lformula,
       pred.vars = pred.vars,
       predFix = predFix,
       simRE = simRE,
       get.avg.missing = get.avg.missing,
       get.miss.vars = get.miss.vars, 
       getMM = getMM)
}

#' @noRd
#' @keywords internal
#' @import stats lme4
proc.lme4<-function(object){
  contrast <- attr(getME(object,'X'),'contrasts')
  categorical <-names(contrast)
  term.labels <- attr(terms(object),"term.labels")
  categorical <- categorical[categorical %in% term.labels]
  offset <- object@resp$offset
  RE <- if ('eta' %in% names(object@resp)) 
    object@resp$eta-getME(object,'X')%*%fixef(object)-offset else 
      object@resp$mu-getME(object,'X')%*%fixef(object)-offset
  RE.nam <- names(getME(object,'flist'))
  mf<-object@frame
  pred.vars <- get.pred.vars(object)
  fixed <- fixef(object)
  fixed.nam <- names(fixed)
  varcov <- as.matrix(vcov(object))
  if (class(object@resp)!="lmerResp") {
    linkinv <- object@resp$family$linkinv
    linkfun <- object@resp$family$linkfun
    family <- object@resp$family$family
    link <- object@resp$family$link
  } else {
    linkinv <- identity
    linkfun <- identity
    family <- 'gaussian'
    link <- 'identity'
  }
  # yname <- deparse(object@call$formula[[2]])
  yname <- names(mf)[attr(terms,"response")]
  # zrobic funkcje dla averaging
  # dziala tylko lme i lme4
  MM <- getME(object,'X')
  getMM <- function(object, newdata=NULL){
    MM <- getME(object,'X')
    if (length(newdata)) {
      MM <- model.matrix(delete.response(terms(object)),
                         data=newdata,
                         contrast.arg=attr(MM,'contrasts'))
    }
    MM
  }
  predFix <- function(object, betas=fixef(object), 
                      newdata=NULL, MM=NULL, .offset=0) {
    opos <- attr(terms(object),'offset')
    if (length(opos)) stop(msgList(10),call. = FALSE)
    if (!length(MM)) stop(msgList(11),call.=FALSE)
    if(length(nrow(betas)) && (nrow(betas)>1)) {
      sapply(seq_len(nrow(betas)), 
             function(k) as.vector(MM%*%betas[k,] + .offset))
    } else MM%*%betas + .offset
  }
  simRE <- simRE.lme4
  lformula<-length(object@call$formula)
  get.avg.missing <- function(object, newdata, coef=fixef(object)) 
    .get.avg.missing(object, newdata, coef)
  get.miss.vars <- .get.miss.vars
  list(categorical = categorical,
       term.labels = term.labels,
       offset = offset,
       RE = RE,
       RE.nam = RE.nam,
       frame = mf,
       MM = MM,
       fixed = fixed,
       fixed.nam = fixed.nam,
       resp.nam = yname,
       varcov = varcov,
       family = family,
       link = link,
       linkinv = linkinv,
       linkfun = linkfun,
       lformula = lformula,
       pred.vars = pred.vars,
       predFix = predFix,
       simRE = simRE,
       get.avg.missing = get.avg.missing,
       get.miss.vars = get.miss.vars, 
       getMM = getMM)
}

#' @noRd
#' @keywords internal
#' @import stats MASS
proc.glmnb<-function(object){
  MM <- stats::model.matrix(object)
  categorical <- names(attr(MM,'contrasts'))
  terms <- stats::terms(object)
  term.labels <- attr(terms,"term.labels")
  categorical <- categorical[categorical %in% term.labels]
  pos<-attr(stats::terms(object), "offset")
  if (length(pos)) {
    offset <- stats::model.frame(object$formula,object$data)[,pos]
  } else offset <- 0
  RE <- NULL
  RE.nam <- NULL
  mf <- stats::model.frame(object)
  pred.vars <- get.pred.vars(object)
  fixed <- stats::coef(object)
  fixed.nam <- names(fixed)
  varcov <- as.matrix (vcov(object)) 
  if (inherits(object,c('glm','negbin'))){
    linkfun<-object$family$linkfun
    linkinv<-object$family$linkinv
    family <- object$family$family
    link <- object$family$link
  } else {
    linkinv <- identity
    linkfun <- identity
    family <- 'gaussian'
    link <-'identity'
  }
  yname <- names(mf)[attr(terms,"response")]
  getMM <- function(object, newdata=NULL){
    MM <- stats::model.matrix(object)
    if (length(newdata)) {
      MM <- model.matrix(stats::delete.response(stats::terms(object)),
                         data=newdata,
                         contrast.arg=attr(MM,'contrasts'))
    }
    MM
  }
  predFix <- function(object, betas=stats::coef(object), 
                      newdata=NULL, MM=NULL, .offset=0) {
    opos <- attr(stats::terms(object),'offset')
    if (length(opos)) stop(msgList(12),call. = FALSE)
    if (!length(MM)) stop(msgList(11),call.=FALSE)
    if(length(nrow(betas)) && (nrow(betas)>1)) {
      sapply(seq_len(nrow(betas)), 
             function(k) as.vector(MM%*%betas[k,] + .offset))
    } else MM%*%betas + .offset
  }
  lformula<-length(stats::formula(object))
  theta<-NULL
  simRE <- function(...) stop(msgList(25), call. = FALSE)
  get.avg.missing <- function(object, newdata, coefi=stats::coef(object)) 
    .get.avg.missing(object, newdata, coefi)
  get.miss.vars <- .get.miss.vars
  list(categorical = categorical,
       term.labels = term.labels,
       offset = offset,
       RE = RE,
       RE.nam = RE.nam,
       frame = mf,
       MM = MM,
       fixed = fixed,
       fixed.nam = fixed.nam,
       resp.nam = yname,
       varcov = varcov,
       family = family,
       link = link,
       linkinv = linkinv,
       linkfun = linkfun,
       lformula = lformula,
       pred.vars = pred.vars,
       predFix = predFix,
       simRE = simRE,
       get.avg.missing = get.avg.missing,
       get.miss.vars = get.miss.vars, 
       getMM = getMM)
}

#' @noRd
#' @keywords internal
collect.model.info<-function(object){
  if (all(c('gam','lme') %in% names(object)) && 
      all(inherits(object, c('gamm','list'), which = TRUE)>0)) {
    proc.gamm(object)
  } else if (inherits(object, c('gam'))) {
    proc.gam(object)
  } else if (all(c('gam','mer') %in% names(object)) && 
             inherits(object, c('list'))) {
    proc.gamm4(object)
  } else if (inherits(object,c('glmerMod','merMod'))){
    proc.lme4(object)
  } else if (inherits(object,'lme')){
    proc.nlme(object)
  } else if (inherits(object,c('lm','glm','negbin'))){
    proc.glmnb(object)
  } else stop(msgList(13),call. = FALSE)
} 

#' @noRd
#' @keywords internal
.predict.classify<-function(objsup, variable, newdata = NULL){
  if (length(newdata)){
    # leave only categorical
    indc <- which(colnames(newdata) %in% objsup$categorical)
    newdata <-droplevels.data.frame(newdata[,indc,drop = FALSE])
  }
  if (length(newdata)) {
    data <- objsup$frame
    data <- data[, colnames(data) %in% colnames(newdata),drop = FALSE]
    colind <- match(colnames(newdata),colnames(data))
    nd_codes <- apply(newdata, 1, paste, collapse='_', sep='_')
    sdata <- data[,colind,drop=FALSE]
    sdata_codes <- apply(sdata, 1, paste, collapse='_', sep='_')
    
    raw <- lapply(nd_codes,function(k) variable[which(sdata_codes==k)])
    names(raw) <- nd_codes
    attr(raw,'codes')<-apply(newdata, 1, 
      function(k) apply(rbind(colnames(newdata),k),2,paste,collapse=''))
    attr(raw,'data')<-newdata
  } else {
    raw <- variable 
  }
  return(raw)
}

#' @noRd
#' @keywords internal
.rebuild.newdata<-function(objsup, newdata){
  m <- match(c(objsup$RE.nam,objsup$resp.nam), colnames(newdata))
  m <- m[which(!is.na(m))]
  if (length(m)) {
    newdata <-newdata[,-m,drop=FALSE]
    warning(msgList(14),call.=FALSE)
  }
  
  # dropping unused
  Xnam <- names(objsup$frame)
  inddro <- which(is.na(match(colnames(newdata), Xnam)))
  if (length(inddro)) newdata <- newdata[,-inddro,drop=FALSE]
  
  data <- objsup$frame
  data <- data[, colnames(data) %in% colnames(newdata),drop = FALSE]
  colind <- which((colnames(newdata) %in% objsup$categorical) & 
                    (colnames(newdata) %in% colnames(data)))
  if (length(colind)) {
    ind <- sapply(seq_len(nrow(newdata)),
             function(k) all(sapply(colnames(newdata)[colind], 
               function(j) newdata[k, j] %in% data[,j])))
    if (any(!ind)) stop(msgList(15),call.=FALSE)
  }
  
  attribu <-lapply(objsup$frame[names(newdata)],attributes)
  for (i in names(attribu)) if (length(attribu[[i]]$levels)) 
    newdata[[i]] <- factor(newdata[[i]], attribu[[i]]$levels)
  newdata
}

#' @noRd
#' @keywords internal
.fixmiss<-function(fixmiss, object, objsup, betas, 
                       newdata, sim.missing, method){
  n<-nrow(betas)
  if (sim.missing) {
    if (method == 'at.each.cat') {
      nd <- newdata 
    } else if (method == 'overall') {
      nd <- NULL 
    } else stop(msgList(16), call. = FALSE)
    
    if (length(nd)){
      # leave only categorical
      indc <- which(colnames(nd) %in% objsup$categorical)
      nd <- droplevels.data.frame(nd[,indc,drop = FALSE])
    }
    if (length(nd)) {
      data <- objsup$frame
      data <- data[, colnames(data) %in% colnames(nd),drop = FALSE]
      colind <- match(colnames(nd),colnames(data))
      nd_codes <- apply(nd, 1, paste, collapse='_', sep='_')
      sdata <- data[,colind,drop=FALSE]
      sdata_codes <- apply(sdata, 1, paste, collapse='_', sep='_')
    } 
    if (length(nd)) {
      gg <- lapply(seq_len(n), function(j) {
        z <- objsup$get.avg.missing(object, newdata, 
                                    as.vector(betas[j,]))
        lapply(nd_codes,function(k) z$diff[which(sdata_codes==k)])
      })
    } else {
      gg <- lapply(seq_len(n), function(j) 
        list(objsup$get.avg.missing(object, newdata, 
                                    as.vector(betas[j,]))$diff))
    }
    fixmissFunc <- gg #function(k) gg[[k]]
  } else fixmissFunc <- list(fixmiss) #<-function(k) fixmiss   
  fixmissFunc
}


#' @noRd
#' @keywords internal
.reff <- function(reff, object, objsup, n, newdata, sim.RE, method){
  if (sim.RE) {
    if (method == 'at.each.cat') {
      nd <- newdata 
    } else if (method == 'overall') {
      nd <- NULL 
    } else stop(msgList(16), call. = FALSE)
    
    if (length(nd)){
      indc <- which(colnames(nd) %in% objsup$categorical)
      nd <- droplevels.data.frame(nd[,indc,drop = FALSE])
    }
    if (length(nd)) {
      data <- objsup$frame
      data <- data[, colnames(data) %in% colnames(nd),drop = FALSE]
      colind <- match(colnames(nd),colnames(data))
      nd_codes <- apply(nd, 1, paste, collapse='_', sep='_')
      sdata <- data[,colind,drop=FALSE]
      sdata_codes <- apply(sdata, 1, paste, collapse='_', sep='_')
    } 
    
    REsim <- objsup$simRE(n, object)
    
    if (length(nd)) {
      gg <- lapply(seq_len(n), function(j) 
        lapply(nd_codes,function(k) REsim[which(sdata_codes==k) ,j]))
    } else {
      gg <- lapply(seq_len(n), function(j) list(REsim[,j]))
    }
    REFunc <- gg 
  } else REFunc <- list(reff) #<-function(k) fixmiss   
  REFunc
}

#' @noRd
#' @keywords internal
.make.predictions<-function(object, objsup, betas, newdata, MM, coffs, as.rate,
                            reff, fixmiss, linkinv, linkfun, sim.missing, 
                            method, sim.RE){
  
  n<-nrow(betas)
  
  fixmiss <- .fixmiss(fixmiss, object, objsup, betas, 
           newdata, sim.missing, method)
  
  reff <- .reff(reff, object, objsup, n, newdata, sim.RE, method)
  
  fixmissFunc <- function(n,k) {
    tmp <- if (length(fixmiss) == 1) fixmiss[[1]] else fixmiss[[n]]
    if (length(tmp) == 1) tmp[[1]] else
      if (length(tmp) > 1) tmp[[k]] else 0
  }
  
  reffFunc <- function(n, k) {
    tmp <- if (length(reff) == 1) reff[[1]] else reff[[n]]
    if (length(tmp) == 1) tmp[[1]] else
      if (length(tmp) > 1) tmp[[k]] else 0
  }
  
  # reffFunc <- function(k) if (length(reff)==1) reff[[1]] else
  #   if (length(reff)>1) reff[[k]] else 0
  
  offsetFunc <- function(k) if (length(coffs) == 1) coffs[[1]] else
    if (length(coffs) >1 ) coffs[[k]] else 0
  
  cexpo <- if (as.rate) {
    if (length(coffs) == 1) list(linkinv(coffs[[1]])) else
    if (length(coffs) > 1) lapply(coffs, linkinv) else list(1)
  } else list(1)
  
  exposuresFunc<-function(k) if (length(cexpo) > 1) cexpo[[k]] else cexpo[[1]]
  
  pred_fixed_base <- sapply(seq_len(n),function(j) 
    objsup$predFix(object, as.vector(betas[j, ]), newdata, MM, 0))
  
  K <- seq_len(nrow(pred_fixed_base))
  J <- seq_len(n)
  
  pred_fixed_biased <- sapply(J, function(j) 
    sapply(K, function(k) 
      mean(linkinv(fixmissFunc(j, k) + 
                     pred_fixed_base[k, j] + 
                     offsetFunc(k))) / mean(exposuresFunc(k))))
  
  pred_fixed_unbiased <- sapply(J, function(j) 
    sapply(K, function(k) 
      mean(linkinv(reffFunc(j, k) + 
                     fixmissFunc(j, k) + 
                     pred_fixed_base[k, j] + 
                     offsetFunc(k)))/ mean(exposuresFunc(k))))
  
  list(biased = pred_fixed_biased,
       unbiased = pred_fixed_unbiased)
} 
  
#' Unbiased predictions for fixed variables
#'
#' @param object a object of class inherited from: \code{\link{glmer}{lme4}}, 
#' \code{\link{gamm}{mgcv}}, \code{\link{gam}{mgcv}}, \code{\link{gamm4}{gamm4}}, \code{\link{glmer.nb}{lme4}}, \code{\link{lmer}{lme4}}, \code{\link{lme}{nlme}}, 
#' \code{\link{glm}{stats}}, \code{\link{lm}{stats}}, or \code{\link{glm.nb}{MASS}}.
#' @param newdata an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param type the type of prediction required, either "\code{link}" or "\code{response}".
#' The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. 
#' see also \code{\link{predict.glm}{stats}}.
#' @param ci.fit a logical indicting if to bootstrap percentile confidence intervals.
#' @param alpha,n.sims a significance level and a number of repetitions for the bootstrap method.
#' @param average.missing a logical indicating if to perform population averaging over the variables missing from the \code{newdata}.
#' @param sim.missing a logical indicating if to simulate uncertenity of the missing variables via the bootstrap method. Miningfull only when \code{average.missing} is \code{TRUE}.
#' @param sim.RE a logical indicating if to simulate uncertenity of the random effects via the bootstrap method. 
#' @param method a method of averaging over random effects. 
#' The default "\code{at.each.cat}" estimate avarage random effects at each cobination of levels for categorical variable(s) present in the \code{newdata}; 
#' the alternative "\code{overall}" estimate average random effects for total population.
#' @param use.offset a logical indicating if to include averaged offset into model predictions.
#' @param as.rate a logical indicating if to return predicitions as rates. Setting this \code{TRUE} makes sense only for count models (e.g., Poisson, neg.bin). 
#' @param return.median a logical indicating if to return median for bootstrap percentile method.
#' @param return.boot a logical indicating if to return bootstraped matrices.
#' 
#' @return a list with model predicitions with the following components:
#'  \item{fit}{ a list with \code{biased} (random effects excluded) and \code{unbiased}(averaged random effects included) predicitions.}
#'  \item{CI.fit}{ a optional list with bootstraped confidence intervals for \code{biased} and \code{unbiased} predicitions.}
#'  \item{boot}{ a optional list with bootstraped matrices for \code{biased} and \code{unbiased} predicitions.}
#'  \item{offset}{ the offset used for model predictions.}
#'  
#' @importFrom MASS mvrnorm
#' @importFrom stats quantile
#' @export
#' @author Maciej J. Danko
#' @examples
#' #############################################################################
#' # EXAMPLE 1
#' #############################################################################
#' 
#' set.seed(3)
#' data1 <- sim_glmer_poisson(formula = ~ X1 * X2, theta = 0.75,
#'                          coef= c(5, 0.5, 0.7, 0.5, -0.7, 0.2, 0.3, 0.3, 0.2),
#'                          n.levels=c(3,3), n.ID = 100, n.pop=1e3)
#' 
#' fit1 <- lme4::glmer(Y~X1*X2+(1|ID),family=poisson(), data=data1)
#' 
#' DATA <- with(data1, tapply(Y, data.frame(X1=X1,X2=X2), mean))
#' posX <- as.vector(barplot(DATA,beside=TRUE,ylim=c(0,1050),
#'                           ylab='Mean counts (response scale)'))
#' text(posX,-40,rep(rownames(DATA),3),xpd=TRUE)
#' axis(1,at=c(-5,100),labels = FALSE)
#' 
#' newdata1 <- expand.grid(X1=levels(data1$X1), X2=levels(data1$X2))
#' 
#' # calculate the biased and the unbiased predictions 
#' # on the response scale
#' cYA <- fixPredict(object=fit1,
#'                   newdata=as.data.frame(newdata1),
#'                   type='response',
#'                   ci.fit = TRUE,
#'                   method = 'at.each.cat')
#' 
#' # using overall method
#' cYB <- fixPredict(object=fit1,
#'                   newdata=as.data.frame(newdata1),
#'                   type='response',
#'                   ci.fit = TRUE,
#'                   method = 'overall')
#' 
#' # plot the predicited responses
#' lines(posX-0.3,cYA$fit$biased,pch=20,type='p',col='red',cex=1.5)
#' lines(posX+0.3,cYA$fit$unbiased,pch=20,type='p',
#'       col=rgb(0,180,0,maxColorValue = 255),cex=1.5)
#' lines(posX,cYB$fit$unbiased,pch=20,type='p',col='darkorange',cex=1.5)
#' 
#' # plot the confidence intervals (the bootstrap percentile method)
#' for (i in seq_along(posX)) lines(c(posX[i],posX[i])+0.3,
#'                                  cYA$CI.fit$unbiased[i,], 
#'                                  col=rgb(0,180,0,maxColorValue = 255))
#' for (i in seq_along(posX)) lines(c(posX[i],posX[i])-0.3,
#'                                  cYA$CI.fit$biased[i,], 
#'                                  col='red')
#' for (i in seq_along(posX)) lines(c(posX[i],posX[i]),
#'                                  cYB$CI.fit$unbiased[i,], 
#'                                  col='darkorange')
#' 
#' legend('topleft',c('Simulated data','RE excluded',
#'                    'RE averaged',
#'                    'RE category-averaged'),
#'        pch=c(NA,19,19,19),pt.cex=c(1,1.5,1.5,1.5),
#'        fill=c(gray.colors(3)[2],NA,NA,NA),
#'        col=c(NA,2,'darkorange',rgb(0,180,0,maxColorValue = 255)), 
#'        bty='n', border=c('black',NA,NA,NA))
#' \donttest{        
#' #############################################################################
#' # EXAMPLE 2
#' #############################################################################
#' data2 <- sim_glmer_data(formula = ~ X1, theta = 0.75,
#'                         coef= c(3, 0.5, -0.7),
#'                         coef.c = 5,
#'                         func.c = function(C1, coef.c) 
#'                           abs(C1)^1.5 * coef.c - 1.5,
#'                         n.levels=3, n.ID = 100, n.pop=1e3)
#' 
#' # fitting the model
#' fit2 <- gamm4::gamm4(Y~ X1 + s(C1), random = ~(1|ID), 
#'                      family = poisson(), data = data2)
#' 
#' # making the predictions
#' newdata2 <- expand.grid(C1=seq(-0.5,0.5,0.05), X1=levels(data1$X1))
#' cY2 <- fixPredict(fit2, newdata2, type = 'response', ci.fit = TRUE)
#' 
#' # plotting the results
#' ina <- which(newdata2$X1=='a'); inb <- which(newdata2$X1=='b')
#' inc <- which(newdata2$X1=='c')
#' xa <- newdata2$C1[ina]; xb <- newdata2$C1[inb]; xc <- newdata2$C1[inc]
#' 
#' CIplot.ci(xa, cY2$fit$unbiased[ina], cY2$CI.fit$unbiased[ina,],
#'           density=40, angle=0,first = TRUE,ylim=c(0,60),col='blue',lwd=4,
#'           xaxt='s', yaxt='s', ylab='Mean counts (response scale)', xlab='C1')
#' CIplot.ci(xa, cY2$fit$unbiased[inb], cY2$CI.fit$unbiased[inb,],
#'           density=40, angle=0,col='red',lwd=4)
#' CIplot.ci(xa, cY2$fit$unbiased[inc], cY2$CI.fit$unbiased[inc,],
#'           density=40, angle=0,col=rgb(0,180,0,maxColorValue = 255),lwd=4)
#' CIplot.ci(xa, cY2$fit$biased[ina], cY2$CI.fit$biased[ina,],
#'           density=20, angle=90,col='blue',lty=2)
#' CIplot.ci(xa, cY2$fit$biased[inb], cY2$CI.fit$biased[inb,],
#'           density=20, angle=90,col='red',lty=2)
#' CIplot.ci(xa, cY2$fit$biased[inc], cY2$CI.fit$biased[inc,],
#'           density=20, angle=90,col=rgb(0,180,0,maxColorValue = 255),lty=2)
#' 
#' legend('top',c('biased','unbiased','a','b','c'), lty=c(2,1,1,1,1),
#'        lwd=c(1,4,2,2,2),col=c(1,1,'blue','red',
#'                               rgb(0,180,0,maxColorValue = 255)),bty='n')
#'                               
#' #############################################################################
#' # EXAMPLE 3
#' #############################################################################
#' data3 <- sim_glmer_data(formula = ~ X1, theta = 1.1,
#'                         coef= c(0.5, 1.5, -1.25),
#'                         coef.c = 2, s=0.3,
#'                         func.c = function(C1, coef.c) C1* coef.c,
#'                         n.levels=3, n.ID = 100, n.pop=1e3,
#'                         family='binomial')
#' 
#' ina <- which(data3$X1=='a'); inb <- which(data3$X1=='b')
#' inc <- which(data3$X1=='c')
#' plot(data3$C1[ina],binomial()$linkinv(data3$linkY)[ina],
#'      ylim=c(0,1), cex=0.2, col = 'blue')
#' lines(data3$C1[inb],binomial()$linkinv(data3$linkY)[inb],
#'       type='p', cex=0.2, col = 'red')
#' lines(data3$C1[inc],binomial()$linkinv(data3$linkY)[inc],
#'       type='p', cex=0.2, col = rgb(0,180,0,maxColorValue = 255))
#' 
#' fit3 <- glmer(Y ~ X1 + C1 + (1|ID), family=binomial(), data=data3, nAGQ = 20)
#' cY3 <- fixPredict(fit3, newdata3, type = 'response', ci.fit = TRUE)
#' 
#' newdata3 <- expand.grid(C1=seq(-0.5,0.5,0.05), X1=levels(data3$X1))
#' 
#' ina <- which(newdata3$X1=='a'); inb <- which(newdata3$X1=='b')
#' inc <- which(newdata3$X1=='c')
#' xa <- newdata3$C1[ina]; xb <- newdata3$C1[inb]; xc <- newdata3$C1[inc]
#' 
#' CIplot.ci(xa, cY3$fit$unbiased[ina], cY3$CI.fit$unbiased[ina,],
#'           density=40, angle=0,first = TRUE,ylim=c(0,1),col='blue',lwd=4,
#'           xaxt='s', yaxt='s', ylab='Probabilyty (response scale)', xlab='C1')
#' CIplot.ci(xa, cY3$fit$unbiased[inb], cY3$CI.fit$unbiased[inb,],
#'           density=40, angle=0,col='red',lwd=4)
#' CIplot.ci(xa, cY3$fit$unbiased[inc], cY3$CI.fit$unbiased[inc,],
#'           density=40, angle=0,col=rgb(0,180,0,maxColorValue = 255),lwd=4)
#' CIplot.ci(xa, cY3$fit$biased[ina], cY3$CI.fit$biased[ina,],
#'           density=20, angle=90,col='blue',lty=2)
#' CIplot.ci(xa, cY3$fit$biased[inb], cY3$CI.fit$biased[inb,],
#'           density=20, angle=90,col='red',lty=2)
#' CIplot.ci(xa, cY3$fit$biased[inc], cY3$CI.fit$biased[inc,],
#'           density=20, angle=90,col=rgb(0,180,0,maxColorValue = 255),lty=2)
#' 
#' legend('top',c('Biased pred. response', 'Unbiased pred. response',
#'                'a','b','c'), lty=c(2,1,1,1,1), lwd=c(1,4,2,2,2), 
#'        col=c(1,1,'blue','red', rgb(0,180,0,maxColorValue = 255)), bty='n')
#'                               
#' #############################################################################
#' # EXAMPLE 4
#' #############################################################################
#' #dontrun
#' #cat(1)
#' }
fixPredict<-function(object,
                     newdata,
                     type = c('link','response'),
                     ci.fit = TRUE,
                     alpha = 0.95,
                     n.sims = 1e3,
                     average.missing = FALSE, 
                     sim.missing = TRUE,
                     sim.RE = FALSE,
                     method = c('at.each.cat','overall'),
                     use.offset = TRUE,
                     as.rate = FALSE, 
                     return.median = FALSE,
                     return.boot = FALSE){
  
  objsup <- collect.model.info(object)
  method <- tolower(method[1])
  type <- tolower(type[1])
  CIlo <- (1 - alpha)/2
  boot <- NULL
  fixmiss <- NULL
  # if (include.rnd.error && simtheta) thetas <- rtheta 
  if (objsup$lformula!=3) stop(msgList(20),call.=FALSE)
  if (!missing(newdata) && length(newdata))
    newdata <- .rebuild.newdata(objsup, newdata)
  if (!length(objsup$get.miss.vars(object, newdata))){
    average.missing <- FALSE
    sim.missing <- FALSE
  }
  
  if (type=='link') {
    linkinv<- identity 
    linkfun<- identity 
    if (as.rate) stop(msgList(18),call.=FALSE)
  } else if (type=='response') {
    linkinv<- objsup$linkinv 
    linkfun<- objsup$linkfun
  } else stop(msgList(19))
  if (as.rate && objsup$link!='log') stop(msgList(23),call.=FALSE)
  if (!average.missing) {
    if (length(newdata)) {
      mv <- objsup$get.miss.vars(object, newdata)
      if (length(mv)) 
        stop(paste(msgList(21),
                   paste(mv,collapse=', '),
                   msgList(22),
                   sep=''), call.=FALSE)
    }
    MM <- try(objsup$getMM(object, newdata), silent=TRUE)
    if (class(MM)=="try-error") MM<-NULL
  } else MM <- NULL
  
  N <- nrow(objsup$frame)
  offset <- objsup$offset
  if (length(offset) && use.offset) {
    if (length(offset)==1) {
      offset <- rep(offset, N)
    } else if (length(offset)!=N) stop(msgList(17), call.=TRUE) 
  } else {
    offset <- rep(0, N)
  }
  objsup$offset <- offset
  
  #if (as.rate && any(offset<=0)) stop(msgList(24),call.=FALSE)
  
  if (length(newdata)){
    if (method == 'at.each.cat') {
      nd <- newdata 
    } else if (method == 'overall') {
      nd <- NULL 
    } else stop(msgList(16), call. = FALSE)
    
    if (length(objsup$RE)){
      reff<- .predict.classify(objsup, objsup$RE, nd)
      reff.data<-attr(reff,'data')
      if (!length(reff.data)) reff<-list(reff)
    } else {
      reff <- NULL
      reff.data <- NULL
    }
    
    if (average.missing) {
      z <- objsup$get.avg.missing(object, newdata)
      MM <- z$model.matrix.incomplete
      fixmiss <- .predict.classify(objsup,z$diff,nd)
      if (!length(attr(fixmiss,'data'))) fixmiss<-list(fixmiss)
    }
    
    if (!length(offset) || sum(abs(offset))==0 || !use.offset) {
      coffs <- NULL
      coffs.data<- NULL
    } else {
      coffs<- .predict.classify(objsup, offset, nd)
      coffs.data<-attr(coffs,'data')
      if (!length(coffs.data)) coffs<-list(coffs)
    }
    
  } else {
    if (length(objsup$RE)){
      reff<- .predict.classify(objsup, objsup$RE)
      reff.data <- NULL
    } else {
      reff <- NULL
      reff.data <- NULL
    }
    coffs <- list(offset)
  }
  
  betas <- t(as.matrix(objsup$fixed)) 
  fit <- .make.predictions(object, objsup, betas, newdata, MM, coffs, as.rate, 
                           reff, fixmiss, linkinv, linkfun, FALSE, method, sim.RE)
  
  if (ci.fit && (n.sims>1)){
    betas <- MASS::mvrnorm(n = n.sims, mu = objsup$fixed, Sigma = objsup$varcov)
    fitCI <- .make.predictions(object, objsup, betas, newdata, MM, coffs, 
                               as.rate, reff, fixmiss, linkinv, linkfun, 
                               sim.missing & average.missing, method, sim.RE)    
    
    if (return.median) {
      CI_fixed_biased <- t(apply(fitCI$biased, 1, stats::quantile, 
                                 probs=c(CIlo, 1-CIlo, 0.5)))
      CI_fixed_unbiased <- t(apply(fitCI$unbiased, 1, stats::quantile,
                                   probs=c(CIlo, 1-CIlo, 0.5)))
    } else {
      CI_fixed_biased <- t(apply(fitCI$biased, 1, stats::quantile,
                                 probs=c(CIlo, 1-CIlo)))
      CI_fixed_unbiased <- t(apply(fitCI$unbiased, 1, stats::quantile,
                                   probs=c(CIlo, 1-CIlo)))
    }
    CI.fit <- list(biased=CI_fixed_biased,
                   unbiased=CI_fixed_unbiased)
    
    if (return.boot) boot <- list(biased = CI_fixed_biased,
                                  unbiased = CI_fixed_unbiased)
  } else CI.fit <- NULL
  
  return(list(fit = list(biased = fit$biased, 
                         unbiased = fit$unbiased),
              flags = list(type = type,
                           ci.fit = ci.fit,
                           alpha = alpha,
                           n.sims = n.sims,
                           average.missing = average.missing, 
                           sim.missing = sim.missing,
                           method = method,
                           use.offset = use.offset,
                           as.rate = as.rate, 
                           return.median = return.median,
                           return.boot = return.boot),
              CI.fit = CI.fit,
              boot = boot,
              offset = offset))
}

