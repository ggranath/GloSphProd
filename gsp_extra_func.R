##########################################################################################
# Extra functions to analyse data and produce the results presetned in                            #
# " Environmental drivers of Sphagnum growth in peatlands across the Holarctic region"   #
#                                                                                        #
# Contact: Gustaf.Granath@gmail.com                                                      #
##########################################################################################


# R2 for random-slop models####
# The function is using code from the MuMIn package
# https://github.com/rforge/mumin
r2.fnc <- function(model, null.model) {
  # Needed functions taken from MuMIn r.squaredGLMM function
  matmultdiag <-
    function(x, y, ty = t(y)) {
      if(ncol(x) != ncol(ty)) stop('non-conformable arguments')
      if(nrow(x) != nrow(ty)) stop('result is not a square matrix')
      return(rowSums(x * ty))
    }
  varRESum <- function(vc, X) {
    n <- nrow(X)
    sum(sapply(vc, function(sig) {
      mm1 <-  X[, rownames(sig), drop = FALSE]
      sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
    }))
  }
  
res <- list() 

# Loop over model and null model. Save variances
# in a list

for (i in 1:2) {
  if(i==1) {object <- model} 
  if(i==2) {object <- null.model}
fam <- family(object)

# fixed effect variance componeent
fe <- fixef(object)
ok <- !is.na(fe)
fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
varFE <- var(fitted)

# variances
vc <- VarCorr(object)

# get random variance and residual variance
mmRE <- do.call("cbind", model.matrix(object, type = "randomListRaw"))
varRE <- varRESum(vc, mmRE[, unique(colnames(mmRE)), drop = FALSE])
se = attr(vc, "sc")^2
res[i] <- list(c(varFE = varFE,varRE = varRE, se = se, vtot = sum(varFE, varRE)))
  }

# print out site and wihin site R2, and marginal and conditional R2

return(list(r2MC = with(data.frame(t(res[[1]])), matrix(c(varFE, vtot) / (vtot + rep(se, each = 2L)),
       ncol = 2L, byrow = TRUE, dimnames = list(names(se), c("R2m", "R2c")))),
    r2bSite = 1-(res[[1]]["varRE"]/res[[2]]["varRE"]),
    r2wSite = 1-(res[[1]]["se"]/res[[2]]["se"]),
    No_sites = NROW(ranef(object)$Site), #Groups
    No_sample = NROW(object@frame))) #observations
}

# Resamp R2 function to get R2 uncertainty by a parametric bootstrap appraoch ####
r2.resamp <- function(mod,  sims=10, random = Site, only.whole=TRUE) {
  #help function
  #https://stackoverflow.com/questions/23381616/r-update-function-how-to-drop-all-the-variables-that-are-related-to-a-pre-speci
  remove_terms <- function(form, term) {
    fterms <- terms(form)
    fac <- attr(fterms, "factors")
    idx <- which(as.logical(fac[term, ]))
    new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
    return(formula(new_fterms))
  }  
  set.seed(1)
  res.list <- list()

  # extract all variables 
  all.v <- all.vars(formula(mod))
  all.v <- all.v[c(-1, -length(all.v))]
  all.v[length(all.v)+1] <- "whole"
  
  # get all terms in the model
  var <- attr(terms(mod), "term.labels")
  var[length(var)+1] <- "whole"
  
  # parametric bootstrap to get distribution
  for (k in 1:sims) {
    print(k)
    
    # simulate data and refir model
    ss <- simulate(mod)
    model <- refit(mod, newresp=ss)
    #mod.new[[1]] <- mod

    # make null model, ie no fixed effects    
    all <- paste(attr(terms(model), "term.labels"),sep="", collapse= " - ")
    null.mod <- update(model, as.formula(paste("~. -",all ))) # null model, no predictors
    
    # df for output
    out = data.frame(row.names = as.character(all.v), b.site = rep(NA,length(all.v)), w.site = rep(NA,length(all.v)), 
                     marginal = rep(NA,length(all.v)), conditional =  rep(NA,length(all.v)), stringsAsFactors=FALSE)
    
    for (i in 1:length(all.v)){
      if(only.whole==FALSE) {term = all.v[i]} else {
            term="whole"
            i=length(all.v)}
      
      if(term!="whole") {
        rem <- grep(all.v[i],var[-length(var)])
        rem <- var[-length(var)][rem]
        form <- as.formula(paste("~. -", paste(rem,sep="", collapse= " - "), sep=""))
        
        # # If species variable, keep main effect
        # if(var[i]=="Species") {form <- as.formula(paste("~. -", 
        #                                                paste(rem[-1],sep="", collapse= " - "), sep=""))}

        red.mod <- update(model,  form) # reduced model,w/o the variable
        
        r2s <- r2.fnc(model, red.mod)
        out[i,1] = r2s$r2bSite
        out[i,2] = r2s$r2wSite
        out[i,3] =  NA#r2s$r2MC[,"R2m"]
        out[i,4] =  NA#r2s$r2MC[,"R2c"]
        rm(r2s)
        
      }
      
      if(term=="whole") {   
        #i = length(var)
        r2s <- r2.fnc(model, null.mod)
        out[i,1] = r2s$r2bSite
        out[i,2] = r2s$r2wSite
        out[i,3] = r2s$r2MC[,"R2m"]
        out[i,4] = r2s$r2MC[,"R2c"]
        rm(r2s)
        
      }
      res.list[[k]] <- out
    }
  }
  return(res.list)
}



# Prediction and plot function ####
# function to plot, works with interactions as well
plot.model.fnc <- function (mod, back = TRUE, x.scale.orig = TRUE, orig.data = NULL, 
                        lam=1, shift = 0, mult = 1, ylab = NULL, vars.plot = NULL) {
  #data used in model
  mod.data <- mod@frame
  
  # function to backtransform power-transformation
  ptrans.back <- function (x,lam) {((exp(log(x)/lam))/mult) - shift} 
  
  if(is.null(vars.plot)) {vars.plot <- attr(terms(mod), "term.labels")}
  vars <- attr(terms(mod), "term.labels")
  
  # fix if there is an interaction and if so remove main term
  check.int <- unlist(strsplit(vars, "Species:"))
  if(any(duplicated(check.int)))  {#add <- check.int[anyDuplicated(check.int)]
    vars = vars[-which(grepl(":", vars))]}
  
  plots <- list()
  
  
  for (i in 3:length(vars)) {
    if(!(vars[i] %in% vars.plot)) {next}
    p.var = c("Species", vars[i])
    xa <- p.var[2] # variable on x-axis
    #if(grepl(":", xa)) {xa <- unlist(strsplit(xa,":"))[2]}
    ya <- colnames(mod.data)[1] # response variable

    p1 <- interactions::interact_plot(model=mod, pred = !! xa, modx = Species, plot.points = T,
                        partial.residuals = T, centered="none",
                        int.type = "confidence", interval = T)
    
    
    p.dat <- ggplot_build(p1)$data[[3]]
    l.dat <- ggplot_build(p1)$data[[2]]
    l.dat$Species <- p1$data$Species
    
    
                # # prepare prediction data
                # pred.dat <- setNames(data.frame(matrix(ncol = length(vars), nrow = 2000)), vars)
                # pred.dat[,which(colnames(pred.dat) %in% p.var)] <- expand.grid(
                #   levels(mod.data$Species), 
                #   seq(min(mod.data[,xa]), max(mod.data[,xa]), length=1000))
                # #pred.dat[,which(colnames(pred.dat) %in% p.var)] <- expand.grid(0:1, seq(-3,3, length=1000))
                # pred.dat[,which(!(colnames(pred.dat) %in% p.var))] <- 0
                # 
                # # make predictions
                # fe.tmp <- fixef(mod)
                # vcov.tmp <- as.matrix(vcov(mod))
                # n.sims <- 1000
                # sigmahat <- rep(1, n.sims)
                # 
                # # Make n.sims draws for each element of the fixed effects
                # betaSim <- abind::abind(lapply(1:n.sims,
                #                                function(x) mvtnorm::rmvnorm(n = 1, mean = fe.tmp, 
                #                                                             sigma = sigmahat[x]*vcov.tmp, method = "chol")), along=1)
                # # Calculate n.sims predictions for each row in PredDat
                # mm <- model.matrix(formula(paste('~ ',strsplit(as.character(formula(mod)[[3]]), "\\(")[[2]])),pred.dat)
                # #mm <- model.matrix(formula(paste('~ ', paste(colnames(pred.dat),collapse=" + "))),pred.dat)
                # mm[,which(attr(mm, "dimnames")[[2]] == "year")] <- 0.5 # take average over the two years
                # fixed <- mm %*% t(betaSim)
                # # For each row (observation) in PredDat calculate the median, upr and lwr 
                # Preds <- data.frame(fit = apply(fixed, 1, median), 
                #                     upr = apply(fixed, 1, quantile, 0.9), 
                #                     lwr = apply(fixed, 1, quantile, 0.1))
                # # backtransform...to be fixed
                # #Preds <- apply(Preds, 2, invlogit)
                # Preds <- cbind(pred.dat[,c("Species", p.var[2])], Preds)
    
    #get raw data values, for site-variables that means site mean and sd
    # Plot on the original predictor scale. Need to switch to correct column.
    if(x.scale.orig== TRUE) {
      orig.xa <- xa
      pred.orig <- c("prec_d", "evap_d", "par_d", "ev_pre",
                     "prec", "temp", "evap", "par", "hwt", "cover", "norain", 
                     "ndep", "Nmean", "NPmean")
      xa_num <- which(c("pr_d", "ev_d", "pa_d","ev_pre",
                        "pr", "tem", "ev", "pa", "wt", "cov", "nor", "ndeL", "Nm", "NPm") %in% xa)
      xa <- pred.orig[xa_num]
      mod.data[[xa]] <- orig.data[[xa]] # put original x-axis scale to the data frame
      l.dat[["x"]] <- (l.dat[["x"]] * sd(orig.data[[xa]])) + mean(orig.data[[xa]]) # backtransform scaled predictor to unscaled
      p.dat[["x"]] <- (p.dat[["x"]] * sd(orig.data[[xa]])) + mean(orig.data[[xa]]) # backtransform scaled predictor to unscaled
    }  

    if(!(xa %in% c("cover", "hwt"))) {
      p.dat <- cbind(p.dat, orig.data[,c("Site", "year", "Species")])
      p.dat$Site_year <- paste(p.dat$Site, p.dat$year, sep="_")
      #form = as.formula(paste("cbind(",ya,", ", xa, ") ~ Site + Species", sep=""))
      #agg.dat <- aggregate(form, mod.data, FUN = function (x) c(mean(x), sd(x)))
      form = as.formula(paste("cbind(","y",", ", "x", ") ~ Site_year + Species", sep=""))
      agg.dat <- aggregate(form, p.dat, FUN = function (x) c(mean(x), sd(x)))
      agg.dat <- data.frame(agg.dat[,1:2], agg.dat[,3], agg.dat[,4][,1])
      colnames(agg.dat)[3:5] <- c(ya, "se", xa)} else {
        p.dat$Species <- orig.data$Species
        agg.dat <-  p.dat[, c("y", "x", "Species")]
        agg.dat$se <- NA
        colnames(agg.dat)[c(1,2,4)] <- c(ya, xa, "se")
      }

    agg.dat$lo = agg.dat[,ya]-agg.dat[,"se"]; agg.dat$hi = agg.dat[,ya]+agg.dat[,"se"] # for within-site SD bars
    
    # backtransform responses if wanted
    if(back==TRUE) {
      agg.dat[, c(ya,"lo", "hi")] <- ptrans.back(agg.dat[, c(ya,"lo", "hi")], lam)
      l.dat[, c("y", "ymax", "ymin")] <- ptrans.back(l.dat[, c("y", "ymax", "ymin")], lam)
    }
    
    # plot
    xlab=NULL
    if(xa=="prec") {xlab = "Precipitation (mm)"} 
    if(xa=="temp") {xlab =  expression(paste("Mean temperature (",degree,"C)"))} 
    if(xa=="evap") {xlab = "Evaporation (mm)"} 
    if(xa=="cover") {xlab = "Vascular plant cover (%)"} 
    if(xa=="par") {xlab = expression(PAR ~ (W ~m^{-2}))} 
    if(xa=="ndeL") {xlab = expression(N ~deposition ~(g ~ m^{-2} ~yr^{-1}))} 
    if(xa=="Nmean") {xlab = "Tissue N concentration (%)"} 
    if(xa=="norain") {xlab = "Mean rain free period (days)"} 
    if(xa=="hwt") {xlab = "HWT (cm)"} 
    if(is.null(ylab)) {ylab = ya} # keep original response as y-axis label if nothing else is given
    if(is.null(xlab)) {xlab = xa} # keep original response as y-axis label if nothing else is given
    
    plots[[i-2]] <- ggplot(agg.dat, aes_string(x=xa, y=ya, ymin="lo", 
                                               ymax="hi", color = "Species", shape = "Species")) +
      geom_point(size=2)+ geom_linerange(show.legend = FALSE, alpha=0.2) +
      geom_ribbon(data=l.dat, aes_string(y="ymin", x="x", ymin = "ymin", ymax = "ymax", fill = "Species"), 
                  alpha = .15, show.legend = FALSE, linetype=0, inherit.aes = FALSE) +
      geom_line(data=l.dat, aes_string(y="y", x="x", color="Species"), 
                alpha =1, size = 1.2,inherit.aes = FALSE) +
      labs(y= ylab, 
           x= xlab) +
      ylim(c(min(agg.dat[,"lo"]), max(agg.dat[,"hi"]))) +
      xlim(c(min(agg.dat[,xa]), max(agg.dat[,xa]))) +
      scale_fill_manual(values=c("black", "red"), name="fill") +
      scale_color_manual(name = "Species", breaks = c("S.fuscum", "S.magellanicum"), 
                         values = c("S.fuscum" = "black", "S.magellanicum" = "red"),
                         labels = c("S. fuscum", "S. magellanicum")) +
      scale_shape_manual(name=c("Species"), breaks = c("S.fuscum", "S.magellanicum"),
                         values = c(1,2), labels = c("S. fuscum", "S. magellanicum")) +
      theme_bw()+
      theme(legend.justification = c(0, 0), 
            #legend.position = c(0.65, 0.03),
            legend.key.size = unit(2, "line"),
            legend.text = element_text(size=8, face = "italic"),
            legend.title = element_text(size=8),
            legend.position="none",
            axis.text = element_text(size=12),
            axis.title = element_text(size=14)) #+
    
  }
  # add legend plot
  dat.leg <- data.frame(species = unique(agg.dat$Species), y=rnorm(10,10), x=rnorm(10,10))
  plots[[i-1]] <- ggplot(dat.leg, aes(x=x, y=y, group=species)) +
    geom_point(aes(shape=species, color=species), size=4)+   
    geom_line(aes(color=species)) +
    ylim(c(-3.0, -0.5)) +
    xlim(c(-2.5, -0.5)) +
    scale_fill_manual(values=c("black", "red"), name="fill") +
    scale_color_manual(name = "Species", breaks = c("S.fuscum", "S.magellanicum"), 
                       values = c("S.fuscum" = "black", "S.magellanicum" = "red"),
                       labels = c("S. fuscum", "S. magellanicum")) +
    scale_shape_manual(name=c("Species"), breaks = c("S.fuscum", "S.magellanicum"),
                       values = c(1,2), labels = c("S. fuscum", "S. magellanicum")) +
    theme_bw()+
    theme(axis.ticks =  element_blank(),
          legend.position = c(0.5, 0.5),
          legend.box = "horizontal",
          axis.text  = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.text.align = 0,
          legend.text = element_text(size=16, face="italic"),
          legend.title = element_text(size=16))
  
  plots <-Filter(Negate(is.null), plots)
  rm(mod.data)
  return(plots) 
  #grid.arrange(plots[[1]], ncol=2)
  #grid.arrange(grobs=plots, ncol=2)
}




