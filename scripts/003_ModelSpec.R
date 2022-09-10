run.R2jags.model <- function(d,
                             ni=1.1e4,
                             nb=1e3,
                             nt=25,
                             nc=1,
                             phy=TRUE){
  
  
  
  if(phy == TRUE){
    model.file = "model_phy.txt"
    my.params <- c("mu.p.0", "p.yr", "sigma.p.sp", "sigma.p.site",
                   "sp.psi.0", "sigma.psi.sp", "psi.yr", "sigma.psi.yr", "psi.site",  "sigma.psi.site", 
                   "psi.temp.env", "psi.temp.trait", "psi.precip.env", "psi.precip.trait",
                   "lambda.temp", "lambda.precip", "sigma.psi.temp", "sigma.psi.precip",
                   "psi.beta.temp", "psi.beta.precip", "psi.area")
  }
  if(phy == FALSE){
    model.file = "model_intercept.txt"
    my.params <- c("mu.p.0", "p.yr", "sigma.p.sp", "sigma.p.site", "mu.psi.0",
                   "sigma.psi.sp", "psi.yr", "sigma.psi.yr", "psi.site",  "sigma.psi.site", 
                   "sigma.psi.sp.temp", "sigma.psi.sp.precip",
                   "sigma.psi.sp.temp", "sigma.psi.sp.precip",
                   "psi.beta.temp", "psi.beta.precip", "psi.area")
  }
  
  res <- jagsUI::jags.basic(data=c(d$my.data, d$my.constants),
                            inits=list(d$my.inits, d$my.inits, d$my.inits),
                            parameters.to.save=my.params,
                            model.file=model.file,
                            n.chains=nc,
                            n.thin=nt,
                            n.iter=ni,
                            n.burnin=nb, parallel=TRUE)
  
  return(res)
}
