run.R2jags.model <- function(dd,
                             ni=1.1e4,
                             nb=1e3,
                             nt=25,
                             nc=1,
                             precip=FALSE){
  
  
  
  if(precip==FALSE){
    model.file = "model_intercept_temp.txt"
    my.params <- c("mu.p.0", "p.yr", "sigma.p.sp", "sigma.p.site", "mu.psi.0",
                   "psi.sp",  "sigma.psi.sp", "psi.site", "sigma.psi.site", 
                   "sigma.psi.sp.temp", "sigma.psi.sp.temp2", "sigma.psi.sp.precip",
                   "psi.beta.temp", "psi.beta.temp2", "psi.beta.precip", "psi.area")
  } else{
    model.file = "model_intercept_precip.txt"
    my.params <- c("mu.p.0", "p.yr", "sigma.p.sp", "sigma.p.site", "mu.psi.0",
                   "psi.sp",  "sigma.psi.sp", "psi.site", "sigma.psi.site", 
                   "sigma.psi.sp.temp", "sigma.psi.sp.temp2", "sigma.psi.sp.precip",
                   "psi.beta.temp", "psi.beta.temp2", "psi.beta.precip", "psi.area")
  }
  
  # Run the model
  my.res <- jagsUI::jags.basic(model.file=model.file,
                                data=c(dd$my.data, 
                                       dd$my.constants),
                                inits=list(dd$my.inits,
                                           dd$my.inits, 
                                           dd$my.inits,
                                           dd$my.inits),
                                n.chains=nc,
                                n.adapt=1e3,
                                n.burnin=nb,
                                n.thin=nt,
                                n.iter=ni,
                               parallel=TRUE,
                               parameters.to.save=my.params)

  return(my.res)
}
