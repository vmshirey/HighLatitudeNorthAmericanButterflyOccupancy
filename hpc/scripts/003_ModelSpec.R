run.R2jags.model <- function(d,
                             ni=5.1e4,
                             nb=1e3,
                             nt=100,
                             nc=1){

  res <- jagsUI::jags.basic(data=d$my.constants,
			    inits=list(d$my.inits),
                            parameters.to.save=d$my.params,
                            model.file="scripts/model.txt",
                            n.chains=nc,
                            n.thin=nt,
                            n.iter=ni,
                            n.burnin=nb,
                            parallel=FALSE)
  return(res)
}
