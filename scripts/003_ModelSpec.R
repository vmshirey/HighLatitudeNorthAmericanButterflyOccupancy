run.R2jags.model <- function(d,
                             ni=1.01e5,
                             nb=1e3,
                             nt=1e2,
                             nc=3){

  res <- jagsUI::jags.basic(data=d$my.constants,
                            inits=list(d$my.inits, d$my.inits, d$my.inits),
                            parameters.to.save=d$my.params,
                            model.file="model.txt",
                            n.chains=nc,
                            n.thin=nt,
                            n.iter=ni,
                            n.burnin=nb,
                            parallel=TRUE)
  return(res)
}