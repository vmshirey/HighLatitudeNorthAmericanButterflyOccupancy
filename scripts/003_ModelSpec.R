run.R2jags.model <- function(d,
                             ni=5.1e4,
                             nb=1e3,
                             nt=1e2,
                             nc=1){
  
  model.jags <- function(){
    
    ## DETECTION COMPONENTS ##
    # Detection intercept
    mu.p.0 ~ dnorm(0,10) # MONITOR
    
    # Effect of occupancy interval on deteection
    p.yr ~ dnorm(0,10) # MONITOR
    
    # Random effect of species on detection
    sigma.p.sp ~ dunif(0,10) # MONITOR
    for(sp in 1:nsp){
      p.sp[sp] ~ dnorm(0,sigma.p.sp) 
    }
    
    # Random effect of site/year on detection
    sigma.p.site ~ dunif(0,20) # MONITOR
    for(site in 1:nsite){
      for(yr in 1:nyr){
        p.site[site,yr] ~ dnorm(0, sigma.p.site) 
      }
    }
    
    ## OCCUPANCY COMPONENTS ##
    # Species-specific occupancy intercepts
    sigma.psi.sp ~ dunif(0,10) # MONITOR
    for(sp in 1:nsp){
      sp.psi.0[sp] ~ dnorm(0,sigma.psi.sp) # MONITOR
    }
    
    # Effect of terrestrial surface area on occupancu
    p.area ~ dnorm(0,10) # MONITOR
    
    # Random effect of site on occupancy
    sigma.psi.site ~ dunif(0,10) # MONITOR
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, sigma.psi.site) # MONITOR
    }
    
    # Random effect of occupancy interval on occupancy
    sigma.psi.yr ~ dunif(0,10) # MONITOR
    for(yr in 1:nyr) {
      psi.yr[yr] ~ dnorm(0, sigma.psi.yr) # MONITOR
    }
    
    # Create species-specific slopes with respect to environmental covariates
    # Phylogeny gets incorporated here:
    psi.temp.env ~ dnorm(0,0.001)
    psi.temp.trait ~ dnorm(0,0.001)
    psi.precip.env ~ dnorm(0,0.001)
    psi.precip.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp){
      mu.beta.temp[sp] <- psi.temp.env + psi.temp.trait*temp.trait[sp]
      mu.beta.precip[sp] <- psi.precip.env + psi.precip.trait*precip.trait[sp]
    }
    
    # Incoporate the phylogenetic signal
    lambda.temp ~ dunif(0,1)
    lambda.precip ~ dunif(0,1)
    
    beta.mat.temp[1:nsp,1:nsp] <- lambda.temp*VCOV[,] + (1-lambda.temp)*ID[,]
    beta.mat.precip[1:nsp,1:nsp] <- lambda.precip*VCOV[,] + (1-lambda.precip)*ID[,]
    
    sigma.psi.temp ~ dunif(0,10)
    sigma.psi.precip ~ dunif(0,10)
    tau.beta.sp.mat.temp[1:nsp,1:nsp] <- inverse((sigma.psi.temp^2)*beta.mat.temp[,])
    psi.beta.temp[1:nsp] ~ dmnorm(mu.beta.temp[], tau.beta.sp.mat.temp[,])
    tau.beta.sp.mat.precip[1:nsp,1:nsp] <- inverse((sigma.psi.precip^2)*beta.mat.precip[,])
    psi.beta.precip[1:nsp] ~ dmnorm(mu.beta.precip[], tau.beta.sp.mat.precip[,])
    
    ## CORE LIKELIHOOD FUNCTIONS ##
    for(site in 1:nsite){
      for(yr in 1:nyr){
        for(sp in 1:nsp){
          
          logit(psi[site,yr,sp]) <-
            sp.psi.0[sp]+
            p.area*gridarea[site]+
            psi.beta.temp[sp]*temp[site,yr]+
            psi.beta.precip[sp]*precip[site,yr]+
            psi.site[site]+
            psi.yr[yr]
          
          logit(p[site,yr,sp]) <-
            mu.p.0+
            p.yr*yr+
            p.sp[sp]+
            p.site[site,yr]
          
        }
      }
    }
    
    # Latent state
    for(site in 1:nsite){
      for(yr in 1:nyr){
        for(sp in 1:nsp){
          
          # Draw latent states
          Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
          
        }
      }
    }
    
    # Likelihood
    for(ind in 1:nind) {
      
      mu.p[ind] <-
        Z[sitev[ind],yrv[ind],spv[ind]]*
        p[sitev[ind],yrv[ind],spv[ind]]
      
      X[ind] ~ dbern(mu.p[ind])
      
    }
  }
  
  my.params <- c("mu.p.0", "p.yr", "sigma.p.sp", "sigma.p.site", # DETECTION COMPONENTS
                 "sp.psi.0", "sigma.psi.sp", "psi.yr", "sigma.psi.yr", "psi.site",  "sigma.psi.site", 
                 "psi.temp.env", "psi.temp.trait", "psi.precip.env", "psi.precip.trait",
                 "lambda.temp", "lambda.precip", "sigma.psi.temp", "sigma.psi.precip",
                 "psi.beta.temp", "psi.beta.precip", "p.area")
  
  res <- R2jags::jags(data=d$my.constants,
                      inits=list(d$my.inits, d$my.inits, d$my.inits),
                      parameters.to.save=my.params,
                      model.file=model.jags,
                      n.chains=nc,
                      n.thin=nt,
                      n.iter=ni,
                      n.burnin=nb,
                      working.directory=NULL)
  
  return(res)
}