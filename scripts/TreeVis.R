bounded.hpp <- function(sample,
                        lower=0,
                        upper=1,
                        mode=TRUE, 
                        HPDcoverage=0.95,
                        codaHPD=TRUE,
                        n=1024) {
  
  reflectedSample <- c(if(lower > -Inf) lower - (sample-lower) else numeric(),
                       sample,
                       if(upper < Inf) upper + (upper-sample) else numeric())
  
  griddedDensityCall <- quote(density(reflectedSample, n=n, bw="SJ"))
  if(lower > -Inf) griddedDensityCall[['from']] <- lower
  if(upper < Inf) griddedDensityCall[['to']] <- upper
  griddedDensity <- eval(griddedDensityCall)
  
  iMax <- which.max(griddedDensity$y)
  ## below is the value of the mode, to within precision of n
  ## (default=1024 points)
  argMaxX <- griddedDensity$x[iMax] 
  if(is.numeric(HPDcoverage)) {
    estimatedPDF <- (griddedDensity$y / sum(griddedDensity$y))
    
    alpha <- 1-HPDcoverage
    
    funThatWeWantToSolveForZero <- function(densityLevel) {
      whichAbove <- which(estimatedPDF >= densityLevel)
      if(length(whichAbove)==0) return(1-alpha)
      if(length(whichAbove)==n) return(-alpha)
      interval <- griddedDensity$x[range(whichAbove)]
      mean(sample < interval[1] | sample > interval[2])  - alpha
    }
    
    HPDsolution <- uniroot(funThatWeWantToSolveForZero,
                           lower=0,
                           upper=estimatedPDF[iMax]) 
    HPDdensityLevel <- HPDsolution$root
    whichAbove <- which(estimatedPDF > HPDdensityLevel)
    HPDinterval <- griddedDensity$x[range(whichAbove)]
  } else {
    HPDinterval <- NULL
  }
  
  list(mode=argMaxX, HPDinterval=HPDinterval)
  
}

# Visualize the tree data, first read in the tree data and species' traits
my_tree <- readRDS("../output/tree_topology.rds")
my_traits <- read.csv("../data/taxa/species_range_clim.csv")
my_traits <- data.frame(binomial=my_tree$tip.label) %>% left_join(my_traits)
my_tree <- as_tibble(my_tree)
my_tree <- left_join(my_tree, my_traits, by=c("label"="binomial"))
my_tree <- as.treedata(my_tree)

# Read in the MCMC sample data
dat <- read_rds("../output/samples/res_200_1.rds")
lambdas <- data.frame(temp=dat$BUGSoutput$sims.list$lambda.temp,
                      precip=dat$BUGSoutput$sims.list$lambda.precip)
lambda.temp.hpp <- bounded.hpp(lambdas$temp)
lambda.precip.hpp <- bounded.hpp(lambdas$precip)

# Read in images for plotting
bol_frig <- readPNG("../output/main/img/bfrigga.png")
bol_frig <- rasterGrob(bol_frig, interpolate=TRUE)

col_nast <- readPNG("../output/main/img/cnastes.png")
col_nast <- rasterGrob(col_nast, interpolate=TRUE)

ere_epip <- readPNG("../output/main/img/eepip.png")
ere_epip <- rasterGrob(ere_epip, interpolate=TRUE)

ica_saep <- readPNG("../output/main/img/isaep.png")
ica_saep <- rasterGrob(ica_saep, interpolate=TRUE)

pap_cani <- readPNG("../output/main/img/pcanadensis.png")
pap_cani <- rasterGrob(pap_cani, interpolate=TRUE)

# Create a plot (x=150, y=110 good point for photo annotations)
traits_tree <- ggtree::ggtree(my_tree)+
  geom_hilight(node=120, fill="#404040", alpha=0.2)+ # Hesperiidae
  geom_hilight(node=204, fill="#404040", alpha=0.2)+ # Lycaenidae
  geom_hilight(node=136, fill="#909090", alpha=0.2)+ # Pieridae
  geom_hilight(node=151, fill="#909090", alpha=0.2)+ # Nymphalidae
  geom_hilight(node=226, fill="#909090", alpha=0.2)+ # Papilionidae
  geom_text(mapping=aes(x=40, y=112), label="Nymphalidae", color="#cacbcc")+
  geom_text(mapping=aes(x=52, y=60), label="Lycaenidae", color="#cacbcc")+
  geom_text(mapping=aes(x=32, y=37), label="Pieridae", color="#cacbcc")+
  geom_text(mapping=aes(x=37, y=22), label="Hesperiidae", color="#cacbcc")+
  geom_text(mapping=aes(x=34, y=6), label="Papilionidae", color="#cacbcc")+
  geom_tree2(color="#cacbcc", layout="rectangular")+
  geom_tippoint(mapping=aes(fill=ave.temp2), shape=21, color="#cacbcc")+
  scale_fill_continuous_divergingx(palette="RdYlBu",
                        name="Average Annual Range-wide \nTemperature (°C)",
                        limits=c(-12,24), breaks=c(-12, 0, 12, 24), rev=TRUE,
                        mid=4.62, l3=21, p3=0.05, p4=0.5, c3=40, h3=258)+
  ggnewscale::new_scale_fill()+
  geom_tippoint(mapping=aes(fill=ave.precip2), shape=22, color="#cacbcc",
                position=position_nudge(x=0.85))+
  scale_fill_continuous_sequential(palette="YlGnBu",
                        name="Average Annual Range-wide \nPrecipitation (mm)",
                        limits=c(100,1500), rev=TRUE,
                        breaks=c(100, 300, 500, 700, 900, 1100, 1300, 1500))+
  geom_tiplab(size=2.5, 
              color="#cacbcc", nudge_x=1, fontface="italic")+
  geom_treescale(fontsize=4, x=0, y=50, color="#cacbcc", width=10)+
  geom_text(mapping=aes(x=5, y=55), label="Branch-length \n(Millions of Years)", color="#cacbcc")+
  annotation_custom(bol_frig, xmin=140, xmax=170, ymin=100, ymax=119)+
  geom_text(mapping=aes(x=154, y=99), label="Boloria frigga", 
            color="#cacbcc", fontface="italic", size=3)+
  annotation_custom(col_nast, xmin=140, xmax=170, ymin=25, ymax=45)+
  geom_text(mapping=aes(x=155, y=26), label="Colias nastes", 
            color="#cacbcc", fontface="italic", size=3)+
  annotation_custom(ere_epip, xmin=140, xmax=170, ymin=75, ymax=95)+
  geom_text(mapping=aes(x=155, y=74), label="Erebia epipsodea",
            color="#cacbcc", fontface="italic", size=3)+
  annotation_custom(ica_saep, xmin=140, xmax=170, ymin=50, ymax=70)+
  geom_text(mapping=aes(x=155, y=49), label="Icaricia saepiolus", 
            color="#cacbcc", fontface="italic", size=3)+
  annotation_custom(pap_cani, xmin=140, xmax=170, ymin=1, ymax=20)+
  geom_text(mapping=aes(x=155, y=0), label="Pterourus canadensis", 
            color="#cacbcc", fontface="italic", size=3)+
  ggplot2::xlim(0,170)+
  theme_cowplot()+
  theme(legend.position="bottom",
        legend.key.width=unit(1.5, "cm"),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_rect(fill="#1e2126"),
        plot.background=element_rect(fill="#1e2126"),
        text = element_text(color="#cacbcc"))

temp_trait_pagel <- ggplot()+
  geom_density(lambdas, mapping=aes(x=temp), fill="#cacbcc", color="#cacbcc", alpha=0.5)+
  geom_pointrange(data.frame(), mapping=aes(x=lambda.temp.hpp$mode, y=1, 
                                  xmin=lambda.temp.hpp$HPDinterval[1], 
                                  xmax=lambda.temp.hpp$HPDinterval[2]),
                  color="#cacbcc")+
  geom_text(mapping=aes(x=0.05, y=4), label="(a)",
            color="#cacbcc", size=4.5)+
  geom_vline(xintercept=0.5, linetype=2, color="#cacbcc", alpha=0.5)+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,4))+
  labs(y="Posterior Density", x=expression("Pagel's"~lambda[' temp.']))+
  theme_cowplot()+
  theme(legend.position="none",
        axis.text = element_text(color="#cacbcc"),
        axis.line = element_line(color="#cacbcc"),
        axis.ticks = element_line(color="#cacbcc"),
        panel.background=element_rect(fill="#1e2126"),
        plot.background=element_rect(fill="#1e2126"),
        text = element_text(color="#cacbcc"))

precip_trait_pagel <- ggplot()+
  geom_density(lambdas, mapping=aes(x=precip), fill="#cacbcc", color="#cacbcc", alpha=0.5)+
  geom_pointrange(data.frame(), mapping=aes(x=lambda.precip.hpp$mode, y=1, 
                                            xmin=lambda.precip.hpp$HPDinterval[1], 
                                            xmax=lambda.precip.hpp$HPDinterval[2]),
                  color="#cacbcc")+
  geom_text(mapping=aes(x=0.05, y=4), label="(b)",
            color="#cacbcc", size=4.5)+
  geom_vline(xintercept=0.5, linetype=2, color="#cacbcc", alpha=0.5)+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,4))+
  labs(y="Posterior Density", x=expression("Pagel's"~lambda[' precip.']))+
  theme_cowplot()+
  theme(legend.position="none",
        axis.text = element_text(color="#cacbcc"),
        axis.line = element_line(color="#cacbcc"),
        axis.ticks = element_line(color="#cacbcc"),
        panel.background=element_rect(fill="transparent"),
        plot.background=element_rect(fill="transparent"),
        text = element_text(color="#cacbcc"))

full_figure_two <- ggdraw()+
  draw_plot(traits_tree)+
  draw_plot(temp_trait_pagel, x=0, y=0.78, width=0.17, height=0.2)+
  draw_plot(precip_trait_pagel, x=0, y=0.58, width=0.17, height=0.2)

full_figure_two
ggsave2("../output/main/Figure2.png", full_figure_two, dpi=350, height=10, width=13)


