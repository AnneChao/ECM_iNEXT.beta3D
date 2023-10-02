# ========================================================================================================== #

# Install and upload libraries ----

library(abind)
library(ggpubr)
library(plotly)
library(lmtest)
library(readxl)
library(ggplot2)
library(openxlsx)
library(magrittr)
library(parallel)
library(reshape2)
library(gridExtra)
library(tidyverse)
library(geosphere)
library(patchwork)
library(future.apply)


# ========================================================================================================== #
# 
# This code includes five parts:
# 
# (1) Figure 2. iNEXT.beta3D standardization for beta diversity between Edge and Interior habitats for two rainforest fragments in Brazil. 
# (2) Figure 3. Beta diversity between Edge and Interior habitats as a function of fragment size under six coverage values.
# (3) Figure 4. The iNEXT.beta3D standardization for two time periods for beetle data.
# (4) Figure 5. Trajectories of temporal beta diversity over time for tree species incidence data among 100 subplots in six second-growth rain forests.
# (5) Figure 6. Temporal and spatial gamma, alpha and beta diversity for fish species based on SWC-IBTS data from 1985 to 2010.
# 
# See "Brief guide" for details. 
# 
# ========================================================================================================== #


library(devtools)
install_github("AnneChao/iNEXT.3D")        # Press 'Enter' to skip update options
library(iNEXT.3D)
install_github("AnneChao/iNEXT.beta3D")    # Press 'Enter' to skip update options
library(iNEXT.beta3D)

source("Source R code.txt")



# ========================================================================================================== #
# Figure 2. iNEXT.beta3D standardization for beta diversity between Edge and Interior habitats for two rainforest fragments in Brazil.
raw_data = read_xlsx("Data rainforest trees.xlsx", sheet = 1)

data_for_beta = lapply(1:12, function(i) {
  tmp = raw_data[c(i, i+12), -(1:4)] %>% t
  colnames(tmp) = c('Edge', 'Interior')
  tmp
})
names(data_for_beta) = raw_data$Site[1:12]


## Figure 2 (a)
output_fig_2a = iNEXTbeta3D(data_for_beta[c('Marim', 'Rebio 2')], diversity = 'TD', nboot = 200, base = 'size')
fig_2a_or_4a(output_fig_2a)


## Figure 2 (b)
output_fig_2b = iNEXTbeta3D(data_for_beta[c('Marim', 'Rebio 2')], diversity = 'TD', nboot = 200, base = 'coverage')
fig_2b_or_4b(output_fig_2b)



# ========================================================================================================== #
# Figure 3. Beta diversity between Edge and Interior habitats as a function of fragment size under six coverage values.
output_fig_3 = list("iNEXTbeta" = iNEXTbeta3D(data_for_beta, "TD", q = c(0,1,2), nboot = 0, 
                                              level = c(0.6, 0.7, 0.8, 0.9, 0.95, 1)),
                    "obsbeta" = obsbeta(data_for_beta, q = c(0,1,2), nboot = 0))

fig_3(output_fig_3)



# ========================================================================================================== #
# Figure 4. The iNEXT.beta3D standardization for two time periods for beetle data.

## Figure 4 (a)
beetle = list('Logged'   = read.xlsx('Data beetles.xlsx', rowNames = T, sheet = 1),
              'Unlogged' = read.xlsx('Data beetles.xlsx', rowNames = T, sheet = 2))

output_fig_4a = iNEXTbeta3D(beetle, datatype = 'abundance', base = "size", nboot = 200)
fig_2a_or_4a(output_fig_4a)


## Figure 4 (b)
output_fig_4b = iNEXTbeta3D(beetle, datatype = 'abundance', base = 'coverage', nboot = 200)
fig_2b_or_4b(output_fig_4b)



# ========================================================================================================== #
# Figure 5. Trajectories of temporal beta diversity over time for tree species incidence data among 100 subplots in six second-growth rain forests.

Cuat.raw   = read.xlsx("Data second-growth trees.xlsx", sheet = 1)
LindEl.raw = read.xlsx("Data second-growth trees.xlsx", sheet = 2)
Tiri.raw   = read.xlsx("Data second-growth trees.xlsx", sheet = 3)
LindSu.raw = read.xlsx("Data second-growth trees.xlsx", sheet = 4)
FEB.raw    = read.xlsx("Data second-growth trees.xlsx", sheet = 5)
JE.raw     = read.xlsx("Data second-growth trees.xlsx", sheet = 6)


inci.raw = list(SGF.data.transf(Cuat.raw   %>% filter(Year >= 2005)), 
                SGF.data.transf(LindEl.raw %>% filter(Year >= 2005)), 
                SGF.data.transf(Tiri.raw   %>% filter(Year >= 2005)), 
                SGF.data.transf(LindSu.raw %>% filter(Year >= 2005)), 
                SGF.data.transf(FEB.raw), 
                SGF.data.transf(JE.raw))

age = data.frame(Assem = c("Cuatro Rios", "Lindero el Peje", "Tirimbina", "Lindero Sur", "Finca el Bejuco", "Juan Enriquez"),
                 Age = c(25,20,15,12,2,2))   ## fix on year 1997


cpu.cores <- detectCores() - 1
cl <- makeCluster(cpu.cores)
clusterExport(cl, varlist = c("inci.raw", "for_fig_5", "obsbeta"), envir = environment())
clusterEvalQ(cl, c(library(tidyverse), library(iNEXT.beta3D), library(iNEXT.3D), library(reshape2), 
                   library(magrittr), library(abind), library(future.apply)))

forests.output = parLapply(cl, inci.raw, function(x) for_fig_5(x, nboot = 200))

stopCluster(cl)

fig_5(forests.output)



# ========================================================================================================== #
# Figure 6. Temporal and spatial gamma, alpha and beta diversity for fish species based on SWC-IBTS data from 1985 to 2010.

fish <- read.csv("Data fish_Lat55-60.csv")
fish = fish %>% mutate(region = ifelse(LatBand %in% c(55.5, 56, 56.5, 57), 'South', 
                                       ifelse(!(LatBand %in% c(57.5, 60)), 'North', 'Other'))) %>% filter(region != 'Other')
groupyear = matrix(1985:2010, nrow = 2)
colnames(groupyear) = paste(groupyear[1,], groupyear[2,], sep = '~')


cpu.cores <- detectCores()-1
cl <- makeCluster(cpu.cores)
clusterExport(cl, varlist = c("rarefysamples", "fish", "groupyear", "obsbeta"), envir = environment())
clusterEvalQ(cl, c(library(tidyverse), library(iNEXT.beta3D), library(iNEXT.3D), library(reshape2)))

simu_output = parLapply(cl, 1:200, function(k) {
  region <- unique(fish$region)
  
  TSrf <- list()
  
  for(i in 1:length(region)){
    data2 <- fish[fish$region == region[i],]
    TSrf[[i]] <- rarefysamples(data2) 
  }
  names(TSrf) <- region
  
  rf <- do.call(rbind, TSrf)
  rf <- data.frame(rf, LatBand = rep(names(TSrf), times = unlist(lapply(TSrf, nrow))))
  rf <- rf[!is.na(rf$Year), ]
  rownames(rf)<-NULL
  data = rf
  
  cov = c(0.99, 0.999, 1)
  
  ## ================== Temporal ================== ##
  beta.temp = lapply(region, function(i) {
    tmp = data %>% filter(LatBand %in% i)
    
    tmp2 = lapply(2:length(unique(tmp$Year)), function(j) {
      g1 = dcast(tmp %>% filter(Year == sort(unique(tmp$Year))[1]), Species ~ LatBand, value.var = 'Abundance')
      
      g2 = dcast(tmp %>% filter(Year == sort(unique(tmp$Year))[j]), Species ~ LatBand, value.var = 'Abundance')
      
      out = full_join(g1, g2, by = 'Species')[,-1]
      
      out[is.na(out)] = 0
      out
    })
    names(tmp2) = sort(unique(tmp$Year))[-1]
    
    return(tmp2)
  })
  
  names(beta.temp) = region
  
  output.temp = lapply(1:length(beta.temp),  function(i) {
    
    result.obsbeta = obsbeta(beta.temp[[i]], q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
    
    result.iNEXTbeta = iNEXTbeta3D(beta.temp[[i]], q = c(0, 1, 2), datatype = 'abundance', level = cov, nboot = 0)
    
    rbind(lapply(result.iNEXTbeta, function(x) rbind(x$gamma %>% rename("Estimate" = "Gamma") %>% mutate(Type = "gamma"),
                                                     x$alpha %>% rename("Estimate" = "Alpha") %>% mutate(Type = "alpha"),
                                                     x$beta %>% rename("Estimate" = "Beta") %>% mutate(Type = "beta")) %>% 
                   select(c("Dataset", "Order.q", "SC", "Type", "Estimate"))) %>%
            do.call(rbind,.),
          
          result.obsbeta %>% do.call(rbind,.) %>% rename("Estimate" = "Diversity") %>% mutate(SC = 'Observed') %>% 
            select(c("Dataset", "Order.q", "SC", "Type", "Estimate"))) %>%
      
      mutate(Latitude = names(beta.temp)[i], .after = "Type")
    
  }) %>% do.call(rbind,.)
  
  output.temp$Dataset = as.numeric(output.temp$Dataset)
  
  ## ================== Spatial ================== ##
  beta.spat = lapply( list( c('South', 'North') ), function(i) {
    tmp1 = data %>% filter(LatBand %in% i[[1]])
    tmp2 = data %>% filter(LatBand %in% i[[2]])
    
    year = unique(data$Year) %>% sort
    
    tmp = lapply(year, function(j) {
      g1 = dcast(tmp1 %>% filter(Year == j), Species ~ LatBand, value.var = 'Abundance')
      g2 = dcast(tmp2 %>% filter(Year == j), Species ~ LatBand, value.var = 'Abundance')
      
      out = full_join(g1, g2, by = 'Species')[,-1]
      
      out[is.na(out)] = 0
      out
    })
    names(tmp) = year
    
    return(tmp)
  })
  
  names(beta.spat) = 'South vs. North'
  
  output.spat = lapply(1:length(beta.spat),  function(i) {
    
    result.obsbeta = obsbeta(beta.spat[[i]], q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
    
    result.iNEXTbeta = iNEXTbeta3D(beta.spat[[i]], q = c(0, 1, 2), datatype = 'abundance', level = cov, nboot = 0)
    
    rbind(lapply(result.iNEXTbeta, function(x) rbind(x$gamma %>% rename("Estimate" = "Gamma") %>% mutate(Type = "gamma"),
                                                     x$alpha %>% rename("Estimate" = "Alpha") %>% mutate(Type = "alpha"),
                                                     x$beta %>% rename("Estimate" = "Beta") %>% mutate(Type = "beta")) %>% 
                   select(c("Dataset", "Order.q", "SC", "Type", "Estimate"))) %>%
            do.call(rbind,.),
          
          result.obsbeta %>% do.call(rbind,.) %>% rename("Estimate" = "Diversity") %>% mutate(SC = 'Observed') %>% 
            select(c("Dataset", "Order.q", "SC", "Type", "Estimate"))) %>%
      
      mutate(Latitude = names(beta.spat)[i], .after = "Type")
    
  }) %>% do.call(rbind,.)
  
  output.spat$Dataset = as.numeric(output.spat$Dataset)
  
  list("temporal" = output.temp, "spatial" = output.spat)
})

stopCluster(cl)


output_fig_6 = list('temporal' = simu_output[[1]]$temporal,
                    'spatial'  = simu_output[[1]]$spatial)

for (i in 2:length(simu_output)) {
  output_fig_6$temporal = full_join(output_fig_6$temporal, 
                                    simu_output[[i]]$temporal, 
                                    by = c('Order.q', 'SC', 'Dataset', 'Type', 'Latitude'))

  output_fig_6$spatial = full_join(output_fig_6$spatial, 
                                   simu_output[[i]]$spatial, 
                                   by = c('Order.q', 'SC', 'Dataset', 'Type', 'Latitude'))
}

output_fig_6$temporal = cbind(output_fig_6$temporal[,1:5], 'Estimate' = apply(output_fig_6$temporal[,-(1:5)], 1, mean))
output_fig_6$spatial  = cbind(output_fig_6$spatial[,1:5],  'Estimate' = apply(output_fig_6$spatial[,-(1:5)],  1, mean))

fig_6a(output_fig_6)
fig_6b(output_fig_6)


