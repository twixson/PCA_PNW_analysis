#######
# PCA exploration
#######
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#
#     THIS SCRIPT HAS A FEW COMPUTATIONALLY EXPENSIVE COMPONENTS
#         (it takes 3-4 hours to run on my old laptop)
#
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!

library(tidyverse)
library(evd)
library(sf)
library(tictoc)
library(cowplot)
library(Matrix)
library(maps)
setwd("PATH_TO_DATA")
precip_data <- readRDS("../extracted_yearly/pnw_precip_data.rds")
  # This contains the pre-processed precipitation data after sub-setting to 
  #   the analysis region but before marginal transformation.  

#####
# Transformed linear operations - for TPDM estimation and reconstruction
#####
f <- function(y){
  x <- ifelse(y > 100, y, log(1+exp(y)))
  return(x)
}

finv <- function(x){
  y <- ifelse(x > 100, x, log(exp(x)-1))
}

tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}

tmult <- function(a,b){
  return(f(a*finv(b)))
}

###
# calculate tpdf
TPDF <- function(vec1, vec2, thresh = 0.975){
  rt <- vec1^2 + vec2^2 # don't sqrt() since you would square it later anyways
  r0 <- quantile(rt, probs = thresh)
  inclusion_indices <- which(rt > r0)
  rt <- rt[inclusion_indices]
  vec1 <- vec1[inclusion_indices]
  vec2 <- vec2[inclusion_indices]
  
  return(mean(vec1 * vec2 / rt))
}



#####
# Check the tail index with Likelihood Ratio Tests
#####
set.seed(9346732)
large_data <- matrix(NA, nrow = 200, ncol = dim(precip_data)[2] - 1)
for(i in 2:dim(precip_data)[2]){
  large_data[,(i-1)] <- sort(precip_data[-c(1,2), i], decreasing = T)[1:200]
}

lat_vec <- precip_data[2, -1]
lon_vec <- precip_data[1, -1]
date_vec <- precip_data[-c(1,2), 1]

get_scales <- function(data, temp_shape) {
  scales <- numeric(dim(data)[2])
  temp_LL <- 0
  for(i in 1:(dim(data)[2])){
    temp_out <- fpot(data[, i], threshold = min(data[,i]), 
                     shape = temp_shape, 
                     std.err = F)
    temp_LL <- temp_LL - unname(logLik(temp_out))
    scales[i] <- temp_out$estimate[1]
    if(i%%5000 == 0){
      print(paste0("iteration: ", i))
    }
  }
  return(list("scales" = scales, "temp_LL" = temp_LL))
}

get_scales_and_shapes <- function(data) {
  scales <- shapes <- shapes_se <- numeric(dim(data)[2])
  for(i in 1:(dim(data)[2])){
    temp_out <- fpot(x = data[, i], threshold = min(data[,i]) - 0.0001)
    scales[i] <- temp_out$estimate[1]
    shapes[i] <- temp_out$estimate[2]
    shapes_se[i] <- temp_out$std.err[2]
    if(i%%5000 == 0){
      print(paste0("iteration: ", i))
    }
  }
  return(list(scales = scales, shapes = shapes, shapes_se = shapes_se))
}

obj_pot <- function(guess, data, scales.){
  if(guess <= 0){
    return(99999999)
  } else {
    log_vals <- rep(NA, dim(data)[2])
    for(i in 1:dim(data)[2]){
      log_vals[i] <- pot_log_likelihood(guess, 
                                        data_vec = data[,i], 
                                        scale. = scales.[i])
    }
    return(sum(log_vals))
  }
}

pot_log_likelihood <- function(guess, data_vec, scale.){
  -sum(dgpd(data_vec[1:200], 
           scale = scale., 
           shape = guess, 
           loc = data_vec[200]-0.0001, 
           log = T))
}

fit_list <- list()
out1 <- get_scales_and_shapes(large_data)
fit_list[["scale_it_0"]] <- out1$scales
fit_list[["shape_it_0"]] <- optimize(obj_pot, lower = 0, upper = 999,
                                     data = large_data, 
                                     scales. = fit_list[["scale_it_0"]])

for(iteration in 1:2){
  shape_it_old <- paste0("shape_it_", iteration-1)
  scale_it_old <- paste0("scale_it_", iteration-1)
  scale_it <- paste0("scale_it_", iteration)
  shape_it <- paste0("shape_it_", iteration)
  
  fit_list[[scale_it]] <- get_scales(large_data, 
                                     temp_shape = fit_list[[shape_it_old]]$minimum)
  fit_list[[shape_it]] <- optimize(obj_pot, lower = 0, upper = 999,
                                   data = large_data, 
                                   scales. = fit_list[[scale_it]]$scales)
}

while(abs(fit_list[[scale_it]]$temp_LL[1] - 
          fit_list[[scale_it_old]]$temp_LL[1]) > 1){
  iteration <- iteration + 1
  shape_it_old <- paste0("shape_it_", iteration-1)
  scale_it_old <- paste0("scale_it_", iteration-1)
  scale_it <- paste0("scale_it_", iteration)
  shape_it <- paste0("shape_it_", iteration)
  
  fit_list[[scale_it]] <- get_scales(large_data, 
                                     temp_shape = fit_list[[shape_it_old]]$minimum)
  fit_list[[shape_it]] <- optimize(obj_pot, lower = 0, upper = 999,
                                   data = large_data, 
                                   scales. = fit_list[[scale_it]]$scales)
}

saveRDS(list(scales = fit_list$scale_it_6$scales, 
             alpha = fit_list$shape_it_6$minimum), 
        "./fitted_scales_and_alpha.rds")

# Scale plot
scale_df <- data.frame("lat" = lat_vec,
                       "lon" = lon_vec,
                       "scales" = fit_list$scale_it_6$scales)
my_sf <- st_as_sf(scale_df, coords = c('lon', 'lat'), crs = 4326)
western_us <- st_as_sf(maps::map("state", fill=T, plot=F,
                                 regions = c("washington", "oregon", "california")),
                       crs=4326)
#Plot it:
ggplot(my_sf) +
  geom_sf(aes(color = scales), shape = 15, size = 1.2) +
  geom_sf(data = western_us, fill = NA, color = "white") +
  xlim(-125, -122) +
  ylim(42, 49) +
  theme(axis.text.x = element_text(angle = 30))

#####
# liklihood ratio test and CI inclusion
#####
restricted_shape <- unname(fit_list[[shape_it]]$minimum)
shape_fixed_LL <- 0 

for(i in 1:(dim(large_data)[2])){
  temp_out <- fpot(large_data[1:200,i], 
                   threshold = large_data[200,i] - 0.0001, 
                   shape = restricted_shape)
  shape_fixed_LL <- shape_fixed_LL - logLik(temp_out)
  if(i%%1000 == 0){
    print(paste0("iteration: ", i))
  }
}

unrestricted_LL <- 0 
inside_CI <- zero_CI <- rep(NA, dim(large_data)[2])
mle_alpha <- rep(NA, dim(precip_data)[2]-1)
se_alpha <- rep(NA, dim(precip_data)[2]-1)

for(i in 1:(dim(large_data)[2])){
  temp_out <- fpot(large_data[1:200,i], 
                   threshold = large_data[200,i] - 0.0001)
  unrestricted_LL <- unrestricted_LL - logLik(temp_out)
  temp_ci <- confint(temp_out, 2)
  if(min(temp_ci) < restricted_shape){
    if(max(temp_ci) > restricted_shape){
      inside_CI[i] <- "inside"
    } else {
      inside_CI[i] <- "below"
    } 
  } else { 
    inside_CI[i] <- "above"
  }
  
  if(min(temp_ci) < 0){
    if(max(temp_ci) > 0){
      zero_CI[i] <- "covers"
    } else {
      zero_CI[i] <- "below"
    } 
  } else { 
    zero_CI[i] <- "above"
  }
  
  mle_alpha[i-1] <- temp_out$estimate[2]
  se_alpha[i-1] <- temp_out$std.err[2]
  if(i%%1000 == 0){
    print(paste0("iteration: ", i))
  }
}

lrt <- -2*(unrestricted_LL[1] - shape_fixed_LL[1])
pchisq(lrt, df = length(inside_CI), lower.tail = F)
# Reject the null here, the power is astronomical due to large sample size

###
# CI inclusion?


CI_inclusion <- data.frame("lat" = lat_vec, 
                           "lon" = lon_vec, 
                           "inside" = as.factor(inside_CI))

my_sf <- st_as_sf(CI_inclusion, coords = c('lon', 'lat'), crs = 4326)
western_us <- st_as_sf(maps::map("state", fill=T, plot=F, 
                                 regions = c("washington", "oregon", "california")), 
                       crs=4326)

#Plot it:
pdf("./shape_param_CI_inclusion_plot.pdf", width = 6, height = 9)

ggplot(my_sf) + 
  geom_sf(aes(color = inside), shape = 15, size = 1.2) + 
  geom_sf(data = western_us, fill = NA, color = "white") +
  scale_color_manual(values = c("orange", "black", "light blue")) + 
  xlim(-125, -122) + 
  ylim(42, 49) + 
  theme(axis.text.x = element_text(angle = 30))
dev.off()



###
# Map of MLE's for \alpha
#precip_data <- readRDS("./all_precip_data.rds")
mle_alpha <- rep(NA, dim(precip_data)[2]-1)
se_alpha <- rep(NA, dim(precip_data)[2]-1)
separate_LL <- 0 

upper <- mle_alpha + 2 * se_alpha
lower <- mle_alpha - 2 * se_alpha

mle_df <- data.frame("index" = 1:length(mle_alpha), 
                     "mle" = mle_alpha, 
                     "se" = se_alpha, 
                     "upper" = upper, 
                     "lower" = lower, 
                     "lat" = precip_data[2, -1], 
                     "lon" = precip_data[1, -1])

ggplot(mle_df[1:5000, ], aes(x = index, y = mle, ymin = lower, ymax = upper)) + 
  geom_linerange() + 
  geom_abline(aes(intercept = 0.05, slope = 0, col = "red"))

my_sf <- st_as_sf(mle_df, coords = c('lon', 'lat'))

pdf("./mle_quantile_map.pdf", width = 3, height = 9)
ggplot(my_sf) + 
  geom_sf(aes(color = mle))
dev.off()




#######
# Location-wise marginal rank transformation to Frechet with tail index 2
#######
transform_marginal <- function(data_vec){
  zeros_ind <- which(data_vec == 0)
  data_vec <- rank(data_vec) / (length(data_vec) + 1)
  data_vec[zeros_ind] <- 1 / (length(data_vec) + 1)
  return(qfrechet(data_vec, shape = 2))
}

tic()
precip_alpha2 <- apply(precip_data[-c(1,2), -c(1)], 2, transform_marginal)
toc()

saveRDS(precip_alpha2, "./precip_data_frechet_alpha2.rds")

# lat and lon df 
lat_lon <- data.frame("lon" = precip_data[1, -c(1)],
                      "lat" = precip_data[2, -c(1)])
dates <- precip_data[-c(1,2), 1]

saveRDS(list("lat_lon" = lat_lon, "date_vec" = dates), 
        "./lat_lon_dates_list.rds")
# rm(precip_data)



##########
# Compute the tpdf and perform eigen decomposition (!!!expensive!!!)
#########
# precip_alpha2 <- readRDS("./precip_data_frechet_alpha2.rds")
# setup tpdm 
tpdm <- matrix(NA, nrow = dim(precip_alpha2)[2], ncol = dim(precip_alpha2)[2])
tic("full tpdm estimation")
tic("finished the 50 th row")
for(i in 1:(dim(tpdm)[1])){
  for(j in i:dim(tpdm)[1]){
    temp_tpdm <- TPDF(precip_alpha2[, i], precip_alpha2[, j])
    tpdm[i,j] <- temp_tpdm
    tpdm[j,i] <- temp_tpdm
  }
  if(i%%50 == 0){
    toc()
    tic(paste0("finished the ", i+50, "th row"))
  } 
}
toc()
toc()
tic("eigen decomp")
#eigen decomp
a <- eigen(tpdm, symmetric = T) 
toc()
tic("near pd comp")
# save.image("./FILENAME.Rdata")
PD_tpdm <- nearPD(tpdm)
toc()
tic("second eigen decomp")
eigen_tpdm <- eigen(PD_tpdm$mat)
toc()





#####
# Eigen reconstruction
##### 

# functions 
plot_eigen <- function(col., data, eigen_vals., lat_lon. = lat_lon){
  temp_df <- data.frame("lat" = lat_lon.$lat, 
                        "lon" = lat_lon.$lon, 
                        "eigen" = data[,col.])
  my_sf <- st_as_sf(temp_df, coords = c("lon", "lat"), crs = 4326)
  
  plot1 <- ggplot(my_sf) + 
    geom_sf(aes(col = eigen), shape = 15, size = 1.5) + 
    scale_color_viridis_c() + # put "limits = c(-0.037, 0.037)" in to get on same scale
    labs(title = paste0("Map of Eigenvector ", col.), 
         subtitle = paste0("Eigenvalue = ", round(eigen_vals.[col.], 2))) + 
    theme(axis.text = element_text(angle = 45))
  
  return(plot1)
}

reconstruction_coefs <- function(data_vec, eigenvecs){
  data_vec[which(data_vec == 0)] <- 0.00001
  t(eigenvecs) %*% finv(data_vec)
}

reconstruct_storm <- function(coefs, eigenvecs, num_eigen){
  temp_val <- coefs[1] * eigenvecs[,1]
  if(num_eigen > 1){
    for(i in 2:num_eigen){
      temp_val <- temp_val + (coefs[i] * eigenvecs[,i])
    }
  }
  return(f(temp_val))
}

plot_sf <- function(data_vec, title., lat_lon. = lat_lon){
  temp_df <- data.frame("lat" = lat_lon.$lat, 
                        "lon" = lat_lon.$lon, 
                        "value" = data_vec)
  my_sf <- st_as_sf(temp_df, coords = c("lon", "lat"), crs = 4326)
  
  plot1 <- ggplot(my_sf) + 
    geom_sf(aes(col = value), shape = 15, size = 1.5) + 
    scale_color_viridis_c() + 
    labs(title = title.) + 
    theme(axis.text = element_text(angle = 45)) 
}

###
# eigenvector plots
###
plots1 <- lapply(1:4, plot_eigen, data = eigen_tpdm$vectors, eigen_vals. = eigen_tpdm$values)
plots2 <- lapply(5:8, plot_eigen, data = eigen_tpdm$vectors, eigen_vals. = eigen_tpdm$values)

pdf("./eigenvector_1to8_plots.pdf", height = 10, width = 7.5)
plot_grid(plotlist = plots1, nrow = 2, ncol = 2, align = "hv")
plot_grid(plotlist = plots2, nrow = 2, ncol = 2, align = "hv")
dev.off()

###
# eigenvector plots after transformation
###
transformed_eigenvectors <- f(eigen_tpdm$vectors)
plots1_transformed <- lapply(1:4, plot_eigen, data = transformed_eigenvectors, 
                             eigen_vals. = eigen_tpdm$values)
plots2_transformed <- lapply(5:8, plot_eigen, data = transformed_eigenvectors, 
                             eigen_vals. = eigen_tpdm$values)

pdf("./transformed_eigenvector_plots.pdf", height = 10, width = 7.5)
plot_grid(plotlist = plots1_transformed, nrow = 2, ncol = 2, align = "hv")
plot_grid(plotlist = plots2_transformed, nrow = 2, ncol = 2, align = "hv")
dev.off()

###
# reconstruct storm plots 
###
storm_ind <- which(dates %in% c("20210111", "20210112", "20210113", "20210114"))
storm_plots <- storm_plots2 <- list()


# look at storm to see which day to use
temp_plot <- plot_sf(precip_alpha2[storm_ind[2], ], title. = "storm_20210112", 
                     lat_lon. = lat_lon)
plot(temp_plot) # the 12th looks good 

jan_12 <- precip_alpha2[storm_ind[2], ]

jan_12_coefs <- reconstruction_coefs(data_vec = jan_12, eigenvecs = eigen_tpdm$vectors)

total_scale <- sum(eigen_tpdm$values)
partial_scale <- cumsum(eigen_tpdm$values)
prop_scale <- round(partial_scale/total_scale, 3)

storm_plots[[1]] <- plot_sf(jan_12, title = "Historical Storm") + 
  labs(subtitle = "January 12, 2021")
storm_plots[[2]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 1),
                            title = "1 Eigenvector") + 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[1]))
storm_plots[[3]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 2), 
                            title = "2 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[2]))
storm_plots[[4]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 4), 
                            title = "4 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[4]))
storm_plots2[[1]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 8), 
                             title = "8 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[8]))
storm_plots2[[2]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 20), 
                             title = "20 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[20]))
storm_plots2[[3]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 100), 
                             title = "100 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[100]))
storm_plots2[[4]] <- plot_sf(reconstruct_storm(jan_12_coefs, eigen_tpdm$vectors, 1000), 
                             title = "1000 Eigenvectors")+ 
  labs(subtitle = paste0("Proportion of scale: ", prop_scale[1000]))


pdf("./reconstructed_storm_20211112.pdf", height = 10, width = 7.5)
plot_grid(plotlist = storm_plots, nrow = 2, ncol = 2, align = "hv")
plot_grid(plotlist = storm_plots2, nrow = 2, ncol = 2, align = "hv")

dev.off()


###
# Get time series plots for first 8 eigen PC's
### 
get_ts_plot <- function(index., eigenvec = eigen_tpdm$vectors, 
                        data = precip_alpha2, dates. = dates){
  ts_vec <- data %*% eigenvec[,index.]
  temp_df <- data.frame("Date" = dates., 
                        "PC_Score" = ts_vec, 
                        "x_index" = 1:length(dates.))
  
  x_ticks <- grep("0331", dates.)
  
  ggplot(temp_df, aes(x = x_index, y = PC_Score)) + 
    geom_line() + 
    geom_point(aes(x = 7302, y = PC_Score[7302]), color = "red", shape = 4) + 
    ylab("PC Score") + 
    xlab("Date") + 
    scale_x_continuous(breaks = x_ticks[seq(5, 40, length.out = 8)], 
                       label = substr(temp_df$Date[x_ticks], start = 1, 
                                      stop = 4)[seq(5, 40, length.out = 8)]) + 
    labs(title = paste0("Principle Component ", index.)) + 
    theme(axis.text.x = element_text(angle = 30))
}

ts_plots <- list()

ts_plots <- lapply(c(1:8), get_ts_plot)


pdf("./PC_ts_plots.pdf", height = 10, width = 7.5)
plot_grid(plotlist = ts_plots, nrow = 4, ncol = 2, align = "hv")
dev.off()

