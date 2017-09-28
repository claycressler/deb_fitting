library(magrittr)
x <- readRDS("TM_parameter_estimates_deb_parasite_E_TypeI.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik

x <- readRDS("TM_parameter_estimates_deb_parasite_E_TypeI_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_E_TypeII.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_E_TypeII_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_E_logistic.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]


x <- readRDS("TM_parameter_estimates_deb_parasite_W_TypeI.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_W_TypeI_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_W_TypeII.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_W_TypeII_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_W_logistic.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]


x <- readRDS("TM_parameter_estimates_deb_parasite_K_constant.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_K_constant_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_K_constant_logistic.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]


x <- readRDS("TM_parameter_estimates_deb_parasite_K_variable.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_K_variable_m.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]

x <- readRDS("TM_parameter_estimates_deb_parasite_K_variable_logistic.RDS")
x %>% unlist %>% matrix(., ncol=length(x[[1]]$params)+2, byrow=TRUE) %>% as.data.frame -> y
colnames(y) <- c(names(x[[1]]$params), 'lik', 'conv')
y <- y[order(y$lik),]
y$aic <- 2*length(x[[1]]$params) + 2*y$lik
y[1,]


