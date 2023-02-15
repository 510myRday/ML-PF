library(MASS)
library(tidyverse)
path <- "C:/Users/Administrator/Desktop/"
file_name <- "real_test_data.csv"
if (file.exists(paste0(path, file_name)) ) {
  read_data <- read.table(paste0(path, file_name), sep = ",")
  exist_num <- min(table(read_data$V29))
} else {
  exist_num <- 0
}
rm(read_data)


# Parameter ----------------------------------------------------------------

attribute_num <- 7
item_num <- 10
person_num <- 500000

Q_matrix <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                     1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                     0, 0, 0, 1, 1, 1, 1, 0, 0, 0), nrow = 10, ncol = 7)

beta_delta_zeta_sigma <- matrix(c(
  -0.653, -5.53, -4.161, -3.428, -0.61, -2.354, -0.955, 0.366, -2.215, -2.805,
  3.941, 5.51, 5.172, 3.158, 1.945, 3.171, 2.326, 0.879, 2.311, 3.121,
  4.223, 4.469, 4.629, 4.778, 3.858, 4.257, 3.739, 4.188, 4.523, 4.376,
  0.524, 0.611, 0.585, 0.411, 0.541, 0.463, 0.499, 0.43, 0.541, 0.63),
  nrow = 10, ncol = 4)

g <- exp(beta_delta_zeta_sigma[, 1]) / (1 + exp(beta_delta_zeta_sigma[, 1]))
s <- 1 / (1 + exp(beta_delta_zeta_sigma[, 1] + beta_delta_zeta_sigma[, 2]))

total_RT_high <- 2400
total_RT_low <- exp(mean(beta_delta_zeta_sigma[,3])) * item_num

lambda_gamma <- matrix(c(1.53, 2.86, -3.353, 2.096, 0.971, 3.963, 0.346, 3.805,
                         1.014, 3.026, 0.37, 3.769, 0.947, 3.853),
                       nrow = 2, ncol = 7)


# Function ----------------------------------------------------------------

run_function <- function(seed, n) {
  
  all_attr_master <- function(attr_num) {
    attr_master <- t(combn(1:2 ^ 4, attr_num)) %% 2 %>% 
      .[!duplicated(.), ]  ## 2^4 only content 9 attribute
    return(attr_master)
  }
  get_person_parameter <- function() {
    person_mean <- c(0, 0)
    person_covariance <- matrix(c(1, -0.23,
                                  -0.23, 0.205), nrow=2, ncol=2)
    theta_tau <- mvrnorm(person_num, person_mean, person_covariance)
    colnames(theta_tau) <- c("θ", "τ")
    dele <- which(theta_tau[, 1] < -3 | theta_tau[, 1] > 3 |
                    theta_tau[, 2] < -3 | theta_tau[, 2] > 3)
    return(theta_tau[-dele, ])
  }
  get_alpha <- function(theta, lambda, gamma){
    p_alpha <- 
      exp(t(tcrossprod(gamma, theta) - lambda)) /
      (1 + exp(t(tcrossprod(gamma, theta) - lambda)))
    
    alpha <- matrix(NA, length(theta), attribute_num)
    pb <- txtProgressBar(1, length(theta), 1, style = 3)
    cat("\n\t\t\t---- Generate_α ----\n")
    for (i in 1:length(theta)) {
      setTxtProgressBar(pb, i)
      for (j in 1:attribute_num) {
        alpha[i, j] <- rbinom(1, 1, p_alpha[i, j])
      }
    }
    colnames(alpha) <- paste0("α", as.character(1:attribute_num))
    close(pb)
    return(alpha)
  }
  get_response_Y <- function() {
    response_Y <- matrix(NA, nrow(alpha), item_num)
    pb <- txtProgressBar(1, nrow(alpha), 1, style = 3)
    cat("\n\t\t\t---- Generate_Y ----\n")
    for (i in 1:nrow(alpha)) {
      setTxtProgressBar(pb, i)
      for (j in 1:item_num) {
        p_response <- g[j] + (1 - s[j] - g[j]) * prod(alpha[i, ] ^ Q_matrix[j, ])
        response_Y[i, j] <- rbinom(1, 1, p_response)
      }
    }
    colnames(response_Y) <- paste0("Y", as.character(1:item_num))
    close(pb)
    return(response_Y)
  }
  get_RT <- function(zeta, tau, sigma) {
    log_RT <- matrix(NA, length(tau), item_num)
    pb <- txtProgressBar(1, length(tau), 1, style = 3)
    cat("\n\t\t   ---- Generate_RT ----\n")
    for (i in 1:length(tau)) {
      setTxtProgressBar(pb, i)
      for (j in 1:item_num) {
        rt <- rnorm(1, zeta[j] - tau[i], sigma[j])
        if(rt < zeta[j] - tau[i] - 3){rt <- zeta[j] - tau[i] - 3}
        if(rt > zeta[j] - tau[i] + 3){rt <- zeta[j] - tau[i] + 3}
        log_RT[i, j] <- rt
      }
    }
    RT <- exp(log_RT) %>% cbind(., rowSums(.))
    name <- paste0("RT", as.character(1:item_num))
    colnames(RT) <- c(name, "total")
    close(pb)
    return(RT)
  }
  Generate_attr_type_and_sum <- function() {
    all_attr_master <- all_attr_master(attribute_num) %>% as.data.frame() %>% 
      unite("merge", everything(), sep = "")
    attr_sum <- rowSums(alpha)
    alpha <- alpha %>%  as.data.frame() %>% 
      unite("merge", everything(), sep = "")
    attr_type <- match(alpha[,1], all_attr_master[,1])
    return(cbind(attr_type, attr_sum))
  }
  get_final_data <- function(attr_type_num) {
    if (attr_type_num == 0) {return()}
    sample <- initial_data %>% as.data.frame() %>%
      filter(attr_type == attr_type_num) %>% 
      .[sample(1:nrow(.), sample_num), ] %>% as.matrix()
    merge <- rbind(sample, get_final_data(attr_type_num - 1))
    return(merge)
  }
  
  print(seed)
  print(n)
  if (seed == 0 | n > 3000) {return()}
  set.seed(seed)
  theta_tau <- get_person_parameter()
  alpha <- get_alpha(theta_tau[, 1], lambda_gamma[1, ], lambda_gamma[2, ])
  Y <- get_response_Y()
  RT <- get_RT(beta_delta_zeta_sigma[, 3], theta_tau[, 2], beta_delta_zeta_sigma[, 4])
  attr_type_and_sum <- Generate_attr_type_and_sum()
  initial_data <- cbind(Y, RT, alpha, attr_type_and_sum) %>% 
    as.data.frame() %>% filter(total > total_RT_low & total < total_RT_high)
  print(table(initial_data[["attr_type"]]))
  if (length(table(initial_data[["attr_type"]])) == 128) {
    sample_num <- min(table(initial_data[["attr_type"]]))
    print(sample_num)
    if (sample_num > 1) {
      final_data <- get_final_data(2 ^ attribute_num)
      write.table(final_data, paste0(path, file_name), sep = ",",
                  append = T, row.names = F, col.names = F)
      n <- n + sample_num
      run_function(seed - 1, n)
    } else {
      run_function(seed - 1, n)
    }
  } else {
    run_function(seed, n)
  }
}


# Run ---------------------------------------------------------------------

run_function(4458, exist_num)

# train_data <- read.table(paste0(path, file_name), sep = ",")
# table(duplicated(train_data$V20))
# train_data <- train_data[-which(duplicated(train_data$V20)), ]
# write.table(final_data, paste0(path, file_name), sep = ",",
#             append = T, row.names = F, col.names = F)

# get_data <- function(attr_type_num) {
#   if (attr_type_num == 0) {return()}
#   sample <- train_data %>% as.data.frame() %>%
#     filter(V29 == attr_type_num) %>% 
#     .[sample(1:nrow(.), 1000 * 3), ] %>% as.matrix()
#   merge <- rbind(sample, get_data(attr_type_num - 1))
#   return(merge)
# }
# train_data <- get_data(128)
# write.table(train_data, paste0(path, file_name), sep = "," , row.names = F,
#             col.names = F)

# set train data and test data --------------------------------------------

library(tidyverse)
path <- "C:/Users/Administrator/Desktop/"
file_name <- "real_test_data.csv"

train_data <- read.table(paste0(path, file_name), sep = ",")

beta_delta_zeta_sigma <- matrix(c(
  -0.653, -5.53, -4.161, -3.428, -0.61, -2.354, -0.955, 0.366, -2.215, -2.805,
  3.941, 5.51, 5.172, 3.158, 1.945, 3.171, 2.326, 0.879, 2.311, 3.121,
  4.223, 4.469, 4.629, 4.778, 3.858, 4.257, 3.739, 4.188, 4.523, 4.376,
  0.524, 0.611, 0.585, 0.411, 0.541, 0.463, 0.499, 0.43, 0.541, 0.63),
  nrow = 10, ncol = 4)

set.seed(1)
aberrant_for_r_guessing <- function(g_data) {
  pb <- txtProgressBar(1, nrow(g_data), 1, style = 3)
  cat("\n\t\t--- Generate_r_guessing ---\n")
  for (i in 1:nrow(g_data)) {
    setTxtProgressBar(pb, i)
    g_data[i, 1:10] <- rbinom(10, 1, 0.25)
    g_data[i, 11:20] <- exp(log(g_data[i, 11:20]) - 3 * beta_delta_zeta_sigma[, 4])
    g_data[i, 21] <- sum(g_data[i, 11:20])
  }
  close(pb)
  return(g_data)
}
normalization <- function(data, mean_vector, sd_vector) {
  if(missing(mean_vector)){mean_vector <- apply(data, 2, mean)}
  if(missing(sd_vector)){sd_vector <- apply(data, 2, sd)}
  for (c in 1:ncol(data)) {
    data[, c] <- (data[, c] - mean_vector[c]) / sd_vector[c]
  }
  return(list(res = data, mean_vector = mean_vector, sd_vector = sd_vector))
}

normal_data <- train_data[seq(1, 384000, by = 3), ] %>% as_tibble() %>% 
  mutate(type = rep("normal", 128000))
speeding_data <- train_data[seq(2, 384000, by = 3), ]
r_guessing_data <- train_data[seq(3, 384000, by = 3), ] %>% 
  aberrant_for_r_guessing() %>% mutate(type = rep("r_guessing", 128000))

output_data <- rbind(normal_data, r_guessing_data)
nor <- output_data[, 1:21] %>% normalization()
train_x <- nor$res %>% as.matrix()
train_mean <- nor$mean_vector
train_sd <- nor$sd_vector

target_fun <- function(behavior_type_vector){
  all_attr_master <- t(combn(1:2 ^ 4, 7)) %% 2 %>% .[!duplicated(.), ]
  ### attribute
  result_attr <- c()
  for (r in 1:nrow(all_attr_master)) {
    attr_str_array <- as.character(all_attr_master[r, ])
    attr_str <- ""
    for (c in 1:ncol(all_attr_master)) {
      attr_str <- paste0(attr_str, attr_str_array[c])
    }
    result_attr <- c(result_attr, attr_str)
  }
  ### combination
  return(
    expand.grid(result_attr, behavior_type_vector) %>%
      unite("target", everything(), sep = "_")
  )
}
target <- c(target_fun("normal")$target, target_fun("r_guessing")$target)

target_y <- c()
pb <- txtProgressBar(1, nrow(output_data), 1, style = 3)
cat("\n\t\t--- Generate_target_y ---\n")
for (i in 1:nrow(output_data)) {
  setTxtProgressBar(pb, i)
  target_y <- c(target_y, paste0(paste(as.character(output_data[i, 22:28]),
                                     collapse = ""), "_", output_data[i, 31]))
}
close(pb)
target_index_train_y <- match(target_y, target)

train_y <- matrix(0, nrow = nrow(output_data), ncol = length(target))
for (r in 1:nrow(output_data)) {
  train_y[r, target_index_train_y[r]] <- 1
}

# write.table(train_x, paste0(path, "train_x.csv"), sep = ",",
#             append = T, row.names = F, col.names = F)
# write.table(train_y, paste0(path, "train_y.csv"), sep = ",",
#             append = T, row.names = F, col.names = F)

write.table(train_x, paste0(path, "test_x.csv"), sep = ",",
            append = T, row.names = F, col.names = F)
write.table(train_y, paste0(path, "test_y.csv"), sep = ",",
            append = T, row.names = F, col.names = F)




# plot --------------------------------------------------------------------
library(tidyverse)
path <- "C:/Users/Administrator/Desktop/"
note <- "(NOTE: \"Y\" means only response data, \"RT\" means only response time data, and \"Y&RT\" means both response data and response time data. The title of each chart means: the 
combination of the number of attributes, the number of items, the quality of items, and the proportion of aberrant items. The \"class\" on the x-axis means that the response behavior 
prediction is correct but the attribute mastery pattern is not necessarily correct. The \"single\" on the x-axis means that the response behavior and attribute mastering pattern are 
both correctly predicted. The \"total\" on the x-axis means that normal or aberrant responses are correctly predicted, ignoring the response behavior and attribute mastery pattern.)"

real_res <- rbind(c("class", "normal", "RT", 1), c("class", "r_guessing", "RT", 1),
                  c("single", "normal", "RT", 0.01), c("single", "r_guessing", "RT", 0.005),
                  c("total", "normal", "RT", 1), c("total", "r_guessing", "RT", 1),
                  c("class", "normal", "Y", 0.963), c("class", "r_guessing", "Y", 0.348),
                  c("single", "normal", "Y", 0.063), c("single", "r_guessing", "Y", 0.003),
                  c("total", "normal", "Y", 0.963), c("total", "r_guessing", "Y", 0.348),
                  c("class", "normal", "Y & RT", 1), c("class", "r_guessing", "Y & RT", 0.999),
                  c("single", "normal", "Y & RT", 0.068), c("single", "r_guessing", "Y & RT", 0.015),
                  c("total", "normal", "Y & RT", 1), c("total", "r_guessing", "Y & RT", 0.999)) %>% 
  as.data.frame()
colnames(real_res) <- c("Result", "type", "Input_features", "Accuracy")
real_res$Accuracy <- as.numeric(real_res$Accuracy)

png(paste0(path, "Simulate real data", ".png"), width = 3840, height = 2160,
    units = "px", res = 270)
ggplot(real_res,
       aes(x = Result, y = Accuracy, group = Input_features, color = Input_features)) +
  geom_point() + geom_line() + facet_wrap(~ type) +
  theme(plot.title = element_text(color="#c04851", size=20, hjust=0.5,
                                  vjust=0.5, angle=360),
        legend.position = "right",
        legend.background = element_rect(fill = "white",colour = "black")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs(x = paste0("Result\n\n", note), y = "Accuracy", title = "Simulate real data",
       col = "Input_features: ")
dev.off()













