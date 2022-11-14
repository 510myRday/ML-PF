library(tidyverse)
library(MASS)

# 导入Q矩阵 ------------------------------------------------------------------
Q_matrix <- read.csv(file.choose(), header = F)[ , 1:25]

Q_3_10 <- as.matrix(Q_matrix[1:10, 1:3]) %>% 
  .[order(rowSums(.)), ]
Q_3_20 <- as.matrix(Q_matrix[1:20, 5:7]) %>% 
  .[order(rowSums(.)), ]
Q_3_40 <- as.matrix(Q_matrix[1:40, 9:11]) %>% 
  .[order(rowSums(.)), ]
Q_6_20 <- as.matrix(Q_matrix[1:20, 13:18]) %>% 
  .[order(rowSums(.)), ]
Q_6_40 <- as.matrix(Q_matrix[1:40, 20:25]) %>% 
  .[order(rowSums(.)), ]
rm(Q_matrix)

# 可调参数 -------------------------------------------------------------------
file_name <- "Q_3_10_high_0.2_"
Q_matrix <- Q_3_10
item_quality <- "high"
abberant_ratio <- 0.2

zeta_low <- 3
zeta_high <- 4.5
cheating_sub_lrt <- 3
sp_gu_sub_lrt <- 3
train_each_num <- 1000
test_each_num <- 500
behavior_type_vector <- c("normal", "speeding", "cheating", "r_guessing")
x_y_path <- "C:/Users/Administrator/Desktop/First/DATA2/"

item_num <- nrow(Q_matrix)
attribute_num <- ncol(Q_matrix)
total_RT_high <- exp(zeta_high) * item_num
total_RT_low <- exp(zeta_low) * item_num

ifelse(attribute_num == 3, person_num <- 15000, person_num <- 150000)

sigma_epsilon <- runif(item_num, 0.4, 0.6)

gamma_vector <- rep(1, attribute_num)
lambda_vector <- rep(0, attribute_num)

# 自定义函数 -----------------------------------------------------------------
#得到所有属性掌握模式
all_attr_master <- function(attr_num) {
  attr_master <- t(combn(1:2 ^ 4, attr_num)) %% 2 %>% 
    .[!duplicated(.), ]  ## 2^4 only content 9 attribute
  return(attr_master)
}

#得到输出层标签类别
target <- function(all_attr_master, behavior_type_vector){
  #attribute
  result_attr <- c()
  for (r in 1:nrow(all_attr_master)) {
    attr_str_array <- as.character(all_attr_master[r, ])
    attr_str <- ""
    for (c in 1:ncol(all_attr_master)) {
      attr_str <- paste0(attr_str, attr_str_array[c], "_")
    }
    result_attr <- c(result_attr, attr_str)
  }
  #combination
  return(
    expand.grid(result_attr, behavior_type_vector) %>%
      unite("target", everything(), sep = "")
    )
}

#得到θ和τ
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

#得到β, δ, ζ, g, s
get_item_parameter <- function() {
  if (item_quality == "low") {
    gs_mean <- 0.25
    cut_value <- 0.4
  } else if (item_quality == "high") {
    gs_mean <- 0.1
    cut_value <- 0.2
  }
  mean_beta <- log(gs_mean / (1 - gs_mean))
  mean_delta <- log((1 - gs_mean) / gs_mean) - mean_beta
  mean_zeta <- 3
  
  item_mean <- c(mean_beta, mean_delta, mean_zeta)
  item_covariance <- matrix(c(3.863, -2.454, -0.444,
                              -2.454, 2.501, 0.25,
                              -0.444, 0.25, 0.242), nrow=3, ncol=3)
  beta_delta_zeta_g_s <- matrix(NA, 1, 5)
  while (nrow(beta_delta_zeta_g_s) < item_num + 1) {
    one <- mvrnorm(1, item_mean, item_covariance)
    g <- exp(one[1]) / (1 + exp(one[1]))
    s <- 1 - exp(one[1] + one[2]) / (1 + exp(one[1] + one[2]))
    if (cut_value > g & cut_value > s &
        one[3] < zeta_high & one[3] > zeta_low) {
      one <- c(one, g, s)
      beta_delta_zeta_g_s <- rbind(beta_delta_zeta_g_s, one)
    }
  }
  beta_delta_zeta_g_s <- beta_delta_zeta_g_s[-1, ] %>% .[order(.[, 3]), ]
  colnames(beta_delta_zeta_g_s) <- c("β", "δ", "ζ", "g", "s")
  rownames(beta_delta_zeta_g_s) <- NULL
  return(beta_delta_zeta_g_s)
}

#得到所有被试的真值α
get_alpha <- function(theta){
  p_alpha <- 
    exp(t(tcrossprod(gamma_vector, theta) - lambda_vector)) /
    (1 + exp(t(tcrossprod(gamma_vector, theta) - lambda_vector)))
  
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

#得到所有被试在所有题目上的得分Y
get_response_Y <- function(g, s, alpha) {
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

#得到所有被试在所有题目上的作答时间和全卷总时间，单位s
get_RT <- function(zeta, tau) {
  log_RT <- matrix(NA, length(tau), item_num)
  pb <- txtProgressBar(1, length(tau), 1, style = 3)
  cat("\n\t\t   ---- Generate_RT ----\n")
  for (i in 1:length(tau)) {
    setTxtProgressBar(pb, i)
    for (j in 1:item_num) {
      rt <- rnorm(1, zeta[j] - tau[i], sigma_epsilon[j])
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

#得到被试alpha在所有掌握模式当中的索引
Generate_attr_type_and_sum <- function(alpha) {
  all_attr_master <- all_attr_master(attribute_num) %>% as.data.frame() %>% 
    unite("merge", everything(), sep = "")
  attr_sum <- rowSums(alpha)
  alpha <- alpha %>%  as.data.frame() %>% 
    unite("merge", everything(), sep = "")
  attr_type <- match(alpha[,1], all_attr_master[,1])
  return(cbind(attr_type, attr_sum))
}

#整理后得到用于抽取样本的数据集
get_initial_data <- function(item_parameter) {
  theta_tau <- get_person_parameter()
  alpha <- get_alpha(theta_tau[, 1])
  Y <- get_response_Y(item_parameter[, 4], item_parameter[, 5], alpha)
  RT <- get_RT(item_parameter[, 3], theta_tau[, 2])
  attr_type_and_sum <- Generate_attr_type_and_sum(alpha)
  merge <- cbind(Y, RT, alpha, attr_type_and_sum) %>% 
    as.data.frame() %>% filter(total > total_RT_low & total < total_RT_high)
  print(table(merge[["attr_type"]]))
  return(merge)
}

#递归函数，从初始数据集中给每种属性掌握模式抽取一定样本，如：500，1000
get_final_data <- function(initial_data, attr_type_num, sample_num) {
  if (attr_type_num == 0) {return()}
  sample <- initial_data %>% as.data.frame() %>%
    filter(attr_type == attr_type_num) %>% 
    .[sample(1:nrow(.), sample_num), ] %>% as.matrix()
  merge <- rbind(sample,
                 get_final_data(initial_data, attr_type_num - 1, sample_num))
  return(merge)
}

#根据不同的作答类型生成对应的数据集, 并递归合并
get_dif_type_data <- function(data_type, each_num) {
  if(length(data_type) == 0){
    return()
  }
  
  initial <- get_initial_data(beta_delta_zeta_g_s)
  final <- get_final_data(initial, 2 ^ attribute_num, each_num)
  
  switch (data_type[1],
          
          speeding = {
            final <- final %>% as_tibble() %>%
              filter(1 < attr_sum) %>% mutate(α_type = "speeding") %>%
              unite("target", contains("α"))
            ab_num <- abberant_ratio * item_num
            for (n in (item_num - ab_num + 1):item_num) {
              final[[paste0("Y", as.character(n))]] <-
                rbinom(nrow(final), 1, 0.25)
              final[[paste0("RT", as.character(n))]] <-
                exp(log(final[[paste0("RT", as.character(n))]]) - 
                      sp_gu_sub_lrt * 0.5)}
          },
          cheating = {
            final <- final %>% as_tibble() %>%
              filter(attr_sum < attribute_num) %>%
              mutate(α_type = "cheating") %>% unite("target", contains("α"))
            ab_num <- abberant_ratio * item_num
            for (n in (item_num - ab_num + 1):item_num) {
              final[[paste0("Y", as.character(n))]] <-
                rbinom(nrow(final), 1, 0.95)
              final[[paste0("RT", as.character(n))]] <-
                exp(log(final[[paste0("RT", as.character(n))]]) - 
                      cheating_sub_lrt * 0.5)}
          },
          r_guessing = {
            final <- final[sample(1:nrow(final), each_num), ] %>% as_tibble()
            for (n in 1:item_num) {
              final[[paste0("Y", as.character(n))]] <-
                rbinom(nrow(final), 1, 0.25)
              final[[paste0("RT", as.character(n))]] <-
                exp(log(final[[paste0("RT", as.character(n))]]) - 
                      sp_gu_sub_lrt * 0.5)}
            final <- final %>% as_tibble() %>%
              mutate(α_type = "r_guessing") %>% unite("target", contains("α"))
          },
          normal = {
            final <- final %>% as_tibble() %>%
              mutate(α_type = "normal") %>% unite("target", contains("α"))
          })
  
  final[["total"]] <- rowSums(dplyr::select(final, contains("RT")))
  data_type <- data_type[-1]
  print(length(data_type))
  return(final %>% as_tibble() %>% dplyr::select(-contains("attr")) %>% 
           rbind(get_dif_type_data(data_type, each_num)))
}

# 主程序 ---------------------------------------------------------------------
#↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
#↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓改一次条件运行以下一次↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

#得到所有组合类别
output_target <- target(all_attr_master(attribute_num), behavior_type_vector)
#重复50次得到50个数据
for (seed in 1:2) {
  set.seed(seed)
  
  #item_parameter
  beta_delta_zeta_g_s <- get_item_parameter()
  
  #train
  train_data <- get_dif_type_data(behavior_type_vector, train_each_num)
  #train_x
  train_x_T <- dplyr::select(train_data, contains("Y") | contains("RT")) %>% 
    as.matrix()
  train_x_F <- dplyr::select(train_data, contains("Y")) %>% as.matrix()
  #train_y
  target_index_train <- match(train_data$target, output_target$target)
  train_y <- matrix(0, nrow(train_data), length(output_target$target))
  for (r in 1:nrow(train_data)) {
    train_y[r, target_index_train[r]] <- 1
  }
  colnames(train_y) <- output_target$target
  #写出train
  write.csv(train_x_T, paste0(x_y_path, as.character(seed), "s_", file_name,
                              "train_x_T", ".csv"))
  write.csv(train_x_F, paste0(x_y_path, as.character(seed), "s_", file_name,
                              "train_x_F", ".csv"))
  write.csv(train_y, paste0(x_y_path, as.character(seed), "s_", file_name,
                            "train_y", ".csv"))
  
  #test
  test_data <- get_dif_type_data(behavior_type_vector, test_each_num)
  for (i in behavior_type_vector) {
    test_type <- i
    #获取作答类型索引
    index <- grep(test_type, test_data[["target"]])
    #test_x
    test_x_T <- dplyr::select(test_data, contains("Y") | contains("RT")) %>% 
      .[index, ] %>% as.matrix()
    test_x_F <- dplyr::select(test_data, contains("Y")) %>% .[index, ] %>% 
      as.matrix()
    #test_y
    target_index_test <- match(test_data$target, output_target$target)
    test_y <- matrix(0, nrow(test_data), length(output_target$target))
    for (r in 1:nrow(test_data)) {
      test_y[r, target_index_test[r]] <- 1
    }
    test_y <- test_y[index, ]
    colnames(train_y) <- output_target$target
    #写出test
    write.csv(test_x_T, paste0(x_y_path, as.character(seed), "s_", file_name,
                               "test_", test_type, "_x_T", ".csv"))
    write.csv(test_x_F, paste0(x_y_path, as.character(seed), "s_", file_name,
                               "test_", test_type, "_x_F", ".csv"))
    write.csv(test_y, paste0(x_y_path, as.character(seed), "s_", file_name,
                             "test_", test_type, "_y", ".csv"))
  }
  print(paste0("***********", "ending ", as.character(seed), "***********"))
}

#↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
#↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑




