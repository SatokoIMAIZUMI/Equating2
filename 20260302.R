setwd("/Users/imaizumisatoko/Documents/研究/AMPC2026") 

rm(list = ls())

# ------------------------------------------------------------------------------
# 1. ライブラリと設定
# ------------------------------------------------------------------------------
packages <- c("mirt", "numDeriv", "ggplot2", "dplyr", "gridExtra")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
}

set.seed(425)


CONFIG <- list(
  N_EXAMINEES = 5000,
  N_REPS = 100,
  N_UNIQUE = 20,
  TRUE_A = 1.2,
  TRUE_B = 0.5,
  D_CONST = 1.702,
  DRIFT_MAG = 1.0,       # ドリフト量(b)
  OUTLIER_PROP = 0.2    # 外れ値の割合 
)

# ------------------------------------------------------------------------------
# 2. 関数定義 (モデル・損失関数・推定)
# ------------------------------------------------------------------------------

# 2PLM Probability
p_func <- function(theta, a, b, D=1.702) {
  1 / (1 + exp(-D * a * (theta - b)))
}

# Equated Probability
pi_func <- function(theta, a, b, A, B, D=1.702) {
  if(abs(A) < 1e-5) A <- sign(A) * 1e-5
  p_func(theta, a/A, A*b + B, D)
}

# Divergence (Element-wise)
calc_div_element <- function(P, pi_val, div_type="DPD", param=0.1) {
  #計算安定性のため
  eps <- 1e-9
  P <- pmax(pmin(P, 1-eps), eps)
  pi_val <- pmax(pmin(pi_val, 1-eps), eps)
  
  if(div_type == "DPD") { 
    beta <- param#tuning parameter
    term1 <- -(1/beta) * (P * pi_val^beta + (1-P) * (1-pi_val)^beta)
    term2 <- (1/(1+beta)) * (pi_val^(1+beta) + (1-pi_val)^(1+beta))
    return(term1 + term2)
  } else if(div_type == "gamma") { 
    gamma <- param
    term1 <- -(1/gamma) * log(P * pi_val^gamma + (1-P) * (1-pi_val)^gamma)
    term2 <- (1/(1+gamma)) * log(pi_val^(1+gamma) + (1-pi_val)^(1+gamma))
    return(term1 + term2)
  } else { 
    # KL (Approx)
    return(P * log(P/pi_val) + (1-P) * log((1-P)/(1-pi_val)))
  }
}

# Objective Function（全項目，全求積点の損失）
obj_func_total <- function(eta, a_new, b_new, P_ref, nodes, weights, div_type, param, D) {
  A <- eta[1]; B <- eta[2]
  loss <- 0
  for(j in 1:length(a_new)) {
    pi_vec <- pi_func(nodes, a_new[j], b_new[j], A, B, D)#等化後正答確率
    div_vec <- calc_div_element(P_ref[,j], pi_vec, div_type, param)#損失
    loss <- loss + sum(div_vec * weights)
  }
  return(loss)
}


# Estimation with Sandwich Variance(SE cal)
estimate_equating_sandwich <- function(a_new, b_new, P_ref, nodes, weights, div_type, param, D=1.702, start_val=c(1, 1)){
  
  
  
  # 1. Point Estimation (点推定)
  opt <- tryCatch({
    optim(par = start_val, fn = obj_func_total,
          a_new = a_new, b_new = b_new, P_ref = P_ref, 
          nodes = nodes, weights = weights, div_type = div_type, param = param, D = D,
          method = "BFGS", control = list(maxit = 5000))
  }, error = function(e) NULL)
  
  if(is.null(opt) || opt$convergence != 0) return(list(conv=FALSE))
  
  hat_A <- opt$par[1]; hat_B <- opt$par[2]
  
  # 2. Sandwich Variance: V = H^-1 G H^-1
  # 解析解を用いて H と G を同時に計算する (高速化)
  
  J <- length(a_new)
  M <- length(nodes)
  
  # 行列の初期化
  H <- matrix(0, 2, 2) # Bread (Hessian)
  G <- matrix(0, 2, 2) # Meat (Outer product of scores)
  
  eps <- 1e-9 # 数値安定化用
  
  for(j in 1:J) {
    aj <- a_new[j]
    bj <- b_new[j]
    
    # Random Item Theory: 項目単位のスコアベクトルを初期化
    Psi_j_sum <- c(0, 0) 
    
    for(m in 1:M) {
      theta <- nodes[m]
      wt <- weights[m]
      P_val <- P_ref[m, j]
      
      # --- 共通計算: pi と z ---
      z <-  (D * aj / hat_A) * (theta - hat_A * bj - hat_B)
      pi_val <- 1 / (1 + exp(-z))
      pi_val <- pmax(pmin(pi_val, 1-eps), eps) # ガード
      p1p <- pi_val * (1 - pi_val)
      
      # --- 勾配用パーツ (1階微分) ---
      dz_dA <- (D * aj / hat_A^2) * (theta - hat_B)
      dz_dB <- (D * aj / hat_A)
      dz_d_eta <- c(dz_dA, dz_dB)
      
      d_pi_d_eta <- p1p * dz_d_eta
      
      # --- ヘッセ行列用パーツ (2階微分) ---
      # zの2階微分
      d2z_dA2 <- -2 * (D * aj / hat_A^3) * (theta - hat_B)
      d2z_dAdB <- - (D * aj / hat_A^2)
      d2z_d_eta2 <- matrix(c(d2z_dA2, d2z_dAdB, d2z_dAdB, 0), 2, 2)
      
      # piの2階微分
      term_grad_prod <- p1p * (1 - 2 * pi_val) * (dz_d_eta %*% t(dz_d_eta))
      term_hess_z    <- p1p * d2z_d_eta2
      d2_pi_d_eta2   <- term_grad_prod + term_hess_z
      
      # --- ダイバージェンスの微分 ---
      dD_dpi <- 0
      d2D_dpi2 <- 0
      
      if (div_type == "DPD") {
        beta <- param
        # 1階
        term_pow1 <- (1 - pi_val)^(beta - 1) + pi_val^(beta - 1)
        dD_dpi <- (pi_val - P_val) * term_pow1
        # 2階
        d_term_pow1 <- (beta - 1) * (pi_val^(beta - 2) - (1 - pi_val)^(beta - 2))
        d2D_dpi2 <- term_pow1 + (pi_val - P_val) * d_term_pow1
        
      } else if (div_type == "gamma") {
        gamma <- param
        # Gamma定義
        pi_g_1  <- pi_val^(gamma - 1); ipi_g_1 <- (1 - pi_val)^(gamma - 1)
        pi_g    <- pi_val * pi_g_1;    ipi_g   <- (1 - pi_val) * ipi_g_1
        
        N1 <- P_val * pi_g_1 - (1 - P_val) * ipi_g_1
        D1 <- P_val * pi_g + (1 - P_val) * ipi_g
        N2 <- pi_g - ipi_g
        D2 <- pi_val * pi_g + (1 - pi_val) * ipi_g
        
        # 1階
        dD_dpi <- - (N1 / D1) + (N2 / D2)
        
        # 2階 (商の微分)
        pi_g_2  <- pi_val^(gamma - 2); ipi_g_2 <- (1 - pi_val)^(gamma - 2)
        dN1 <- (gamma - 1) * (P_val * pi_g_2 + (1 - P_val) * ipi_g_2)
        dD1 <- gamma * N1
        term1_2nd <- - (dN1 * D1 - N1 * dD1) / (D1^2)
        
        dN2 <- gamma * (pi_g_1 + ipi_g_1)
        dD2 <- (1 + gamma) * N2
        term2_2nd <- (dN2 * D2 - N2 * dD2) / (D2^2)
        
        d2D_dpi2 <- term1_2nd + term2_2nd
        
      } else { 
        # KL
        dD_dpi <- (pi_val - P_val) / p1p
        d2D_dpi2 <- (p1p - (pi_val - P_val)*(1 - 2*pi_val)) / (p1p^2)
      }
      
      # --- H行列への加算 (Chain Rule) ---
      # A_element = (d2D/dpi2) * (d_pi)(d_pi)^T + (dD/dpi) * (d2_pi)
      term1 <- d2D_dpi2 * (d_pi_d_eta %*% t(d_pi_d_eta))
      term2 <- dD_dpi * d2_pi_d_eta2
      H <- (H + (term1 + term2) * wt)
      
      # --- 勾配ベクトル (Psi) の蓄積 ---
      # psi_jm = dD/dpi * d_pi/d_eta * w
      psi_jm <- (dD_dpi * d_pi_d_eta) * wt　#これだともしかして重みが2乗されるかも
      Psi_j_sum <- Psi_j_sum + psi_jm
    }
    
    # --- G行列への加算 (Random Item Theory) ---
    # 項目単位のスコアの合計の外積をとる
    G <- (G + (Psi_j_sum %*% t(Psi_j_sum)))
  }
  
  #分散の計算
  H <- H/J
  G <- G/J
  
  # 3. 分散共分散行列の計算
  #共通項目数少ないときのランク落ち(frag立てる）------------------------------
  
  inv_H <- try(solve(H), silent=TRUE)
  
  if(inherits(inv_H, "try-error")) {
    warning(sprintf("H の逆行列計算失敗 (J=%d)", J))  # ← 追加
    CovMat <- matrix(NA, 2, 2)
    SE <- c(NA, NA)
  } else {
    # 条件数チェックを追加
    cond_num <- kappa(H)
    if(cond_num > 1e10) warning(sprintf("H が ill-conditioned (kappa=%.2e, J=%d)", cond_num, J))
    
    CovMat <- inv_H %*% G %*% inv_H / J
    SE <- sqrt(diag(CovMat))
  }
  return(list(A=hat_A, B=hat_B, SE_A=SE[1], SE_B=SE[2], CovMat=CovMat, conv=TRUE))
}

# ------------------------------------------------------------------------------
# 3. シミュレーション実行コア (1反復分)
# ------------------------------------------------------------------------------
run_simulation_core <- function(J, M, div_type, param, has_outlier=FALSE){
  
  N <- CONFIG$N_EXAMINEES
  N_uniq <- CONFIG$N_UNIQUE
  valid_data <- FALSE
  retry_count <- 0
  
  while(!valid_data && retry_count < 10) {# 全問正解，不正解がでたら再生成
    retry_count <- retry_count + 1
    
    # Param Generation
    a_com_ref <- rlnorm(J, 0, 0.3); b_com_ref <- rnorm(J, 0, 1)
    a_com_new_true <- a_com_ref * CONFIG$TRUE_A
    b_com_new_true <- (b_com_ref - CONFIG$TRUE_B) / CONFIG$TRUE_A
    
    ## 外れ値-------------------------------------------------
    if(has_outlier){#外れ値の導入
      n_drift <- floor(J * CONFIG$OUTLIER_PROP)
      if(n_drift > 0){
        idx <- sample(1:J, n_drift)
        b_com_new_true[idx] <- b_com_new_true[idx] + CONFIG$DRIFT_MAG
        a_com_new_true[idx] <- a_com_new_true[idx]*1.2
      }
    }
    #独自項目の生成
    a_uniq_ref <- rlnorm(N_uniq, 0, 0.3); b_uniq_ref <- rnorm(N_uniq, 0, 1)
    a_uniq_new <- rlnorm(N_uniq, 0, 0.3); b_uniq_new <- rnorm(N_uniq, 0, 1)
    #全項目を結合
    a_ref_total <- c(a_com_ref, a_uniq_ref); b_ref_total <- c(b_com_ref, b_uniq_ref)
    a_new_total <- c(a_com_new_true, a_uniq_new); b_new_total <- c(b_com_new_true, b_uniq_new)
    
    # Data Generation
    #能力値の生成
    th_ref <- rnorm(N); th_new <- rnorm(N)
    #反応データ生成：正答確率の計算（受験者×項目），一様乱数との比較，正誤の決定
    resp_ref <- 1 * (t(sapply(th_ref, function(t) p_func(t, a_ref_total, b_ref_total))) > matrix(runif(N*(J+N_uniq)), N, J+N_uniq))
    resp_new <- 1 * (t(sapply(th_new, function(t) p_func(t, a_new_total, b_new_total))) > matrix(runif(N*(J+N_uniq)), N, J+N_uniq))
    colnames(resp_ref) <- colnames(resp_new) <- paste0("I", 1:(J+N_uniq))
    #全問正解，不正解の項目チェック
    if(all(apply(resp_ref, 2, var) > 0) && all(apply(resp_new, 2, var) > 0)) valid_data <- TRUE
  }
  if(!valid_data) return(NULL)
  
  # IRT Calibration--------------------------------
  # mirt_res <- try({
  #   invisible(capture.output({
  #     mod_ref <- mirt(as.data.frame(resp_ref), 1, itemtype="2PL", verbose=FALSE, technical=list(NCYCLES=500))
  #     mod_new <- mirt(as.data.frame(resp_new), 1, itemtype="2PL", verbose=FALSE, technical=list(NCYCLES=500))
  #   }))
  #   list(ref=mod_ref, new=mod_new)
  # }, silent=TRUE)
  # 
  # if(inherits(mirt_res, "try-error")) return(NULL)
  # #共通項目のパラメータ抽出
  # cfs_ref <- coef(mirt_res$ref, IRTpars=TRUE, simplify=TRUE)$items
  # cfs_new <- coef(mirt_res$new, IRTpars=TRUE, simplify=TRUE)$items
  # 
  # # Extract Common Items
  # est_a_ref_com <- cfs_ref[1:J, "a"] / CONFIG$D_CONST; est_b_ref_com <- cfs_ref[1:J, "b"]
  # est_a_new_com <- cfs_new[1:J, "a"] / CONFIG$D_CONST; est_b_new_com <- cfs_new[1:J, "b"]
  
  #オラクル条件---------------------------
  est_a_ref_com <- a_com_ref
  est_b_ref_com <- b_com_ref
  est_a_new_com <- a_com_new_true
  est_b_new_com <- b_com_new_true
  # Equating
  nodes <- seq(-4, 4, length.out=M)
  weights <- dnorm(nodes); weights <- weights/sum(weights)#正規分布の重み，正規化
  #基準テストの確率行列
  P_ref <- matrix(0, M, J)
  for(j in 1:J) P_ref[,j] <- p_func(nodes, est_a_ref_com[j], est_b_ref_com[j], CONFIG$D_CONST)
  #等化係数推定
  res <- estimate_equating_sandwich(est_a_new_com, est_b_new_com, P_ref, nodes, weights, div_type, param, CONFIG$D_CONST)
  
  if(res$conv) return(data.frame(res[c("A","B","SE_A","SE_B")], J=J, M=M)) else return(NULL)
}

# ==============================================================================
# 初期値依存性の確認 
# ==============================================================================
# cat("\n【特別検証】初期値依存性のチェックを実行します...\n")
# 
# check_initial_values <- function() {
#   # 1. テスト用のクリーンなデータ（外れ値なし）を1セット生成
#   J_test <- 50
#   M_test <- 301
#   nodes <- seq(-4, 4, length.out=M_test)
#   weights <- dnorm(nodes); weights <- weights/sum(weights)
#   
#   a_ref <- rlnorm(J_test, 0, 0.3)
#   b_ref <- rnorm(J_test, 0, 1)
#   a_new <- a_ref * CONFIG$TRUE_A
#   b_new <- (b_ref - CONFIG$TRUE_B) / CONFIG$TRUE_A
#   
#   P_ref <- matrix(0, M_test, J_test)
#   for(j in 1:J_test) P_ref[,j] <- p_func(nodes, a_ref[j], b_ref[j], CONFIG$D_CONST)
#   
#   # 2. 初期値のグリッドを作成 
#   A_init_seq <- c(0.5, 1.0, 1.5, 2.0)
#   B_init_seq <- c(-1.0, 0.0, 2.0)
#   init_grid <- expand.grid(A_init = A_init_seq, B_init = B_init_seq)
#   
#   results <- list()
#   
#   # 3. 各初期値の組み合わせで最適化を実行
#   for(i in 1:nrow(init_grid)) {
#     init_val <- c(init_grid$A_init[i], init_grid$B_init[i])
#     
#     # DPD(0.1)を例として実行
#     res <- estimate_equating_sandwich(a_new, b_new, P_ref, nodes, weights, 
#                                       div_type="DPD", param=0.1, 
#                                       D=CONFIG$D_CONST, start_val=init_val)
#     
#     if(res$conv) {
#       results[[i]] <- data.frame(
#         Init_A = init_val[1], Init_B = init_val[2],
#         Est_A = res$A, Est_B = res$B,
#         Conv = TRUE
#       )
#     } else {
#       results[[i]] <- data.frame(
#         Init_A = init_val[1], Init_B = init_val[2],
#         Est_A = NA, Est_B = NA,
#         Conv = FALSE
#       )
#     }
#   }
#   
#   res_df <- bind_rows(results)
#   
#   # コンソールに結果を表示
#   print(res_df)
#   
#   # 4. 可視化: 横軸を初期値、縦軸を推定値Aとするプロット
#   res_df <- res_df %>% mutate(Init_Label = paste0("(", Init_A, ",", Init_B, ")"))
#   
#   p_init <- ggplot(res_df, aes(x=Init_Label, y=Est_B)) +
#     geom_point(size=4, color="blue", alpha=0.7) +
#     geom_hline(yintercept=CONFIG$TRUE_A, linetype="dashed", color="red", size=1) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title="Initial Value Dependence Check (Target A = 1.2)", 
#          subtitle="If stable, all points should align horizontally on the red dashed line.",
#          x="Initial Values (A, B)", y="Estimated A")
#   
#   print(p_init)
# }
# 
# # 実行
# check_initial_values()
# ==============================================================================
# ------------------------------------------------------------------------------
# 4. 実験メインループ (Monte Carlo Simulation)
# ------------------------------------------------------------------------------
results_list <- list()
counter <- 1

J_vec <- c(5, 10, 20, 500, 1000)#共通項目数
M_vec <- c(301)# 求積点数
conditions <- expand.grid(J=J_vec, M=M_vec)#全組み合わせ

methods <- list(
  #list(name="KL (Approx)", type="DPD", p=0.01),
  list(name="DPD (0.1)",   type="DPD", p=0.1),
  list(name="DPD (0.3)",   type="DPD", p=0.3),
  list(name="DPD (0.5)",   type="DPD", p=0.5),
  list(name="DPD (1.0)",   type="DPD", p=1.0),
  #list(name="DPD (1.0)",   type="DPD", p=1.0),
  #list(name="Gamma (0.1)", type="gamma", p=0.1),
  list(name="Gamma (0.3)", type="gamma", p=0.3),
  list(name="Gamma (0.5)", type="gamma", p=0.5),
  list(name="Gamma (1.0)", type="gamma", p=1.0)
)

outlier_conds <- c(TRUE,FALSE)

cat("【Step 1/2】統計的性質の検証シミュレーションを開始します...\n")
cat(sprintf("  総条件数: %d x 3手法 x 2状況 (反復各%d回)\n", nrow(conditions), CONFIG$N_REPS))

for(has_outlier in outlier_conds) {
  cond_label <- if(has_outlier) "With_Drift" else "No_Outlier"
  cat(sprintf("\n=== Condition: %s ===\n", cond_label))
  
  for(met in methods) {
    cat(sprintf("  Method: %s ", met$name))
    for(k in 1:nrow(conditions)) {
      J_curr <- conditions$J[k]; M_curr <- conditions$M[k]
      for(i in 1:CONFIG$N_REPS){
        res <- run_simulation_core(J_curr, M_curr, met$type, met$p, has_outlier)
        if(!is.null(res)){
          res$Method <- met$name; res$Condition <- cond_label
          results_list[[counter]] <- res; counter <- counter + 1
        }
      }
      cat(".")
    }
    cat(" Done\n")
  }
}

all_results <- bind_rows(results_list)

# ------------------------------------------------------------------------------
# 5. 集計と可視化 (統計的性質)
# ------------------------------------------------------------------------------
cat("\n【Step 2/2】結果の集計とグラフ作成中...\n")

library(dplyr)

# 全ての統計量（SE_Ratioを含む）を再計算
summary_stats <- all_results %>%
  group_by(Condition, Method, J, M) %>%
  summarise(
    Mean_A = mean(A),
    Bias_A = mean(A - CONFIG$TRUE_A),
    RMSE_A = sqrt(mean((A - CONFIG$TRUE_A)^2)),
    Mean_SE = mean(SE_A, na.rm=TRUE),
    SD_Est = sd(A),
    
    
    SE_Ratio = mean(SE_A, na.rm=TRUE) / sd(A),
    
    Coverage = mean(abs(A - CONFIG$TRUE_A) < 1.96 * SE_A, na.rm=TRUE),
    
    Mean_B = mean(B),
    Bias_B = mean(B - CONFIG$TRUE_B),
    RMSE_B = sqrt(mean((B - CONFIG$TRUE_B)^2)),
    Mean_SE_B = mean(SE_B, na.rm=TRUE),
    SD_Est_B = sd(B),
    
    
    SE_Ratio_B = mean(SE_B, na.rm=TRUE) / sd(B),
    
    Coverage_B = mean(abs(B - CONFIG$TRUE_B) < 1.96 * SE_B, na.rm=TRUE),
    
    .groups = 'drop'
  )

# 列ができたか確認
head(summary_stats)

##可視化----------------------------
# Plot 1: 一致性 (Consistency - RMSE)
p1_b_A_no <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
             aes(x=factor(J), y=Bias_A, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (BIAS Convergence)", subtitle="Condition: No Outlier", x="Items (J)", y="Bias of A")

p1_b_B_no <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
                    aes(x=factor(J), y=Bias_B, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (BIAS Convergence)", subtitle="Condition: No Outlier", x="Items (J)", y="Bias of B")

p1_r_A_no <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
                    aes(x=factor(J), y=RMSE_A, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (RMSE Convergence)", subtitle="Condition: No Outlier", x="Items (J)", y="RMSE of A")

p1_r_B_no <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
                    aes(x=factor(J), y=RMSE_B, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (RMSE Convergence)", subtitle="Condition: No Outlier", x="Items (J)", y="RMSE of B")

p1_b_A <- ggplot(summary_stats %>% filter(Condition == "With_Drift"), 
                 aes(x=factor(J), y=Bias_A, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (BIAS A Convergence)", subtitle="Condition: With_Drift", x="Items (J)", y="BIAS of A")

p1_b_B <- ggplot(summary_stats %>% filter(Condition == "With_Drift"), 
                 aes(x=factor(J), y=Bias_B, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (BIAS B Convergence)", subtitle="Condition: With_Drift", x="Items (J)", y="BIAS of B")

p1_r_A <- ggplot(summary_stats %>% filter(Condition == "With_Drift"), 
                 aes(x=factor(J), y=RMSE_A, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (RMSE A Convergence)", subtitle="Condition: With_Drift", x="Items (J)", y="RMSE of A")

p1_r_B <- ggplot(summary_stats %>% filter(Condition == "With_Drift"), 
                 aes(x=factor(J), y=RMSE_B, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (RMSE B Convergence)", subtitle="Condition: With_Drift", x="Items (J)", y="RMSE of B")

# Plot 2: 分散推定 (Variance - SE Ratio)
p2 <- ggplot(summary_stats %>% filter(Condition == "With_Drift"), 
             aes(x=factor(J), y=SE_Ratio, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + geom_hline(yintercept=1.0, linetype="dashed") +
  facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="2. SE Accuracy (Sandwich/Empirical)", subtitle="Target=1.0 (No Outlier)", x="Items (J)", y="SE Ratio")

# Plot 3: 信頼区間 (Coverage) - 外れ値の影響比較
p3 <- ggplot(summary_stats, aes(x=factor(J), y=Coverage, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + geom_hline(yintercept=0.95, linetype="dashed") +
  facet_grid(M ~ Condition, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="3. CI Coverage Probability", subtitle="Target=0.95. Comparison of Robustness.", x="Items (J)", y="Coverage")

# Plot 4: 漸近正規性 (Normality - QQ Plot)
dat_norm <- all_results %>% filter(Condition=="No_Outlier", J==80, M==301)
p4 <- ggplot(dat_norm, aes(sample=A)) + stat_qq(aes(color=Method)) + stat_qq_line() +
  facet_wrap(~ Method) + theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="4. Asymptotic Normality", subtitle="Condition: No Outlier, J=80, M=81", x="Theoretical", y="Sample")

dat_norm <- all_results %>% filter(Condition=="With_Drift", J==500, M==301)
p4_d <- ggplot(dat_norm, aes(sample=A)) + stat_qq(aes(color=Method)) + stat_qq_line() +
  facet_wrap(~ Method) + theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title=" Asymptotic Normality", subtitle="Condition: With Drift, J=500, M=301", x="Theoretical", y="Sample")

p4_d_B <- ggplot(dat_norm, aes(sample=B)) + stat_qq(aes(color=Method)) + stat_qq_line() +
  facet_wrap(~ Method) + theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title=" Asymptotic Normality", subtitle="Condition: With Drift, J=500, M=301", x="Theoretical", y="Sample")

grid.arrange(p1_b_A, p1_b_B,ncol=2)
grid.arrange(p1_r_A, p1_r_B,ncol=2)
grid.arrange(p1_b_A_no, p1_b_B_no, ncol=2)
grid.arrange(p1_r_A_no, p1_r_B_no, ncol = 2)

library(ggplot2)
library(dplyr)

# 信頼帯プロット------------------------------------------------
# ==============================================================================
# 1. 準備：データ生成（外れ値あり）
# ==============================================================================
set.seed(1234)

# 設定
J_test <- 50             # 項目数
M_test <- 301             # 求積点数
TRUE_A <- 1.2            # 真のA
TRUE_B <- 0.5            # 真のB
DRIFT_MAG <- 2.0         # ドリフト（外れ値）の大きさ
OUTLIER_PROP <- 0.1      # 外れ値の割合 

# パラメータ生成
a_ref <- rlnorm(J_test, 0, 0.3)
b_ref <- rnorm(J_test, 0, 1)

# 真の変換関係を持つパラメータ
a_new_true <- a_ref * TRUE_A
b_new_true <- (b_ref - TRUE_B) / TRUE_A

# 【重要】観測データとして「外れ値（パラメタドリフト）」を混入させる
a_new_obs <- a_new_true
b_new_obs <- b_new_true
idx_drift <- sample(1:J_test, floor(J_test * OUTLIER_PROP))
b_new_obs[idx_drift] <- b_new_obs[idx_drift] + DRIFT_MAG # bパラメータをずらす

# 求積点と重み
nodes <- seq(-4, 4, length.out=M_test)
weights <- dnorm(nodes)
weights <- weights / sum(weights)

# 基準テストの確率分布 P_ref（求積点×項目数の行列）
P_ref <- matrix(0, M_test, J_test)
for(j in 1:J_test) P_ref[,j] <- p_func(nodes, a_ref[j], b_ref[j])


# ==============================================================================
# 2. 推定の実行（ロバスト vs 非ロバスト）
# ==============================================================================
methods_list <- list(
  list(name="DPD(0.01) ", type="DPD", param=0.0001), # KL近似
  list(name="DPD(0.5)",    type="DPD", param=0.5) ,
  list(name="gamma(0.5)",    type="gamma", param=0.5)
  
)


plot_data <- data.frame()
theta_seq <- seq(-4, 4, length.out = 100) # スコアの範囲 (-4から+4まで)

for(met in methods_list) {
  # 推定実行
  res <- estimate_equating_sandwich(
    a_new = a_new_obs, 
    b_new = b_new_obs, 
    P_ref = P_ref, 
    nodes = nodes, 
    weights = weights, 
    div_type = met$type, 
    param = met$param
  )
  
  
  if(res$conv && !any(is.na(res$CovMat))) {
    
    # 共分散行列の要素を取り出し
    Var_A <- res$CovMat[1,1]
    Var_B <- res$CovMat[2,2]
    Cov_AB <- res$CovMat[1,2]
    
    # デルタ法によるSE計算
    # Var(Ax + B) = x^2 * Var(A) + 2*x*Cov(A,B) + Var(B)
    var_y <- (theta_seq^2 * Var_A) + (2 * theta_seq * Cov_AB) + Var_B
    se_y <- sqrt(var_y)
    
    # 変換スコア (点推定)
    y_hat <- res$A * theta_seq + res$B
    
    # データフレームに格納
    tmp_df <- data.frame(
      Method = met$name,
      Theta_New = theta_seq,      # 横軸：変換前のスコア
      Theta_Ref = y_hat,          # 縦軸：変換後のスコア
      Lower = y_hat - 1.96 * se_y, # 95%信頼区間 下限
      Upper = y_hat + 1.96 * se_y  # 95%信頼区間 上限
    )
    plot_data <- bind_rows(plot_data, tmp_df)
  }
}

# ==============================================================================
# 3. 可視化（信頼帯プロット）
# ==============================================================================
g <- ggplot(plot_data, aes(x = Theta_New, y = Theta_Ref)) +
  # 信頼帯 (リボン)
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.3) +
  
  # 推定された変換線 (実線)
  geom_line(aes(color = Method), size = 1.2) +
  
  # 真の変換線 (点線：正解)
  geom_abline(slope = TRUE_A, intercept = TRUE_B, 
              linetype = "dashed", color = "black", size = 0.8) +
  
  # デザイン調整
  facet_wrap(~ Method) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Equating Score Transformation with 95% Confidence Bands",
    subtitle = paste0("True Relation: A=", TRUE_A, ", B=", TRUE_B, 
                      " (With ", OUTLIER_PROP*100, "% Outliers)"),
    x = "New Form Score (theta)",
    y = "Equated Reference Score",
    caption = "Dashed Line: True Relationship / Shaded Area: 95% Confidence Band"
  )

print(g)

# 保存ファイル名
timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
filename <- paste0("workspace_", timestamp, ".RData")

# ------- 現在の作業空間を保存-------------
save.image(file = filename)