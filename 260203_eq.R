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

set.seed(123)

# ★論文用設定★
CONFIG <- list(
  N_EXAMINEES = 1000,
  N_REPS = 20,
  N_UNIQUE = 20,
  TRUE_A = 1.2,
  TRUE_B = 0.5,
  D_CONST = 1.702,
  DRIFT_MAG = 2.5,       # ドリフト量
  OUTLIER_PROP = 0.3     # 外れ値の割合 
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
  eps <- 1e-9
  P <- pmax(pmin(P, 1-eps), eps)
  pi_val <- pmax(pmin(pi_val, 1-eps), eps)
  
  if(div_type == "DPD") { 
    beta <- param
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
    pi_vec <- pi_func(nodes, a_new[j], b_new[j], A, B, D)
    div_vec <- calc_div_element(P_ref[,j], pi_vec, div_type, param)
    loss <- loss + sum(div_vec * weights)
  }
  return(loss)
}


# Estimation with Sandwich Variance
estimate_equating_sandwich <- function(a_new, b_new, P_ref, nodes, weights, div_type, param, D=1.702){
  start_val <- c(CONFIG$TRUE_A, CONFIG$TRUE_B)#初期値
  
  # 1. Point Estimation
  opt <- tryCatch({
    optim(par = start_val, fn = obj_func_total, 
          a_new = a_new, b_new = b_new, P_ref = P_ref, 
          nodes = nodes, weights = weights, div_type = div_type, param = param, D = D,
          method = "Nelder-Mead", control = list(maxit = 5000))
  }, error = function(e) NULL)
  
  if(is.null(opt) || opt$convergence != 0) return(list(conv=FALSE))
  
  hat_A <- opt$par[1]; hat_B <- opt$par[2]
  
  # 2. Sandwich Variance: V = H^-1 G H^-1
  #へシアン行列(目的関数の2回微分）2×2
  H <- numDeriv::hessian(obj_func_total, x=c(hat_A, hat_B), 
                         a_new=a_new, b_new=b_new, P_ref=P_ref, 
                         nodes=nodes, weights=weights, div_type=div_type, param=param, D=D)
  
  
  obj_func_single <- function(eta, j, idx_m, a_new, b_new, P_ref, nodes, weights, div_type, param, D) {
    # 1つの項目j、1つの求積点mだけの損失
    A <- eta[1]; B <- eta[2]
    # idx_m を使ってアクセス
    pi_val <- pi_func(nodes[idx_m], a_new[j], b_new[j], A, B, D)
    div_val <- calc_div_element(P_ref[idx_m, j], pi_val, div_type, param)
    return(div_val * weights[idx_m])
  }
  
  J <- length(a_new)
  M <- length(nodes)
  G <- matrix(0, 2, 2)
  
  for(j in 1:J) {
    for(m in 1:M) {
      
      psi_jm <- numDeriv::grad(obj_func_single, x=c(hat_A, hat_B),
                               j=j, idx_m=m, 
                               a_new=a_new, b_new=b_new, P_ref=P_ref,
                               nodes=nodes, weights=weights, div_type=div_type, param=param, D=D)
      
      G <- G + (psi_jm %*% t(psi_jm))#勾配の外積和
    }
  }
  
  inv_H <- try(solve(H), silent=TRUE)#Hの逆行列
  if(inherits(inv_H, "try-error")) {
    CovMat <- matrix(NA, 2, 2)
  } else {
    CovMat <- inv_H %*% G %*% inv_H#共分散行列
  }
  
  SE <- sqrt(diag(CovMat))#標準誤差の計算
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
    
    if(has_outlier){#外れ値の導入
      n_drift <- floor(J * CONFIG$OUTLIER_PROP)
      if(n_drift > 0){
        idx <- sample(1:J, n_drift)
        b_com_new_true[idx] <- b_com_new_true[idx] + CONFIG$DRIFT_MAG
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
  
  # IRT Calibration
  mirt_res <- try({
    invisible(capture.output({
      mod_ref <- mirt(as.data.frame(resp_ref), 1, itemtype="2PL", verbose=FALSE, technical=list(NCYCLES=500))
      mod_new <- mirt(as.data.frame(resp_new), 1, itemtype="2PL", verbose=FALSE, technical=list(NCYCLES=500))
    }))
    list(ref=mod_ref, new=mod_new)
  }, silent=TRUE)
  
  if(inherits(mirt_res, "try-error")) return(NULL)
  #共通項目のパラメータ抽出
  cfs_ref <- coef(mirt_res$ref, IRTpars=TRUE, simplify=TRUE)$items
  cfs_new <- coef(mirt_res$new, IRTpars=TRUE, simplify=TRUE)$items
  
  # Extract Common Items
  est_a_ref_com <- cfs_ref[1:J, "a"] / CONFIG$D_CONST; est_b_ref_com <- cfs_ref[1:J, "b"]
  est_a_new_com <- cfs_new[1:J, "a"] / CONFIG$D_CONST; est_b_new_com <- cfs_new[1:J, "b"]
  
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

# ------------------------------------------------------------------------------
# 4. 実験メインループ (Monte Carlo Simulation)
# ------------------------------------------------------------------------------
results_list <- list()
counter <- 1

J_vec <- c(5, 20, 80)#共通項目数
M_vec <- c(21, 41, 101)# 求積点数
conditions <- expand.grid(J=J_vec, M=M_vec)#全組み合わせ

methods <- list(
  list(name="KL (Approx)", type="DPD", p=0.01),
  list(name="DPD (0.3)",   type="DPD", p=0.3),
  list(name="DPD (0.5)",   type="DPD", p=0.5),
  list(name="DPD (1.0)",   type="DPD", p=1.0),
  list(name="Gamma (0.3)", type="gamma", p=0.3),
  list(name="Gamma (0.5)", type="gamma", p=0.5)
)

outlier_conds <- c(FALSE, TRUE)

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
    .groups = 'drop'
  )

# 列ができたか確認
head(summary_stats)
# Plot 1: 一致性 (Consistency - RMSE)
p1 <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
             aes(x=factor(J), y=RMSE_A, group=Method, color=Method)) +
  geom_line(size=1) + geom_point(size=3) + facet_wrap(~ M, labeller = label_both) +
  theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="1. Consistency (RMSE Convergence)", subtitle="Condition: No Outlier", x="Items (J)", y="RMSE of A")

# Plot 2: 分散推定 (Variance - SE Ratio)
p2 <- ggplot(summary_stats %>% filter(Condition == "No_Outlier"), 
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
dat_norm <- all_results %>% filter(Condition=="No_Outlier", J==80, M==81)
p4 <- ggplot(dat_norm, aes(sample=A)) + stat_qq(aes(color=Method)) + stat_qq_line() +
  facet_wrap(~ Method) + theme_bw() + scale_color_brewer(palette="Set1") +
  labs(title="4. Asymptotic Normality", subtitle="Condition: No Outlier, J=80, M=81", x="Theoretical", y="Sample")

grid.arrange(p1, p2, p3, p4, ncol=2)


library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. 準備：データ生成（外れ値あり）
# ==============================================================================
set.seed(123)

# 設定
J_test <- 40             # 項目数
M_test <- 41             # 求積点数
TRUE_A <- 1.2            # 真のA
TRUE_B <- 0.5            # 真のB
DRIFT_MAG <- 3.0         # ドリフト（外れ値）の大きさ
OUTLIER_PROP <- 0.3      # 外れ値の割合 

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

# 基準テストの確率分布 P_ref
P_ref <- matrix(0, M_test, J_test)
for(j in 1:J_test) P_ref[,j] <- p_func(nodes, a_ref[j], b_ref[j])


# ==============================================================================
# 2. 推定の実行（ロバスト vs 非ロバスト）
# ==============================================================================
methods_list <- list(
  list(name="KL ", type="DPD", param=0.01), # KL近似
  list(name="Robust (DPD)",    type="DPD", param=0.3)   # 提案手法
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