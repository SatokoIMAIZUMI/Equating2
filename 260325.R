
rm(list = ls())

  library(ggplot2)
  library(dplyr)
  library(tidyr)


# ==============================================================================
# §0  設定
# ==============================================================================

# ---- 真の等化係数 ----
TRUE_A <- 1.2
TRUE_B <- 0.5

# ---- 尺度因子 ----
D_CONST <- 1.702

# ---- 求積点（N(0,1) を [-4, 4] で離散化）----
M_NODES  <- 301
NODES    <- seq(-4, 4, length.out = M_NODES)
WEIGHTS  <- dnorm(NODES); WEIGHTS <- WEIGHTS / sum(WEIGHTS)

# ---- シミュレーション条件 ----
N_REPS <- 300                          # 繰り返し
J_VEC  <- c(10, 20, 30, 100, 1000)  # 共通項目数


BANK <- list(
  MU_LOG_A = 0.0,    # log(a) の平均  → E[a] = exp(0 + 0.3²/2) ≈ 1.05
  SD_LOG_A  = 0.3,   # log(a) の SD
  MU_B      = 0.0,   # b の平均
  SD_B      = 1.0    # b の SD
)

# ---- 外れ値（パラメータドリフト）設定 ----
CONT_RATE  <- 0.03   # 汚染割合
CONT_DRIFT <- 3.0    # b のドリフト量（汚染項目の b にこれを加算）
CONT_SD    <- 0.5    # ドリフトの SD

# ---- 比較手法 ----
METHODS <- list(
  list(name = "Haebara",      type = "Haebara", param = 0.0),
  list(name = "KL",           type = "KL",      param = 0.0),
  list(name = "DPD(b=0.3)",   type = "DPD",     param = 0.3),
  list(name = "DPD(b=0.5)",   type = "DPD",     param = 0.5),
  list(name = "DPD(b=1.0)",   type = "DPD",     param = 1.0),
  list(name = "g-div(g=0.3)", type = "gamma",   param = 0.3),
  list(name = "g-div(g=0.5)", type = "gamma",   param = 0.5),
  list(name = "g-div(g=1.0)", type = "gamma",   param = 1.0)
)

METHOD_COLORS <- c(
  "Haebara"      = "#E41A1C", # 赤
  "KL"           = "#377EB8", # 青
  "DPD(b=0.3)"   = "#4DAF4A", # 緑
  "DPD(b=0.5)"   = "#984EA3", # 紫
  "DPD(b=1.0)"   = "#1B9E77", # ティール（青緑）
  "g-div(g=0.3)" = "#FF7F00", # オレンジ
  "g-div(g=0.5)" = "#A65628", # 茶色
  "g-div(g=1.0)" = "#F781BF"  # ピンク
)

# ==============================================================================
# §1  IRTモデル関数
# ==============================================================================
p_irt <- function(theta, a, b) {
  1 / (1 + exp(-D_CONST * a * (theta - b)))
}

pi_irt <- function(theta, a, b, A, B) {
  z <- D_CONST * (a / A) * (theta - A * b - B)
  1 / (1 + exp(-z))
}

# ==============================================================================
# §2  ベルヌーイ分布間ダイバージェンス
# ==============================================================================
bern_div <- function(P, pi, type, param) {
  eps <- 1e-10
  P  <- pmax(pmin(P,  1 - eps), eps)
  pi <- pmax(pmin(pi, 1 - eps), eps)
  switch(type,
         DPD = {
           b <- param
           -(1/b)*(P*pi^b + (1-P)*(1-pi)^b) + (1/(1+b))*(pi^(1+b) + (1-pi)^(1+b))
         },
         gamma = {
           g <- param
           -(1/g)*log(P*pi^g + (1-P)*(1-pi)^g) + (1/(1+g))*log(pi^(1+g) + (1-pi)^(1+g))
         },
         KL      = { P*log(P/pi) + (1-P)*log((1-P)/(1-pi)) },
         Haebara = { (P - pi)^2 }
  )
}

# ==============================================================================
# §3  目的関数
# ==============================================================================
obj_fn <- function(eta, a_new, b_new, P_mat, nodes, weights, type, param) {
  A <- eta[1]; B <- eta[2]
  if (abs(A) < 1e-6) return(1e10)
  J <- length(a_new)
  total <- 0
  for (j in 1:J) {
    z_v   <- D_CONST * (a_new[j] / A) * (nodes - A * b_new[j] - B)
    pi_v  <- 1 / (1 + exp(-z_v))
    total <- total + sum(bern_div(P_mat[, j], pi_v, type, param) * weights)
  }
  total
}

# ==============================================================================
# §4  スコア関数 ψ_j とヘッセ寄与 H_j（解析的計算）
# ==============================================================================
score_H_item <- function(a, b, A, B, P_j, nodes, weights, type, param) {
  eps <- 1e-10
  M   <- length(nodes)
  
  z_v   <- D_CONST * (a / A) * (nodes - A * b - B)
  pi_v  <- pmax(pmin(1 / (1 + exp(-z_v)), 1 - eps), eps)
  p1p   <- pi_v * (1 - pi_v)
  P_v   <- pmax(pmin(P_j, 1 - eps), eps)
  
  dz_A  <- D_CONST * a * (B - nodes) / A^2
  dz_B  <- -D_CONST * a / A
  dpi_A <- p1p * dz_A
  dpi_B <- p1p * dz_B
  
  d2z_AA <- 2 * D_CONST * a * (nodes - B) / A^3
  d2z_AB <- D_CONST * a / A^2
  
  d2pi_AA <- p1p*(1-2*pi_v)*dz_A^2       + p1p*d2z_AA
  d2pi_AB <- p1p*(1-2*pi_v)*dz_A*dz_B    + p1p*d2z_AB
  d2pi_BB <- p1p*(1-2*pi_v)*dz_B^2
  
  if (type == "DPD") {
    bv  <- param
    f1  <- pi_v^(bv-1); f2 <- (1-pi_v)^(bv-1)
    dD1 <- (pi_v - P_v) * (f1 + f2)
    dD2 <- (f1 + f2) + (pi_v - P_v)*(bv-1)*(f1/pi_v - f2/(1-pi_v))
  } else if (type == "gamma") {
    gv  <- param
    Ag  <- P_v*pi_v^gv + (1-P_v)*(1-pi_v)^gv
    Bg  <- P_v*pi_v^(gv-1) - (1-P_v)*(1-pi_v)^(gv-1)
    Cg  <- pi_v^(1+gv) + (1-pi_v)^(1+gv)
    Dg  <- pi_v^gv - (1-pi_v)^gv
    dD1 <- -(Bg/Ag) + (Dg/Cg)
    Ag_pr <- gv*Bg
    Bg_pr <- (gv-1)*(P_v*pi_v^(gv-2) + (1-P_v)*(1-pi_v)^(gv-2))
    Cg_pr <- (1+gv)*Dg
    Dg_pr <- gv*(pi_v^(gv-1) + (1-pi_v)^(gv-1))
    dD2   <- -(Bg_pr*Ag - Bg*Ag_pr)/Ag^2 + (Dg_pr*Cg - Dg*Cg_pr)/Cg^2
  } else if (type == "KL") {
    dD1 <- (pi_v - P_v) / p1p
    dD2 <- (p1p - (pi_v - P_v)*(1-2*pi_v)) / p1p^2
  } else {
    dD1 <- 2*(pi_v - P_v); dD2 <- rep(2.0, M)
  }
  
  psi_j <- c(sum(dD1 * dpi_A * weights),
             sum(dD1 * dpi_B * weights))
  
  H_AA <- sum((dD2*dpi_A^2     + dD1*d2pi_AA) * weights)
  H_AB <- sum((dD2*dpi_A*dpi_B + dD1*d2pi_AB) * weights)
  H_BB <- sum((dD2*dpi_B^2     + dD1*d2pi_BB) * weights)
  H_j  <- matrix(c(H_AA, H_AB, H_AB, H_BB), 2, 2)
  
  list(psi = psi_j, H = H_j)
}

# ==============================================================================
# §5  メイン推定関数（点推定 + サンドイッチSE）
# ==============================================================================
equate_fit <- function(a_new, b_new, P_ref, nodes, weights, type, param,
                       start = c(1.0, 0.0)) {
  J <- length(a_new)
  
  fit <- tryCatch(
    optim(par=start, fn=obj_fn,
          a_new=a_new, b_new=b_new, P_mat=P_ref,
          nodes=nodes, weights=weights, type=type, param=param,
          method="BFGS", control=list(maxit=5000, reltol=1e-12)),
    error = function(e) list(convergence=99)
  )
  if (fit$convergence != 0) {
    fit <- tryCatch(
      optim(par=start, fn=obj_fn,
            a_new=a_new, b_new=b_new, P_mat=P_ref,
            nodes=nodes, weights=weights, type=type, param=param,
            method="Nelder-Mead", control=list(maxit=10000)),
      error = function(e) list(convergence=99)
    )
  }
  if (fit$convergence != 0) return(list(conv=FALSE))
  
  A_hat <- fit$par[1]; B_hat <- fit$par[2]
  
  
  H_sum <- matrix(0, 2, 2)
  G_sum <- matrix(0, 2, 2)
  
  for (j in 1:J) {
    out   <- score_H_item(a_new[j], b_new[j], A_hat, B_hat,
                          P_ref[, j], nodes, weights, type, param)
    H_sum <- H_sum + out$H
    G_sum <- G_sum + tcrossprod(out$psi)
  }
  
  H_bar <- H_sum / J
  G_bar <- G_sum / J
  
  inv_H <- tryCatch(solve(H_bar), error = function(e) NULL)
  if (is.null(inv_H)) return(list(conv=FALSE))
  if (kappa(H_bar) > 1e8)
    warning(sprintf("H ill-conditioned (J=%d, kappa=%.2e)", J, kappa(H_bar)))
  
  V_hat <- (inv_H %*% G_bar %*% t(inv_H)) / J
  SE    <- sqrt(pmax(diag(V_hat), 0))
  
  list(conv=TRUE, A=A_hat, B=B_hat, SE_A=SE[1], SE_B=SE[2], V=V_hat)
}

run_once <- function(J, has_outlier, seed) {
  set.seed(seed)
  
  
  a_ref <- rlnorm(J, meanlog = BANK$MU_LOG_A, sdlog = BANK$SD_LOG_A)
  b_ref <- rnorm(J, mean = BANK$MU_B, sd = BANK$SD_B)
  
  
  a_new <- TRUE_A * a_ref
  b_new <- (b_ref - TRUE_B) / TRUE_A
  
  BASE_NOISE_SD <- 0.05
  a_new <- a_new + rnorm(J, mean = 0, sd = BASE_NOISE_SD)
  b_new <- b_new + rnorm(J, mean = 0, sd = BASE_NOISE_SD)
  
  # --- 外れ値: 一部項目の b パラメータをドリフト ---
  if (has_outlier) {
    n_cont <- max(1, round(J * CONT_RATE))
    idx    <- sample(J, n_cont)
    b_new[idx] <- b_new[idx] + rnorm(n_cont, CONT_DRIFT, CONT_SD)
  }
  
  # --- 基準テスト確率行列 P_ref (M × J) ---
  P_ref <- matrix(0, M_NODES, J)
  for (j in 1:J) P_ref[, j] <- p_irt(NODES, a_ref[j], b_ref[j])
  

  results <- list()
  for (met in METHODS) {
    res <- equate_fit(a_new, b_new, P_ref, NODES, WEIGHTS, met$type, met$param)
    if (res$conv) {
      results[[met$name]] <- data.frame(
        A=res$A, B=res$B, SE_A=res$SE_A, SE_B=res$SE_B,
        Method=met$name
      )
    }
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# ==============================================================================
# §7  全シミュレーション実行
# ==============================================================================
cat("=== シミュレーション開始（ランダム項目理論版）===\n")
cat(sprintf("項目バンク: a ~ LogN(%.1f, %.1f²), b ~ N(%.1f, %.1f²)\n",
            BANK$MU_LOG_A, BANK$SD_LOG_A, BANK$MU_B, BANK$SD_B))
cat(sprintf("N_REPS=%d, J_VEC=%s, 手法数=%d\n",
            N_REPS, paste(J_VEC, collapse=","), length(METHODS)))
cat(sprintf("外れ値: %.0f%%の項目でbドリフト(平均+%.1f)\n\n",
            CONT_RATE*100, CONT_DRIFT))

t_start    <- proc.time()
all_results <- list()

for (cond in c("No_Outlier", "With_Outlier")) {
  has_out <- (cond == "With_Outlier")
  cat(sprintf("【%s】\n", cond))
  
  for (J in J_VEC) {
    cat(sprintf("  J=%3d : ", J))
    for (r in seq_len(N_REPS)) {
      seed_val <- r + J * 1000 + has_out * 1e7
      res <- run_once(J, has_out, as.integer(seed_val))
      if (!is.null(res)) {
        res$J <- J; res$Condition <- cond; res$Rep <- r
        all_results[[length(all_results)+1]] <- res
      }
    }
    cat(sprintf("完了 (n=%d)\n", N_REPS))
  }
  cat("\n")
}

all_df <- bind_rows(all_results) %>%
  mutate(Method = factor(Method, levels = sapply(METHODS, `[[`, "name")))

elapsed <- (proc.time() - t_start)["elapsed"]
cat(sprintf("=== 完了 (%.1f秒 = %.1f分) ===\n", elapsed, elapsed/60))

# ==============================================================================
# §8  集計
# ==============================================================================
smry <- all_df %>%
  group_by(Condition, Method, J) %>%
  summarise(
    n_conv     = n(),
    # --- A ---
    Bias_A     = mean(A - TRUE_A),
    SD_A       = sd(A),
    RMSE_A     = sqrt(mean((A - TRUE_A)^2)),
    MeanSE_A   = mean(SE_A, na.rm=TRUE),
    
    SE_Ratio_A = mean(SE_A, na.rm=TRUE) / sd(A),
    Coverage_A = mean(abs(A - TRUE_A) < 1.96 * SE_A, na.rm=TRUE),
    # --- B ---
    Bias_B     = mean(B - TRUE_B),
    SD_B       = sd(B),
    RMSE_B     = sqrt(mean((B - TRUE_B)^2)),
    MeanSE_B   = mean(SE_B, na.rm=TRUE),
    SE_Ratio_B = mean(SE_B, na.rm=TRUE) / sd(B),
    Coverage_B = mean(abs(B - TRUE_B) < 1.96 * SE_B, na.rm=TRUE),
    .groups = "drop"
  )

# ==============================================================================
# §9  可視化
# ==============================================================================
to_long2 <- function(df, col_A, col_B) {
  bind_rows(
    df %>% mutate(Param="A (slope)",     Value=.data[[col_A]]),
    df %>% mutate(Param="B (intercept)", Value=.data[[col_B]])
  )
}

# ---- Figure 0: 項目バンク母集団分布の可視化 ----

set.seed(42)
J_demo <- 1000
demo_items <- data.frame(
  a = rlnorm(J_demo, BANK$MU_LOG_A, BANK$SD_LOG_A),
  b = rnorm(J_demo, BANK$MU_B, BANK$SD_B)
)
fig0_a <- ggplot(demo_items, aes(x=a)) +
  geom_histogram(fill="steelblue", color="white", bins=40) +
  geom_vline(xintercept=exp(BANK$MU_LOG_A), linetype="dashed", color="red") +
  theme_bw(base_size=12) +
  labs(title="Figure 0a: 項目バンク母集団 — 識別力 a",
       subtitle=sprintf("a ~ LogNormal(%.1f, %.1f²), E[a]≈%.2f",
                        BANK$MU_LOG_A, BANK$SD_LOG_A,
                        exp(BANK$MU_LOG_A + BANK$SD_LOG_A^2/2)),
       x="識別力 a", y="度数")

fig0_b <- ggplot(demo_items, aes(x=b)) +
  geom_histogram(fill="coral", color="white", bins=40) +
  geom_vline(xintercept=BANK$MU_B, linetype="dashed", color="red") +
  theme_bw(base_size=12) +
  labs(title="Figure 0b: 項目バンク母集団 — 困難度 b",
       subtitle=sprintf("b ~ Normal(%.1f, %.1f²)", BANK$MU_B, BANK$SD_B),
       x="困難度 b", y="度数")

# ---- Figure 1: 一致性 — Bias ----
fig1 <- smry %>% to_long2("Bias_A", "Bias_B") %>%
  ggplot(aes(J, Value, color=Method, group=Method)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray50") +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  facet_grid(Param ~ Condition, scales="free_y") +
  scale_x_log10(breaks=J_VEC) +
  scale_color_manual(values=METHOD_COLORS) + theme_bw(base_size=12) +
  labs(title="Figure 1: 一致性 — Bias",
       subtitle="J↑ → Bias→0（一致性の確認）",
       x="共通項目数 J（対数スケール）", y="Bias", color="Method")

# ---- Figure 2: 一致性 — RMSE（両軸対数）----
fig2 <- smry %>% to_long2("RMSE_A", "RMSE_B") %>%
  ggplot(aes(J, Value, color=Method, group=Method)) +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  facet_grid(Param ~ Condition, scales="free_y") +
  scale_x_log10(breaks=J_VEC) + scale_y_log10() +
  scale_color_manual(values=METHOD_COLORS) + theme_bw(base_size=12) +
  labs(title="Figure 2: 一致性 — RMSE（両軸対数）",
       
       x="共通項目数 J（対数スケール）", y="RMSE（対数スケール）")

# ---- Figure 3: SE比（ランダム項目理論の検証）----
fig3 <- smry %>% to_long2("SE_Ratio_A", "SE_Ratio_B") %>%
  ggplot(aes(J, Value, color=Method, group=Method)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=1.1, color="gray30") +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  facet_grid(Param ~ Condition) +
  coord_cartesian(ylim=c(0, 2.5)) +
  scale_x_log10(breaks=J_VEC) +
  scale_color_manual(values=METHOD_COLORS) + theme_bw(base_size=12) +
  labs(title="Figure 3: SE比 = E[SE_hat] / SD(η_hat)（サンドイッチ推定の検証）",
       subtitle="目標値 = 1.0（破線）",
       x="共通項目数 J（対数スケール）",
       y="SE比（サンドイッチSE / 経験的SD）")

# ---- Figure 4: 被覆確率 ----
fig4 <- smry %>% to_long2("Coverage_A", "Coverage_B") %>%
  ggplot(aes(J, Value, color=Method, group=Method)) +
  geom_hline(yintercept=0.95, linetype="dashed", linewidth=1.1, color="gray30") +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  facet_grid(Param ~ Condition) +
  coord_cartesian(ylim=c(0.6, 1.02)) +
  scale_x_log10(breaks=J_VEC) +
  scale_color_manual(values=METHOD_COLORS) + theme_bw(base_size=12) +
  labs(title="Figure 4: 95%信頼区間の被覆確率",
       subtitle="目標値 = 0.95（破線）",
       x="共通項目数 J（対数スケール）", y="Coverage")

# ---- Figure 5: 漸近正規性 QQプロット（J最大, No_Outlier）----
J_large <- max(J_VEC)
fig5 <- all_df %>%
  filter(Condition=="No_Outlier", J==J_large) %>%
  group_by(Method) %>%
  mutate(Z_A=scale(A)[,1], Z_B=scale(B)[,1]) %>% ungroup() %>%
  pivot_longer(c(Z_A, Z_B), names_to="Param", values_to="Z") %>%
  mutate(Param=recode(Param, Z_A="A (slope)", Z_B="B (intercept)")) %>%
  ggplot(aes(sample=Z)) +
  stat_qq(aes(color=Method), size=0.8, alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.8) +
  facet_grid(Param ~ Method) +
  scale_color_manual(values=METHOD_COLORS, guide="none") +
  theme_bw(base_size=11) +
  labs(title=sprintf("Figure 5: 漸近正規性 QQプロット（J=%d, 外れ値なし）", J_large),
       subtitle="標準化推定量が対角線上に乗れば漸近正規性が成立",
       x="理論分位数 N(0,1)", y="標本分位数（標準化）")


# ---- 全図表示 ----
print(fig0_a); print(fig0_b)
print(fig1); print(fig2); print(fig3); print(fig4); print(fig5)

# ==============================================================================
# §10  数値テーブル（論文掲載用）
# ==============================================================================
fmt <- function(x) round(x, 4)

cat("\n")
cat(rep("=",70), "\n", sep="")
cat("Table 1: 統計的性質（外れ値なし条件）\n")
cat("  SE比 → 1.0 に収束: ランダム項目理論によりサンドイッチSEが正確であることの証左\n")
cat(rep("=",70), "\n", sep="")
smry %>%
  filter(Condition=="No_Outlier") %>%
  select(Method, J, Bias_A, RMSE_A, SE_Ratio_A, Coverage_A,
         Bias_B, RMSE_B, SE_Ratio_B, Coverage_B) %>%
  mutate(across(where(is.numeric), fmt)) %>%
  arrange(J, Method) %>% print(n=Inf)

cat("\n")
cat(rep("=",70), "\n", sep="")
cat("Table 2: ロバスト性（外れ値あり条件）\n")
cat(rep("=",70), "\n", sep="")
smry %>%
  filter(Condition=="With_Outlier") %>%
  select(Method, J, Bias_A, RMSE_A, SE_Ratio_A, Coverage_A,
         Bias_B, RMSE_B, SE_Ratio_B, Coverage_B) %>%
  mutate(across(where(is.numeric), fmt)) %>%
  arrange(J, Method) %>% print(n=Inf)

# ==============================================================================
# §11  保存
# ==============================================================================
save.image(file=paste0("equating_random_item_", format(Sys.time(),"%Y%m%d_%H%M"), ".RData"))
cat("\n作業空間を保存しました。\n")