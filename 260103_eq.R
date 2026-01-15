setwd("/Users/imaizumisatoko/Documents/研究/AMPC2026")

rm(list = ls())
set.seed(1234)
library(tidyverse)
library(mirt)
library(ggmirt)
library(Rmpfr)
library(openxlsx)
library(xtable)
library(readr)
library(dplyr)

#2plm
icc2PL <- function(a, b, th, clip_eps = 1e-6) {
  p <- outer(th, seq_along(a), function(t, i) plogis(a[i] * (t - b[i])))#すべての項目，能力の組み合わせで計算(plogisはロジスティック関数計算)
  # クリップ処理で0や1を避ける
  p <- pmax(pmin(p, 1 - clip_eps), clip_eps)#能力×項目の行列
  return(p)
}



# ------------------------------------------------------------------------------
# 1. 補助関数 (Utility Functions)
# ------------------------------------------------------------------------------

# 安全なクリッピング (0,1区間)
clip_prob <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

# ベルヌーイ分布のPE距離 (Probabilistic Error)
PEB <- function(p1, p2) {
  p2 <- clip_prob(p2)
  (p2 - p1)^2 / (p2 * (1 - p2))
}

# ベルヌーイ分布のKLダイバージェンス
BKL <- function(p1, p2) {
  p1 <- clip_prob(p1)
  p2 <- clip_prob(p2)
  p1 * log(p1 / p2) + (1 - p1) * log((1 - p1) / (1 - p2))
}

# ベルヌーイ分布のJSダイバージェンス
JSB <- function(p1, p2) {
  m <- 0.5 * (p1 + p2)
  0.5 * (BKL(p1, m) + BKL(p2, m))
}

# Hellinger距離
HEB <- function(p1, p2) {
  sqrt(1 - (sqrt(p1 * p2) + sqrt((1 - p1) * (1 - p2))))
}

# ------------------------------------------------------------------------------
# 2. 最適化コア関数 (Core Optimization Logic)
# ------------------------------------------------------------------------------
# 全ての手法で共通する「初期値グリッド探索 -> optim -> 結果整形」を担う関数

run_equating_optim <- function(obj_func, a1, b1, a2, b2, 
                               A_range, B_range, trA, trB) {
  
  # 初期値グリッドの作成
  init_grid <- expand.grid(fA = A_range, fB = B_range)
  results <- list()
  
  for (i in seq_len(nrow(init_grid))) {
    fA <- init_grid$fA[i]
    fB <- init_grid$fB[i]
    
    # optim実行 (エラーハンドリング付き)
    res_try <- try(
      optim(par = c(fA, fB), fn = obj_func, method = "BFGS", control = list(reltol = 1e-8)),
      silent = TRUE
    )
    
    # 失敗時はスキップ
    if (inherits(res_try, "try-error") || is.null(res_try$par) || any(is.na(res_try$par))) {
      next
    }
    
    # パラメータ変換 (最適化空間から元の空間へ)
    # ※ 元コードの定義に従い A = 1/x[1], B = -x[2]/x[1] とする
    est_A <- 1 / res_try$par[1]
    est_B <- -res_try$par[2] / res_try$par[1]
    
    # 結果格納
    results[[length(results) + 1]] <- data.frame(
      fA_init  = fA,
      fB_init  = fB,
      A        = est_A,
      B        = est_B,
      dA       = (est_A - trA)^2,
      dB       = (est_B - trB)^2,
      dth      =  mean((est_A * est.th2[,1] + est_B - th1) ** 2),
      cf_final = res_try$value
    )
  }
  
  # 有効な結果がない場合
  if (length(results) == 0) {
    warning("All initial values failed to converge.")
    return(NULL)
  }
  
  results_df <- do.call(rbind, results)
  
  # 最小損失の結果を返す
  best_row <- which.min(results_df$cf_final)
  return(results_df[best_row, ])
}

# ------------------------------------------------------------------------------
# 3. モーメント法 (Moment Methods)
# ------------------------------------------------------------------------------

# MS (Mean-Sigma)
MS <- function(b1, b2) {
  A <- sd(b2) / sd(b1)
  B <- mean(b2) - A * mean(b1)
  dth <- mean((A * est.th2[,1] + B - th1)^2)
  return(list(A = A, B = B, dA = (A - trA)^2, dB = (B - trB)^2, dth = dth))
}

# MM (Mean-Mean)
MM <- function(a1, b1, a2, b2) {
  A <- mean(a1) / mean(a2)
  B <- mean(b2) - A * mean(b1)
  dth <- mean((A * est.th2[,1] + B - th1)^2)
  return(list(A = A, B = B, dA = (A - trA)^2, dB = (B - trB)^2,dth = dth))
}

# ------------------------------------------------------------------------------
# 4. 特性曲線手法 (IRF Methods)
# ------------------------------------------------------------------------------
# 各関数は「目的関数(cf)」を定義し，run_equating_optim を呼ぶだけにする

# Haebara
HA <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    # Haebara: 確率の差の二乗和
    Q1 <- sum(rowSums((tra1$p1 - tra1$p1star)^2) * H1)
    Q2 <- sum(rowSums((tra1$p2 - tra1$p2star)^2) * H2)
    return(Q1 + Q2)
  }
  
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# SL (Stocking-Lord)
SL <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    # SL: テスト特性曲線(TCC)の差の二乗和
    # transの戻り値の構造に依存 (リストのインデックスを使用していた部分)
    # tra1[[1]]=p1, [[2]]=p1star, [[3]]=p2, [[4]]=p2star と仮定
    Q1 <- sum((rowSums(tra1[[1]]) - rowSums(tra1[[2]]))^2 * H1)
    Q2 <- sum((rowSums(tra1[[3]]) - rowSums(tra1[[4]]))^2 * H2)
    return(Q1 + Q2)
  }
  
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# RA (Density Ratio / Relative Accuracy?)
RA <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums((1 - tra1[[2]] / tra1[[1]])^2) * H1)
    Q2 <- sum(rowSums((1 - tra1[[4]]/ tra1[[3]])^2) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# IC2 (Mean Squared Log Error)
IC2 <- function(a1, b1, a2, b2, ym, H1, H2, 
                A_range = seq(0.5, 1.5, by = 0.25),
                B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums((log1p(tra1[[1]]) - log1p(tra1[[2]]))^2) * H1)
    Q2 <- sum(rowSums((log1p(tra1[[3]]) - log1p(tra1[[4]]))^2) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# PE (Posterior/Probabilistic Error?)
PE <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(PEB(tra1[[1]], tra1[[2]])) * H1)
    Q2 <- sum(rowSums(PEB(tra1[[3]], tra1[[4]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# PER (Reverse PE)
PER <- function(a1, b1, a2, b2, ym, H1, H2,
                A_range = seq(0.5, 1.5, by = 0.25),
                B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(PEB(tra1[[2]], tra1[[1]])) * H1)
    Q2 <- sum(rowSums(PEB(tra1[[4]], tra1[[3]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# PES (Symmetric PE)
PES <- function(a1, b1, a2, b2, ym, H1, H2, 
                A_range = seq(0.5, 1.5, by = 0.25),
                B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    # PEB(P, Q) + PEB(Q, P)
    Q1 <- sum(rowSums(PEB(tra1[[2]], tra1[[1]]) + PEB(tra1[[1]], tra1[[2]])) * H1)
    Q2 <- sum(rowSums(PEB(tra1[[4]], tra1[[3]]) + PEB(tra1[[3]], tra1[[4]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# KL (Kullback-Leibler)
KL <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(BKL(tra1[[1]], tra1[[2]])) * H1)
    Q2 <- sum(rowSums(BKL(tra1[[3]], tra1[[4]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# KLR (Reverse KL)
KLR <- function(a1, b1, a2, b2, ym, H1, H2,
                A_range = seq(0.5, 1.5, by = 0.25),
                B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(BKL(tra1[[2]], tra1[[1]])) * H1)
    Q2 <- sum(rowSums(BKL(tra1[[4]], tra1[[3]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# JS (Jensen-Shannon)
JS <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(JSB(tra1[[1]], tra1[[2]])) * H1)
    Q2 <- sum(rowSums(JSB(tra1[[3]], tra1[[4]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# HE (Hellinger)
HE <- function(a1, b1, a2, b2, ym, H1, H2, 
               A_range = seq(0.5, 1.5, by = 0.25),
               B_range = seq(-1, 1, by = 0.25)) {
  
  cf <- function(x) {
    A <- x[1]; B <- x[2]
    tra1 <- trans(a1 = a1, b1 = b1, a2 = a2, b2 = b2, A = A, B = B)
    Q1 <- sum(rowSums(HEB(tra1[[1]], tra1[[2]])) * H1)
    Q2 <- sum(rowSums(HEB(tra1[[3]], tra1[[4]])) * H2)
    return(Q1 + Q2)
  }
  run_equating_optim(cf, a1, b1, a2, b2, A_range, B_range, trA, trB)
}

# ------------------------------------------------------------------------------
# 5. 同時推定 (Concurrent Calibration)
# ------------------------------------------------------------------------------

CC <- function(r1_cc, r2_cc,  group1_name = "G1", group2_name = "G2") {
  require(mirt)
  
  # データ整形
  n_items <- ncol(r1_cc)
  # 列名を統一
  colnames(r1_cc) <- paste0("I", 1:n_items)
  colnames(r2_cc) <- paste0("I", 1:n_items)
  
  comb_data <- rbind(r1_cc, r2_cc)
  group_id <- factor(c(rep(group1_name, nrow(r1_cc)), rep(group2_name, nrow(r2_cc))))
  
  # multipleGroup 推定
  fit <- tryCatch({
    multipleGroup(
      data = comb_data,
      model = 1,
      itemtype = "2PL",
      group = group_id,
      invariance = c(paste0("I", 1:n_items), "free_means", "free_vars")
    )
  }, error = function(e) {
    warning("mirt estimation failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  # パラメータ取得
  params <- coef(fit, simplify = TRUE)
  
  # G1の平均・分散 (G1が参照群としてN(0,1)固定の場合、params$G1$means等は0,1になるはず)
  # もしG1も推定している場合はここから取得
  mu1 <- params[[group1_name]]$means
  sig1 <- sqrt(params[[group1_name]]$cov)
  
  mu2 <- params[[group2_name]]$means
  sig2 <- sqrt(params[[group2_name]]$cov)
  
  # 等化係数の計算 
  A <- as.numeric(sig1 / sig2)
  B <- as.numeric(mu1 - A * mu2)
  dth <- mean((A * est.th2[,1] + B - th1)^2)
  return(list(A = A, B = B, dA = (A - trA)^2, dB = (B - trB)^2,dth=dth))
}


# Function to generate response patterns
genPattern = function(a, b, th){
  I <- length(th)#受験者数
  J <- length(a)#項目数
  # 1. すべての受験者と項目のペアに対する正答確率を行列で一気に計算
  prob_matrix <- icc2PL(a, b, th)
  # 2. 一様乱数の行列を生成
  rand_matrix <- matrix(runif(I * J), nrow = I, ncol = J)
  # 3. 確率行列と乱数行列を比較して、0/1の応答行列を一気に生成
  r <- ifelse(prob_matrix > rand_matrix, 1, 0)
  return(r)
}

#全員正解（不正解）の項目を取り除く
rem_all01 <- function(resp) {
  resp <- as.matrix(resp)#でーたを行列に
  
  # 列名がない場合は整数IDを設定
  if (is.null(colnames(resp))) {
    colnames(resp) <- as.character(seq_len(ncol(resp)))#1列目から順にラベリング
  }
  
  # 全員正解 or 不正解の列を検出
  remove_idx <- apply(resp, 2, function(col) all(col == 0) || all(col == 1))
  
  kept_items <- colnames(resp)[!remove_idx]
  removed_items <- colnames(resp)[remove_idx]
  
  resp_clean <- resp[, kept_items, drop = FALSE]
  
  # 数値名に変換できる場合だけ数値化
  safe_as_int <- function(x) {
    out <- suppressWarnings(as.integer(x))
    if (any(is.na(out))) {
      # "Item3" → 3 などに変換
      out <- as.integer(gsub("[^0-9]", "", x))
    }
    return(out)
  }
  
  return(list(
    resp = resp_clean,
    removed_items = safe_as_int(removed_items),
    kept_items = safe_as_int(kept_items)
  ))
}



#項目パラ,能力パラの推定
est2plm <- function(x) {#反応行列が引数
  
  # mirtモデルを推定
  fit2PL <- mirt(as.data.frame(x), 1, itemtype = "2PL", verbose = FALSE)#1次元性，2PLM
  
  # 項目パラメタを抽出して変数に保存
  item_params <- coef(fit2PL, IRTpars = TRUE, simplify = TRUE)$items#項目数×パラメータ数
  
  
  # 能力パラメタを抽出して変数に保存
  theta_scores <- fscores(fit2PL, method = "EAP")
  
  # 両方の結果を名前付きリストにまとめて返す
  return(list(items = item_params, scores = theta_scores))
}





#共通項目選定(ブロック法）
select_anchor <- function(item_param, n_ci){
  n_items <- nrow(item_param)
  
  item_sort <- item_param %>% 
    mutate(origID = ItemID) %>% 
    arrange(desc(b))
  
  block_size <- n_items / n_ci
  
  item_block <- item_sort %>% 
    mutate(block = gl(n = n_ci, k = block_size, length = n_items))
  
  select_anchor <- item_block %>% 
    group_by(block) %>% 
    slice_max(order_by = a, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
  return(select_anchor)
}

# base: group1
trans=function(a1, b1, a2, b2, A, B){
  # パラメータ変換
  a2star <- a2 / A
  b2star <- A * b2 + B
  #base; group2
  a1star <- a1 * A
  b1star <- b1 / A - B / A
  
  # ベクトル化された確率
  p1 <- icc2PL(a1,b1,ym)　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　                                                                                                                                                                                                                                      　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
  p1star <- icc2PL(a2star, b2star, ym)
  p2 <- icc2PL(a2,b2,ym)
  p2star <- icc2PL(a1star, b1star, ym)
  
  return(list(p1=p1, p1star=p1star, p2=p2, p2star=p2star))
}




# - シミュレーション設定------

nsim <- 10 #シミュレーション回数


#set of theta param
# 能力値の母平均と母標準偏差
mu1 <-  0 # 0/ -1/ 0/ 0/1 /-1
sig1 <-  1 # 1/ 1/ 1/ 1/ 2/ 1
mu2 <-  0.5 # 0.5/ 1/ 0/ 4.5/ 1.7/ -0.7
sig2 <-  1.2# 1.2/ 1/ 3/ 2.5/ 2.4/ 1.2




#総項目数
n.test1 <- 50 #項目数
n.test2 <- 50 #項目数

#共通項目
m_values <- c(25) # common item 10,25

# 受験者数
i_values <- c(3000) 



#結果格納ベクトル
result_list <- list()
init_list <- list()

results <- list()

remidx <- list()
remidx_lists <- list()

results_modi <- list()

#等化のシミュレーション
for(i in i_values){#受験者数
  for(m in m_values){#共通項目
    anchor.items <- 1:m
    resMSs <- data.frame()
    resMMs <- data.frame()
    resHAs <- data.frame()
    resSLs <- data.frame()
    resRAs <- data.frame()
    
    resIC2s <- data.frame()
    resPEs <- data.frame()
    resPERs <- data.frame()
    resPESs <- data.frame()
    resKLs <- data.frame()
    resKLRs <- data.frame()
    resJSs <- data.frame()
    resHEs <- data.frame()
    resCCs <- data.frame()
    
    
    
    
    resMSs_list <- vector("list", nsim) 
    resMMs_list <- vector("list", nsim) 
    resHAs_list <- vector("list", nsim) 
    resSLs_list <- vector("list", nsim) 
    resRAs_list <- vector("list", nsim) 
    #resICs_list <- vector("list", nsim) 
    resIC2s_list <- vector("list", nsim) 
    resPEs_list <- vector("list", nsim) 
    resPERs_list <- vector("list", nsim) 
    resPESs_list <- vector("list", nsim) 
    resKLs_list <- vector("list", nsim) 
    resKLRs_list <- vector("list", nsim) 
    resJSs_list <- vector("list", nsim) 
    resHEs_list <- vector("list", nsim) 
    resCCs_list <- vector("list", nsim) 
    
    for(s in 1:nsim){#繰り返し数（テスト実施回数）
      
      #test1
      #項目パラメタの真値
      #正規分布の平均と標準偏差
      norm_mu <- log(0.5)
      norm_sd <- 0.5
      
      #対数正規分布の期待値と分散
      log_mu <- exp(norm_mu+norm_sd^2/2)
      log_var <- (exp(norm_sd^2)-1)*exp(2*norm_mu+norm_sd^2)
      
      
      
      true_params1 <- data.frame(
        ItemID = 1:n.test1,
        a1 = rlnorm(n.test1, norm_mu, norm_sd)  ,#識別力
        b1 = rnorm(n.test1, 0, 1)#困難度普通
      )
      
      #ability parameter
      th1 <- rnorm(i, mu1, sig1)#真値
      
      
      # --- (1) Form1反応データ生成 ---
      r1 <- genPattern(true_params1$a1, true_params1$b1, th1)
      colnames(r1) <- as.character(true_params1$ItemID)
      
      # 全員正解/不正解項目を削除
      tmp1 <- rem_all01(r1)
      r1 <- tmp1$resp
      kept_idx1 <- tmp1$kept_items
      rem_idx1 <- tmp1$removed_items
      
      # --- (2) Form1パラメータ推定 ---
      est1_raw <- est2plm(r1)
      estIt1_raw <- est1_raw$items
      estIt1 <- data.frame(
        ItemID = kept_idx1,
        a = estIt1_raw[, "a"],
        b = estIt1_raw[, "b"]
        
      )
      
      # --- (3) 共通項目選定 ---
      anchor_item_selection <- select_anchor(estIt1, m)
      anchor_IDs <- anchor_item_selection$ItemID
      
      
      # --- (4) Form2反応データ生成 ---
      #グループ2の真の能力分布
      th2 <- rnorm(i, mu2, sig2)
      
      # 共通項目の真値をForm1から引き継ぐ
      true_anchor_params <- true_params1[anchor_IDs, ]
      
      # Form2の独自項目を生成
      a2_unique <- rlnorm(n.test2 - m, norm_mu, norm_sd)
      b2_unique <- rnorm(n.test2 - m, 0, 1)
      
      # 真値の構成（共通 + 独自）
      a2_true <- c(true_anchor_params$a1, a2_unique)
      b2_true <- c(true_anchor_params$b1, b2_unique)
      
      # Form2反応データ生成
      r2 <- genPattern(a2_true, b2_true, th1)
      
      # 列名を設定（共通＋独自）
      colnames(r2) <- c(
        as.character(anchor_IDs),
        as.character((max(anchor_IDs) + 1):(max(anchor_IDs) + n.test2 - m))
      )
      
      # --- (5) 全員正解/不正解項目を削除 ---
      tmp2 <- rem_all01(r2)
      r2 <- tmp2$resp
      kept_idx2 <- tmp2$kept_items
      rem_idx2 <- tmp2$removed_items
      
      # --- (6) 共通項目の整合確認 ---
      # 両フォームで残った共通項目だけを抽出（文字列ベースで厳密に）
      valid_anchor_IDs <- intersect(
        as.character(anchor_IDs),
        intersect(colnames(r1), colnames(r2))
      )
      valid_anchor_IDs <- as.integer(valid_anchor_IDs)
      
      # 共通項目が削除されていた場合の処理
      if (length(valid_anchor_IDs) < length(anchor_IDs)) {
        cat(sprintf("共通項目が %d/%d 残存\n",
                    length(valid_anchor_IDs), length(anchor_IDs)))
        
        # Form1・Form2 両方から削除された共通項目を除外
        r1 <- r1[, intersect(colnames(r1), as.character(valid_anchor_IDs)), drop = FALSE]
        r2 <- r2[, intersect(colnames(r2), as.character(valid_anchor_IDs)), drop = FALSE]
      }
      
      # --- 安全チェックを追加 ---
      if (length(valid_anchor_IDs) == 0 ||
          is.null(r1) || is.null(r2) ||
          ncol(r1) == 0 || ncol(r2) == 0) {
        cat("⚠️ 共通項目がすべて削除されました。このシミュレーションをスキップします。\n")
        next  # ← forループ内で使う場合
      }
      
      
      # --- (7) Form2パラメータ推定 ---
      est2_raw <- est2plm(r2)
      estIt2_raw <- est2_raw$items
      
      # IDと推定結果を対応づけ
      estIt2 <- data.frame(
        ItemID = as.numeric(rownames(estIt2_raw)),
        a = estIt2_raw[, "a"],
        b = estIt2_raw[, "b"]
        
      )
      
      # --- 共通項目部分の推定値抽出 ---
      est.a1 <- estIt1 %>% filter(ItemID %in% valid_anchor_IDs) %>% pull(a)
      est.b1 <- estIt1 %>% filter(ItemID %in% valid_anchor_IDs) %>% pull(b)
      est.a2 <- estIt2 %>% filter(ItemID %in% valid_anchor_IDs) %>% pull(a)
      est.b2 <- estIt2 %>% filter(ItemID %in% valid_anchor_IDs) %>% pull(b)
      
      # --- 共通項目の反応データ ---
      r1_cc <- r1[, as.character(valid_anchor_IDs), drop = FALSE]
      r2_cc <- r2[, as.character(valid_anchor_IDs), drop = FALSE]
      
      
      
      # --- 能力推定値それぞれで0,1尺度 ---
      est.th1 <- est1_raw$scores
      est.th2 <- est2_raw$scores
      　
      
      
      
      
      # --- 状況確認ログ ---
      cat(sprintf(
        "m=%g, i=%g, sim=%g : valid_anchor_IDs=%d, estIt2_raw=%d\n",
        m, i, s, length(valid_anchor_IDs), nrow(estIt2_raw)
      ))
      
      
    
      # 求積点と重み
      ym <- seq(-4, 4, length.out = 50)
      w <- dnorm(ym, 0, 1)
      w <- w / sum(w)
      H1 <- w
      H2 <- w
      
      
      p1 <- icc2PL(est.a1, est.b1, ym)
      p2 <- icc2PL(est.a2, est.b2, ym)
      
      
      trA <- sd(est.th1)  / sd(est.th2)
      trB <- mean(est.th1) - trA * mean(est.th2)
      
      #等化の実行
  
      resMS = MS(b1 = est.b1, b2 = est.b2)
      resMSs_list[[s]] <- c(A=resMS$A, B=resMS$B, dA=resMS$dA, dB=resMS$dB, dth = resMS$dth)
      
      
      resMM = MM(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2)
      resMMs_list[[s]] <- c(A=resMM$A, B=resMM$B, dA=resMM$dA, dB=resMM$dB, dth = resMM$dth)
      
      
      resHA = HA(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2, A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resHAs_list[[s]] <- c(fA_init = resHA$fA_init, fB_init = resHA$fB_init ,A=resHA$A, B=resHA$B, dA=resHA$dA, dB=resHA$dB, dth = resHA$dth)
      
      
      resSL = SL(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resSLs_list[[s]] <- c(fA_init = resSL$fA_init, fB_init = resSL$fB_init ,A=resSL$A, B=resSL$B, dA=resSL$dA, dB=resSL$dB, dth = resSL$dth)
      
      
      resRA = RA(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2 ,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5) )
      resRAs_list[[s]] <- c(fA_init = resRA$fA_init, fB_init = resRA$fB_init ,A=resRA$A, B=resRA$B, dA=resRA$dA, dB=resRA$dB, dth = resRA$dth)
      
      
    
      resIC2 = IC2(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2 ,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resIC2s_list[[s]] <- c(fA_init = resIC2$fA_init, fB_init = resIC2$fB_init ,A=resIC2$A, B=resIC2$B, dA=resIC2$dA, dB=resIC2$dB, dth = resIC2$dth)
      
      
      
      resPE = PE(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resPEs_list[[s]] <- c(fA_init = resPE$fA_init, fB_init = resPE$fB_init ,A=resPE$A, B=resPE$B, dA=resPE$dA, dB=resPE$dB, dth = resPE$dth)
      
   
      resPER = PER(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resPERs_list[[s]] <- c(fA_init = resPER$fA_init, fB_init = resPER$fB_init ,A=resPER$A, B=resPER$B, dA=resPER$dA, dB=resPER$dB, dth = resPER$dth)
      
      
      resPES = PES(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resPESs_list[[s]] <- c(fA_init = resPES$fA_init, fB_init = resPES$fB_init ,A=resPES$A, B=resPES$B, dA=resPES$dA, dB=resPES$dB, dth = resPES$dth)
      
      
      resKL = KL(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2 ,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resKLs_list[[s]] <- c(fA_init = resKL$fA_init, fB_init = resKL$fB_init ,A=resKL$A, B=resKL$B, dA=resKL$dA, dB=resKL$dB, dth = resKL$dth)
      
      
      resKLR = KLR(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resKLRs_list[[s]] <- c(fA_init = resKLR$fA_init, fB_init = resKLR$fB_init ,A=resKLR$A, B=resKLR$B, dA=resKLR$dA, dB=resKLR$dB, dth = resKLR$dth)
      
      
      resJS = JS(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5) )
      resJSs_list[[s]] <- c(fA_init = resJS$fA_init, fB_init = resJS$fB_init ,A=resJS$A, B=resJS$B, dA=resJS$dA, dB=resJS$dB, dth = resJS$dth)
      
      
      resHE = HE(a1 = est.a1, b1 = est.b1, a2 = est.a2, b2 = est.b2, ym=ym, H1=H1, H2=H2,A_range = seq(1, 2, by = 0.5),B_range = seq(-3, -2, by = 0.5))
      resHEs_list[[s]] <- c(fA_init = resHE$fA_init, fB_init = resHE$fB_init ,A=resHE$A, B=resHE$B, dA=resHE$dA, dB=resHE$dB, dth = resHE$dth)
      
     
      resCC = CC(r1_cc,r2_cc)
      resCCs_list[[s]] <- c(A=resCC$A, B=resCC$B, dA=resCC$dA, dB=resCC$dB, dth = resCC$dth )
      
      # （岡田追記ここから）――評価に必要な材料を保存
      if (!exists("eval_payloads_list")) eval_payloads_list <- vector("list", nsim)
      eval_payloads_list[[s]] <- list(
        anchor_ids = valid_anchor_IDs,   # 各レプリケートで残ったアンカー集合
        a1 = est.a1, b1 = est.b1,        # Form1 側アンカーの推定 a,b
        a2 = est.a2, b2 = est.b2,        # Form2 側アンカーの推定 a,b
        est_th1 = est.th1,                   # 推定能力（群1）
        est_th2 = est.th2,                   # 推定能力（群2）
        true_th1 = th1,       # ★真値1
        true_th2 = th2 ,       # ★真値2
        true_a1 = true_params1$a1,
        true_b1 = true_params1$b1,
        true_a2 = a2_true,
        true_b2 = b2_true
       
      )
      
      remidx[[paste0("s", s)]] <- c(rem_idx1,rem_idx2)
      
      
    }#繰り返し
    
    remidx_lists[[paste0("m", m, "_i", i)]] <- remidx
    
    #データフレーム化
    resMSs <- do.call(rbind, resMSs_list)
    resMMs <- do.call(rbind, resMMs_list)
    resHAs <- do.call(rbind, resHAs_list)
    resSLs <- do.call(rbind, resSLs_list)
    resRAs <- do.call(rbind, resRAs_list)
    #resICs <- do.call(rbind, resICs_list)
    resIC2s <- do.call(rbind, resIC2s_list)
    resPEs <- do.call(rbind, resPEs_list)
    resPERs <- do.call(rbind, resPERs_list)
    resPESs <- do.call(rbind, resPESs_list)
    resKLs <- do.call(rbind, resKLs_list)
    resKLRs <- do.call(rbind, resKLRs_list)
    resJSs <- do.call(rbind, resJSs_list)
    resHEs <- do.call(rbind, resHEs_list)
    resCCs <- do.call(rbind, resCCs_list)
    
    
    results[[paste0("m", m, "_i", i, "_", format(Sys.time(), "%Y%m%d_%H%M"))]]  <- list(
      resMSs = resMSs,
      resMMs = resMMs,
      resHAs = resHAs,
      resSLs = resSLs,
      resRAs = resRAs,
      #resICs = resICs,
      resIC2s = resIC2s,
      resPEs = resPEs,
      resPERs = resPERs,
      resPESs = resPESs,
      resKLs = resKLs,
      resKLRs = resKLRs,
      resJSs = resJSs,
      resHEs = resHEs,
      resCCs = resCCs,
      eval_payloads = eval_payloads_list   # ←（岡田追記）評価用の材料
    )
    
    
    
    
    dfres <- rbind(
      data.frame(Method="MS", AB=c("A", "B", "dA", "dB","dth"), M=apply(resMSs, 2, mean, na.rm = TRUE), SD=apply(resMSs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="MM", AB=c("A", "B", "dA", "dB","dth"), M=apply(resMMs, 2, mean , na.rm = TRUE), SD=apply(resMMs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="HA", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resHAs, 2, mean, na.rm = TRUE), SD=apply(resHAs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="SL", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resSLs, 2, mean, na.rm = TRUE), SD=apply(resSLs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="RA", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resRAs, 2, mean, na.rm = TRUE), SD=apply(resRAs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      #data.frame(Method="IC", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","ce", "eva_th"), M=apply(resICs, 2, mean, na.rm = TRUE), SD=apply(resICs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="IC2", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resIC2s, 2, mean, na.rm = TRUE), SD=apply(resIC2s, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="PE", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB", "dth"), M=apply(resPEs, 2, mean, na.rm = TRUE), SD=apply(resPEs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="PER", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resPERs, 2, mean, na.rm = TRUE), SD=apply(resPERs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="PES", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resPESs, 2, mean, na.rm = TRUE), SD=apply(resPESs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="KL", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resKLs, 2, mean, na.rm = TRUE), SD=apply(resKLs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="KLR", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resKLRs, 2, mean, na.rm = TRUE), SD=apply(resKLRs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="JS", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB","dth"), M=apply(resJSs, 2, mean, na.rm = TRUE), SD=apply(resJSs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="HE", AB=c("fA_init", "fB_init" , "A", "B", "dA", "dB", "dth"), M=apply(resHEs, 2, mean, na.rm = TRUE), SD=apply(resHEs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim)),
      data.frame(Method="CC", AB=c("A", "B", "dA", "dB","dth"), M=apply(resCCs, 2, mean, na.rm = TRUE), SD=apply(resCCs, 2, sd, na.rm = TRUE)*sqrt((nsim-1)/nsim))
      
    )
    
    true_vals <- c(A=trA, B=trB, dth=0)
    
    init_df <- dfres %>%
      filter(AB %in% c("fA_init", "fB_init")) %>%
      mutate(AB = ifelse(AB=="fA_init","A","B"),
             stat = "Init",
             value = M) %>% 
      dplyr::select(Method, AB, stat, value)
    
    
    
    # BIASの計算
    bias_df <- dfres %>%
      filter(AB %in% c("A", "B")) %>%
      mutate(stat = "BIAS",
             value = M - true_vals[AB]) %>%
      select(Method, AB, stat, value)
    
    # RMSEの計算
    rmse_df <- dfres %>%
      filter(AB %in% c("dA", "dB")) %>%
      mutate(AB = ifelse(AB=="dA", "A", "B"),
             stat = "RMSE",
             value = sqrt(M)) %>%
      select(Method, AB, stat, value)
    
    # SDのデータフレーム作成
    sd_df <- dfres %>%
      filter(AB %in% c("A", "B")) %>%
      mutate(stat = "SD",
             value = SD) %>%
      select(Method, AB, stat, value)
    
    
    
    
    # 全結合
    final_df <- bind_rows(sd_df, rmse_df, bias_df) %>%
      arrange(Method, AB, stat)
    
    result_list[[paste0("i", i, "_m", m)]] <- final_df
    
    init_wide <- init_df %>%
      pivot_wider(
        id_cols = Method,
        names_from = AB,
        values_from = value
      ) %>%
      # stat 列も追加
      mutate(stat = "Init") %>%
      select(Method, stat, A, B) %>%
      arrange(Method)
    
    init_list[[paste0("i", i, "_m", m)]] <- init_wide
    
    
  }#共通項目
  
}#共通受験者

#実行時間
goal<- proc.time()


#結果の出力
# 処理関数を定義(method並び替え，ABごとに列)
# Method の表示順を指定（factor で順序づけ）
method_order <- c("MS", "MM", "HA", "SL","CC", "RA",  "IC2",  "KL", "KLR", "JS", "HE","PE", "PER", "PES")
reshape_df <- function(df) {
  df %>%
    mutate(Method = factor(Method, levels = method_order)) %>%
    filter(stat %in% c( "BIAS", "RMSE", "SD")) %>%
    select(Method, AB, stat, value) %>%
    pivot_wider(names_from = AB, values_from = value) %>%
    arrange(Method, stat)
}

# リスト内の全データフレームに適用
result_list <- lapply(result_list, reshape_df)

# 名前付きリストなら名前も保持
#names(reshaped_list) <- names(result_list)

# タイムスタンプ作成（例: 20250719_1630）
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# 各データフレームをそれぞれ名前付きCSVとして保存
for (name in names(result_list)) {
  
  file_name <- paste0(name, "_", timestamp, ".csv")
  write.csv(result_list[[name]], file = file_name, row.names = FALSE)
}


# 保存ファイル名
timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
filename <- paste0("workspace_", timestamp, ".RData")

# ------- 現在の作業空間を保存-------------
save.image(file = filename)


#---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# データ抽出
df <- result_list[["i3000_m25"]]

# --- 変更点: ここで RA を除外する条件を追加 ---
exclude_list <- c("RA")
# RMSE のみ抽出 (RAを除外)
df_rmse <- df %>% 
  filter(stat == "RMSE", !Method %in% exclude_list)

# BIASのみ抽出 (RAを除外)
df_bias <- df %>% 
  filter(stat == "BIAS", !Method %in% exclude_list)

# ----------------------------------------------

df_bias <- df_bias %>%
  mutate(
    A = abs(A),
    B = abs(B)
  )

df_rmse_long <- df_rmse %>%
  pivot_longer(cols = c(A, B), names_to = "Condition", values_to = "RMSE")

df_bias_long <- df_bias %>%
  pivot_longer(cols = c(A, B), names_to = "Condition", values_to = "BIAS")

# 白にしたいMethod 
white_methods <- c("RA","IC2", "KL", "KLR", "JS", "HE","PE", "PER", "PES")

# fill色を指定するベクトルを作成（白かグレー）
df_bias_long$fill_color <- ifelse(df_bias_long$Method %in% white_methods, "white", "grey70")
df_rmse_long$fill_color <- ifelse(df_rmse_long$Method %in% white_methods, "white", "grey70")

# プロット（RMSE）
ggplot(df_rmse_long, aes(x = Method, y = RMSE, fill = fill_color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_identity() + 
  facet_wrap(~ Condition) +
  labs(title = "RMSE by Method I=3000,anchor=10", x = "Method", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) # X軸ラベルが見やすいように少し回転

# プロット（BIAS）
ggplot(df_bias_long, aes(x = Method, y = BIAS, fill = fill_color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_identity() + 
  facet_wrap(~ Condition) +
  labs(title = "BIAS by Method I=3000,anchor=10", x = "Method", y = "|BIAS|") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) # X軸ラベルが見やすいように少し回転
#tex-------------------------------------



# ---------評価1(dth RMSE)-----------------
# 1. 実データが入っているリストそのものを指定
target_data <- results[[1]] 

# 2. 【重要】target_data（リスト）の要素名から "res" を探す
methods_names <- grep("^res", names(target_data), value = TRUE)

# 確認用：ここで名前が表示されるはずです
print(methods_names) 

# 3. 各手法の平均を計算
dth_means <- sapply(methods_names, function(name) {
  # target_data から各手法の行列/リストを取り出す
  val <- target_data[[name]]
  
  # リスト形式なら行列に変換
  if (is.list(val)) val <- do.call(rbind, val)
  
  # dth列の平均を計算
  if (!is.null(val) && "dth" %in% colnames(val)) {
    mean(val[, "dth"], na.rm = TRUE)
  } else {
    NA
  }
})

# 4. グラフ用データ作成（NAを除去）
plot_data <- data.frame(
  Method = gsub("res", "", names(dth_means)),
  Mean_dth = as.numeric(dth_means)
)
plot_data <- plot_data[!is.na(plot_data$Mean_dth), ]

# 5. 描画
par(mfrow = c(1, 1), mar = c(8, 4, 4, 2))
bp <- barplot(plot_data$Mean_dth, 
              names.arg = plot_data$Method,
              col = "gray",
              las = 2,
              main = "Mean dth Comparison",
              ylab = "Average dth")

text(x = bp, y = plot_data$Mean_dth, label = round(plot_data$Mean_dth, 4), pos = 3, xpd = TRUE)






res3_1 <- genPattern(true_params1$a1, true_params1$b1, th2)
res3_2 <- genPattern(a2_true/trA, trA * b2_true + trB, th2)

est_3_1 <- est2plm(res3_1)
est_3_2 <- est2plm(res3_2)

est.th3_1 <- est_3_1$scores
est.th3_2 <- est_3_2$scores
eq <- trA*est.th2+trB





# --------plink()で等化------------

ip1 <- data.frame(
  a = est.a1,
  b = est.b1
)

ip2 <- data.frame(
  a = est.a2,
  b = est.b2
)

para.on <- list(ip2, ip1)
cat.2 <- rep(x = 2, times = 25)
cat.1 <- rep(x = 2, times = 25)
cat.on <- list(cat.2, cat.1)

pm.old <- as.poly.mod(n = 25)
pm.new <- as.poly.mod(n = 25)

pm.on <- list(pm.old, pm.new)

common.on <- cbind(1:10, 1:10)

pars.on <- as.irt.pars(x = para.on, common = common.on, cat = cat.on, poly.mod = pm.on)

out.on <- plink(x = pars.on, rescale = "MS", symmetric = TRUE)

summary(out.on)


