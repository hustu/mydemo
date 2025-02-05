#' @title High dimensional mediation analysis.
#'
#' @description
#' High dimensional mediation analysis with confounding adjustment.
#' Methods of confounding adjustment include regression adjustment,
#' propensity score adjustment, inverse probability weighting, overlapping weighting.
#'
#' @param X Binary exposure.
#' @param M High dimensional mediators.
#' @param Y Continuous outcome.
#' @param COV Covariates, also called confounders.
#' @param k SIS parameter.
#' @param CONFmethod Method of confounding adjustments, including RA,PS,IPW,OW.
#' @param MTmethod Method of mediation testing, including HIMA, HIMA2.
#' @param tips Whether to output message.
#' @param simplify Only output mediator which p value less than 0.05.
#'
#' @return A dataframe contain mediators' point estimate and p value.
#' @export
#'
#' @examples
#' library(mydemo)
#' dat <- sim_data(n = 300, p = 2000)
#' X <- dat$X
#' Y <- dat$Y
#' M <- dat$M
#' COV <- dat$COV
#' himaconf(X, M, Y, COV, CONFmethod = "OW", MTmethod = "HIMA2")
himaconf <- function(X, M, Y, COV, k=2,
                     CONFmethod = c("RA", "PS", "IPW", "OW"),
                     MTmethod = c("HIMA", "HIMA2"),
                     tips = TRUE, simplify = TRUE) {

  n <- length(X)    # number of samples
  p <- dim(M)[2]    # number of mediators
  q <- dim(COV)[2]  # number of covariates

  # ================================ Step 0: Calculate PS ================================ #
  pscore <- glm(X ~ COV, family = "binomial")[["fitted.values"]]    # propensity scores for exposure X
  wt4x <- ifelse(X==1, 1/pscore, 1/(1-pscore))                      # ipw weights for exposure X
  ow4x <- ifelse(X==1, 1-pscore, pscore)                            # ow weights for exposure X

  # ================================ Step 1: SIS step ================================ #
  if (tips) { message("Step 1: Sure Independent Screening ...", "  (", Sys.time(), ")") }
  d_0 <- k*round(n/log(n))     # SIS筛选出的变量数
  b_SIS <- c()                 # SIS计算出的p值
  for (i in 1:p) {
    if (CONFmethod == "RA") {
      b_SIS[i] <- betaRA(M[,i], X=X, Y=Y, COV=COV, pscore.c = pscore)
    } else if (CONFmethod %in% c("PS", "PSOW1")) {
      b_SIS[i] <- betaPS(M[,i], X=X, Y=Y, COV=COV, pscore.c = pscore)
    } else if (CONFmethod %in% c("IPW", "PSIPW")) {
      b_SIS[i] <- betaIPW(M[,i], X=X, Y=Y, COV=COV, pscore.c = pscore)
    } else if (CONFmethod %in% c("OW", "PSOW")) {
      b_SIS[i] <- betaOW(M[,i], X=X, Y=Y, COV=COV, pscore.c = pscore)
    } else {print('CONFmethod input error!')}
  }

  ID_SIS <- which(b_SIS <= sort(b_SIS)[d_0])   # SIS筛选出的变量序号
  d = length(ID_SIS)                           # SIS筛选出的变量数量

  # --------------- Estimate Total Effect --------------- #
  if (CONFmethod %in% c("IPW")) {
    TE <- summary(glm(Y ~ X, weights = wt4x))$coefficient[2,1]
  } else if (CONFmethod %in% c("OW")) {
    TE <- summary(glm(Y ~ X, weights = ow4x))$coefficient[2,1]
  } else if (CONFmethod %in% c("PS")) {
    TE <- summary(glm(Y ~ X + pscore))$coefficient[2,1]
  } else if (CONFmethod %in% c("RA")) {
    TE <- summary(glm(Y ~ X + COV))$coefficient[2,1]
  } else { print('CONFmethod input error!') }


  if (MTmethod == "HIMA2") {
    # =================== Step 2: DLasso Penalized Estimation =================== #
    if (tips) { message("Step 2: DLasso Penalized Estimation ...", "   (", Sys.time(), ")") }
    if (CONFmethod %in% c("PS")) {
      MXZ <- cbind(M[, ID_SIS], X, pscore)
      DLASSO_fit <- hdi::lasso.proj(x=MXZ, y=Y, family = "gaussian", Z = NULL)
    } else if (CONFmethod %in% c("RA", "IPW", "OW")) {
      MXZ <- cbind(M[, ID_SIS], X, COV)
      DLASSO_fit <- hdi::lasso.proj(x=MXZ, y=Y, family = "gaussian", Z = NULL)
    } else { print('CONFmethod input error!') }
    # The estimator for beta
    est_beta <- DLASSO_fit$bhat[1:d]
    se_beta <- DLASSO_fit$se[1:d]
    p_beta <- DLASSO_fit$pval[1:d]


    # ========================= Step 3: Mediation test ========================= #
    if (tips) { message("Step 3: Mediation test ...", "     (", Sys.time(), ")") }

    # --------------- Estimate alpha --------------- #
    est_alpha <- se_alpha <- p_alpha <- c()
    # XZ <- data.frame(X=X, COV=COV)
    for (i in 1:d){
      if (CONFmethod %in% c("RA")) {
        fit_a <- summary(glm(M[,ID_SIS[i]] ~ X + COV, family = "gaussian"))$coefficients
      } else if (CONFmethod %in% c("PS"))  {
        fit_a <- summary(glm(M[,ID_SIS[i]] ~ X + pscore, family = "gaussian"))$coefficients
      } else if (CONFmethod %in% c("IPW")) {
        fit_a <- summary(glm(M[,ID_SIS[i]] ~ X, family = "gaussian", weights = wt4x))$coefficients
      } else if (CONFmethod %in% c("OW"))  {
        fit_a <- summary(glm(M[,ID_SIS[i]] ~ X, family = "gaussian", weights = ow4x))$coefficients
      } else { print('CONFmethod input error!') }
      # The estimator for alpha
      est_alpha[i] <- fit_a[2, 1]
      se_alpha[i] <- fit_a[2, 2]
      p_alpha[i] <- fit_a[2, 4]
    }

    # --------------- Estimate beta --------------- #
    # De-biased Lasso Estimates

    # --------------- Mediation test --------------- #
    # Joint Significance Test whit mixture distribution
    PA <- cbind(p_alpha, p_beta)
    P_raw <- apply(PA, 1, max)
    N0 <- dim(PA)[1]*dim(PA)[2]
    input_pvalues <- PA + matrix(runif(N0, 0, 10^(-10)), dim(PA)[1], 2)
    nullprop <- HDMT::null_estimation(input_pvalues)
    fdrcut <- HDMT::fdr_est(nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
                      nullprop$alpha1, nullprop$alpha2,
                      input_pvalues,
                      exact=0)

    IDE <- est_beta*est_alpha             # mediation/indirect effect
    # direct effect
    DE <- DLASSO_fit$bhat[d+1]
    DE_se <- DLASSO_fit$se[d+1]
    DE_p <- DLASSO_fit$pval[d+1]

    output <- data.frame(alpha = est_alpha, alpha_se = se_alpha, alpha_p = p_alpha,
                         beta = est_beta, beta_se = se_beta, beta_p = p_beta,
                         IDE = IDE, P_fdr = fdrcut, P_test = P_raw, Med_prop = IDE/TE*100,
                         DE = DE, DE_se = DE_se, DE_p = DE_p, TE = TE)
    if (simplify) {
      ID_fdr <- which(fdrcut <= 0.05)
      output <- output[ID_fdr, ]
    }     # 简化输出,只输出p<0.05的结果

  } else if (MTmethod == "HIMA") {
    # ====================== Step 2: MCP Penalized Estimation ====================== #
    if (tips) { message("Step 2: MCP Penalized Estimation ...", "   (", Sys.time(), ")") }
    if (CONFmethod %in% c("RA", "OW", "IPW")) {
      MXZ <- cbind(M[, ID_SIS], X, COV)
      MCPfit <- ncvreg::ncvreg(MXZ, Y, family = "gaussian", penalty = "MCP",
                       penalty.factor = c(rep(1, d), rep(0, 1+q)))
    } else if (CONFmethod %in% c("PS")) {
      MXZ <- cbind(M[, ID_SIS], X, pscore)
      MCPfit <- ncvreg::ncvreg(MXZ, Y, family = "gaussian", penalty = "MCP",
                       penalty.factor = c(rep(1, d), rep(0, 1+1)))
    } else {
      print('Method input error')
    }
    lam <- MCPfit$lambda[which.min(BIC(MCPfit))]   # lambda for MCP penalized regression
    est <- coef(MCPfit, lambda = lam)[2:(d + 1)]
    ID_non <- which(est != 0)
    ID_test <- ID_SIS[ID_non]     # index of non-zero variables in M
    d1 <- length(ID_test)         # SIS筛选出的变量数量

    # --------------- Estimate alpha --------------- #
    est_alpha <- se_alpha <- p_alpha <- c()
    XZ <- data.frame(X=X, COV=COV)
    for (i in 1:d1){
      if (CONFmethod == "RA") {
        fit_a <- summary(glm(M[,ID_test[i]] ~ X + COV, family = "gaussian"))$coefficients
      } else if (CONFmethod %in% c("PS"))  {
        fit_a <- summary(glm(M[,ID_test[i]] ~ X + pscore, family = "gaussian"))$coefficients
      } else if (CONFmethod %in% c("IPW")) {
        fit_a <- summary(glm(M[,ID_test[i]] ~ X, family = "gaussian", weights = wt4x))$coefficients
      } else if (CONFmethod %in% c("OW")){
        fit_a <- summary(glm(M[,ID_test[i]] ~ X, family = "gaussian", weights = ow4x))$coefficients
      } else {
        print("CONFmethod input error!")
      }
      est_alpha[i] <- fit_a[2, 1]
      se_alpha[i] <- fit_a[2, 2]
      p_alpha[i] <- fit_a[2, 4]
    }

    # --------------- Estimate beta --------------- #
    if (CONFmethod %in% c("PS")) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV = pscore)
      fit_b <- summary(glm(Y ~ ., family = "gaussian", data = YMX))$coefficients
    } else if (CONFmethod %in% c("RA", "OW", "IPW")) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV = COV)
      fit_b <- summary(glm(Y ~ ., family = "gaussian", data = YMX))$coefficients
    } else {
      print("CONFmethod input error!")
    }
    # the estimator for beta
    est_beta <- fit_b[2:(d1+1), 1]
    se_beta <- fit_b[2:(d1+1), 2]
    p_beta <- fit_b[2:(d1+1), 4]

    # --------------- Mediation test --------------- #
    # Joint Significance Test BH/bon adjustment
    # Use the maximum value as p value
    P_BH_beta <- p.adjust(p_beta, "BH", d)   # BH/fdr adjustd p-value
    P_BH_alpha <- p.adjust(p_alpha, 'BH', d)
    P_joint_BH <- apply(rbind(P_BH_beta, P_BH_alpha), 2, max)

    P_bon_beta <- p.adjust(p_beta, "bonferroni", d)   # bonferroni adjustd p-value
    P_bon_alpha <- p.adjust(p_alpha, "bonferroni", d)
    P_joint_bon <- apply(rbind(P_bon_beta, P_bon_alpha), 2, max)

    IDE <- est_beta*est_alpha             # mediation/indirect effect
    DE <- fit_b[d1+2, 1]                  # direct effect
    DE_se <- fit_b[d1+2, 2]
    DE_p <- fit_b[d1+2, 4]

    output <- data.frame(alpha = est_alpha, alpha_SE = se_alpha, p_alpha = p_alpha,
                         beta = est_beta, beta_SE = se_beta, p_beta = p_beta,
                         IDE = IDE, P_joint_BH = P_joint_BH, P_joint_bon = P_joint_bon,
                         Med_prop = IDE/TE*100,
                         DE = DE, DE_SE = DE_se, DE_p=DE_p, TE = TE)
    if (simplify) {
      ID_sig <- which(P_joint_bon<0.05 | P_joint_BH<0.05)
      output <- output[ID_sig,]
    }      # 简化输出,只输出p<0.05的结果

  } else {
    print('MTmethod input error!')
  }

  if (tips) { message("Done!", "     (", Sys.time(), ")") }
  return(output)
}
