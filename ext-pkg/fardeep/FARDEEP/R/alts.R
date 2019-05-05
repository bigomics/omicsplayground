#' Using the basic idea of least trimmed square to detect and remove outliers before estimating the coefficients.
#' Adaptive least trimmed square.
#'
#' @param x input matrix of predictors with n rows and p columns.
#' @param y input vector of dependent variable with length n.
#' @param alpha1  parameter used to adjust the upper bound of outliers. Take value from 0 to 1, default 0.1.
#' @param alpha2  parameter used to adjust the lower bound of outliers. Take value larger than 1, default 1.5.
#' @param k  parameter used to determine the boundary of outliers in the following step of algorithm. Take value from 1 to 10, default 6.
#' @param nn  whether coefficients are non-negative,default TRUE.
#' @param intercept whether intercept is included in model, default TRUE.
#' @return beta:  estimation of coefficients.
#' @return number_outlier: number of outliers.
#' @return outlier_detect:  index of detected outliers.
#' @return X.new:  good observed points for independent variables.
#' @return Y.new:  good observed points for dependent variables.
#' @return k:  modified k (if the input value is not appropriate).
#' @author YuningHao, Ming Yan, Yu Lei and YuyingXie
#' @references YuningHao, Ming Yan, Yu Lei and YuyingXie. Fast and Robust Deconvolution of Tumor Infiltrating Lymphocyte from Expression Profiles using Least Trimmed Squares.
#' @examples
#' library(FARDEEP)
#' samp = sample.sim(n = 500, p = 20, sig = 1, a1 = 0.1, a2 = 0.2, nn = TRUE, intercept = TRUE)
#' result = alts(samp$x, samp$y, alpha1 = 0.1, alpha2 = 1.5, k = 6, nn = TRUE, intercept = TRUE)
#' coef = result$beta
#' @export
alts = function (x,  y, alpha1 = 0.1, alpha2 = 1.5, k = 6, nn = TRUE, intercept = TRUE){

  if(nn){
    if(intercept){
      m1 = nnls::nnls (cbind(1, x), y)
    }else{
      m1 = nnls::nnls (x, y)
    }
  }else{
    if(intercept){
      m1 = stats::lm (y ~ x)
    }else{
      m1 = stats::lm (y ~ x - 1)
    }
  }
  res = abs (stats::resid(m1))
  n   = length (y)
  order_id = order (res, decreasing = F)
  res_int  = res[order_id]
  Y_int = y[order_id]
  X_int = x[order_id, ]
  index1 = min (which (res_int > 1 * stats::median (res_int)))
  k_up   = n - index1
  y_out_up = Y_int[c ((n - k_up + 1) : n)]
  k.low.ex = alpha1 * k_up
  k_low    = ceiling (k.low.ex)
  out_id   = (1:n > (n - k_low))
  kep_id   = (1:n < (n - k_low + 1))
  Y_new    = Y_int[kep_id]
  X_new    = X_int[kep_id, ]
  repeat{
    if(nn){
      if(intercept){
        m_new = nnls::nnls (cbind(1, X_new), Y_new)
        beta_new = m_new$x
        Y_hat    = cbind(1, x) %*% beta_new
      }else{
        m_new = nnls::nnls (X_new, Y_new)
        beta_new = m_new$x
        Y_hat    = x %*% beta_new
      }
    }else{
      if(intercept){
        m_new = stats::lm (Y_new ~ X_new)
        beta_new = m_new$coefficients
        Y_hat    =  cbind (1, x) %*% beta_new
      }else{
        m_new = stats::lm (Y_new ~ X_new - 1)
        beta_new = m_new$coefficients
        Y_hat    = x %*% beta_new
      }
    }
    res_new  = abs (y - Y_hat)
    order_id = order (res_new, decreasing = F)
    res_new_ord = res_new[order_id]
    Y_ord    = y[order_id]
    X_ord    = x[order_id, ]
    index1   = min (which (res_new_ord > k * stats::median (res_new_ord)))
    temp     = n - index1
    if(temp <= 0) {
      while (temp <= 0) {
        k = k - 1
        index1 = min (which (res_int > k * median (res_int)))
        temp   = n - index1
      }
    }
    k_up     = min (temp, k_up)
    k.low.ex = alpha2 * k.low.ex
    k_low    = min (ceiling(k.low.ex), k_up)
    out_id   = (1:n > (n - k_low))
    kep_id   = (1:n < (n - k_low + 1))
    Y_new    = Y_ord[kep_id]
    X_new    = X_ord[kep_id, ]
    if (prod(out_id == ((1:n) > n-k_up)) == 1){
      break
    }
  }
  if(nn){
    if(intercept){
      model = nnls::nnls (cbind(1, X_new), Y_new)
    }else{
      model = nnls::nnls (X_new, Y_new)
    }
    coefficients = model$x
  }else{
    if(intercept){
      model = stats::lm (Y_new ~ X_new)
    }else{
      model = stats::lm (Y_new ~ X_new - 1)
    }
    coefficients = model$coefficients
  }
  number_outlier = sum(out_id)
  outlier_id     = order_id[out_id]
  result         = list(beta = coefficients, number_outlier = number_outlier, outlier_detect = outlier_id,
                        X.new = X_new, Y.new = Y_new, k = k)
  return(result)
}
