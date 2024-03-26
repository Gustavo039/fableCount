globalVariables(c("p", "q", "constant"))
library(fpp3)
library(tscount)
library(tidyverse)
library(rlang)


glarma_call = function(y = y, ic = ic, type = type,
                       p = 0, d = 0, q = 0,
                       method = method, residuals = residuals,
                       automatic = automatic, trace = trace, xreg = NULL){
  if(is.null(xreg) == F)
      xreg = cbind(rep(1, length(y)), xreg)
     else
       xreg = matrix(rep(1, length(y)), ncol = 1)
  colnames(xreg)[1] = 'Intercept'
  params = c(p, d, q)
  if(automatic){
    warning('Automatic parameter selection is not avaliable yet')
  }
  if(params[1] == 0) past_obs = NULL
  else past_obs = 1:params[1]
  if(params[3] == 0) past_mean = NULL
  else past_mean = 1:params[3]
  gl_model = tryCatch(expr = glarma::glarma(y = y,
                               X = xreg,
                               thetaLags = past_mean,
                               phiLags = past_obs,
                               type = type,
                               method = method,
                               residuals = residuals) ,
                      error = function(err){
                        stop('An error has occurred in the model estimation, try changing the model parameters')
                        return(NA)
                      })

print(gl_model)
return(list(params = params, gl_model = gl_model))
}


#' Estimate a GLARMA model
#'
#' Estimate Generalized Linear Autoregressive Moving Average  model
#' with Poisson or Negative Binomial distribution.
#' Also is provide a automatic parameter algorithm selection for the Autorregressive and Moving Avarege params
#'
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic Character, can be 'AIC','BIC'. The information criterion used in selecting the model.
#' @param distr Character, can be 'poisson' or 'nbinom'. The probabilty distribution used for the generalized model
#' @param trace Logical. If the automatic parameter algorithm is runnig, print the path to the best model estimation
#' @param check_residuals Logical. If the automatic parameter algorithm is runnig, it checks if the residuals respect the model assumptions
#'
#' @section Specials:
#'
#' \subsection{pq}{
#' Also called p:Autoregressive and q:Moving Avarages,
#' the pq can be define by the user,
#' or if it's omited  the automatic parameter selection algorithm is trigered.
#'
#' The automatic parameter selection algorithm gonna fit the best model based on the information criterion
#'
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an INGARCH model without explicitly using the `xreg()` special.
#' Common exogenous regressor specials as specified in [`common_xregs`] can also be used.
#' These regressors are handled using [stats::model.frame()],
#' and so interactions and other functionality behaves similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`],
#' where including `1` will add a constant and `0` or `-1` will remove the constant.
#' If left out, the inclusion of a constant will be determined by minimising `ic`.
#'
#' If a xreg is provided, the model forecast is not avaliable
#'
#' \preformatted{
#' xreg(..., fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)\cr
#'   `fixed`    \tab A named list of fixed parameters for coefficients. The names identify the coefficient, and should match the name of the regressor. For example, `fixed = list(constant = 20)`.
#' }
#' }
#'
#' @return A model specification.
#' @export
GLARMA = function(formula,
                   ic = c('AIC', 'BIC'),
                   distr = c('Poi', 'NegBin'),
                   method = c('FS', 'NR'),
                   residuals = c('Pearson', 'Score'),
                   trace = F) {

  ic = match.arg(ic)
  distr = match.arg(distr)
  method = match.arg(method)
  residuals = match.arg(residuals)


  model_GLARMA = new_model_class("GLARMA",
                                  train = train_GLARMA,
                                  specials = specials_GLARMA,
                                  check = function(.data) {
                                    if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                  }
  )

  # Return a model definition which stores the user's model specification
  new_model_definition(model_GLARMA, {{formula}},
                       ic = ic,
                       distr = distr,
                       method = method,
                       residuals = residuals,
                       trace = trace)
}

specials_GLARMA = new_specials(
  pdq = function(p = 'not choosen', d = 'not choosen', q = 'not choosen',
                p_init = 2, q_init = 2,
                fixed = list()) {

    if(!all(grepl("^(ma|ar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for pq() must begin with ar or ma, followed by a lag number.")
    }
    as.list(environment())
  },
  common_xregs,
  xreg = function(..., fixed = list()) {
    dots = enexprs(...)
    env = map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] = NULL
    env = if (!is_empty(env)) get_env(env[[1]]) else base_env()

    constants = map_lgl(dots, inherits, "numeric")
    constant_forced = any(map_lgl(dots[constants], `%in%`, 1))

    model_formula = new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )

    # Mask user defined lag to retain history when forecasting
    env = env_bury(env, lag = lag)

    xreg = model.frame(model_formula, data = env, na.action = stats::na.pass)
    tm = terms(xreg)
    constant = as.logical(tm %@% "intercept")
    xreg = model.matrix(tm, xreg)

    if (constant) {
      xreg = xreg[, -1, drop = FALSE]
    }

    list(
      constant = if (!constant || constant_forced) constant else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg,
      fixed = fixed
    )
  },
  .required_specials = "pdq",
  .xreg_specials = names(common_xregs)
)


train_GLARMA = function(.data, specials, ic,
                         distr, residuals,
                         method, trace, ...){
  mv = tsibble::measured_vars(.data)
  if(length(mv) > 1) stop("GLARMA is a univariate model.")
  y = .data[[mv]]


  automatic = F
  if(specials$pdq[[1]]$p |>
     is.character()) automatic = T
  xreg = specials$xreg
  glarma_model = glarma_call(y = y,
                               p = specials$pdq[[1]]$p,
                               d = specials$pdq[[1]]$d,
                               q = specials$pdq[[1]]$q,
                               ic = ic,
                               type = distr,
                               method = method,
                               residuals = residuals,
                               automatic = automatic,
                               trace = trace,
                               xreg = xreg[[1]]$xreg)

  # Compute fitted values and residuals
   fit = glarma_model$gl_model |> fitted()
   e = y - fit

  # Create S3 model object
  # It should be small, but contain everything needed for methods below
   structure(
     list(
       coef = glarma_model$params,
       gl_model = glarma_model$gl_model,
       distr = distr,
       residuals_type = residuals,
       method = method,
       n = length(y),
       y_name = mv,
       fitted = fit,
       residuals = e,
       sigma2 = var(e, na.rm = TRUE)
     ),
    class = "GLARMA"
  )
}

model_sum.GLARMA = function(x){
  if(x$gl_model$r == 1) out = sprintf("GLARMA(%i, %i)", x$coef[1],x$coef[3])
  else out = sprintf("GLARMA(%i, %i) w/ covariates", x$coef[1], x$coef[3])
  out
}



split_tibble = function(input_tibble) {
  # Calculate the number of rows needed
  num_rows = nrow(input_tibble) %/% 4

  # Initialize an empty list to store rows
  rows = list()

  # Loop through each row and extract the data for 4 columns
  for (i in 1:num_rows) {
    start_row = (i - 1) * 4 + 1
    end_row = min(i * 4, nrow(input_tibble))
    rows[[i]] = input_tibble[start_row:end_row, 1]
  }

  # Pad rows with NA values if necessary to ensure consistent column lengths
  max_length = max(sapply(rows, length))
  rows = sapply(rows, function(row) c(row, rep(NA, max_length - length(row))))

  # Create a tibble from the list of rows
  result_tibble = as_tibble(do.call(cbind, rows))

  return(result_tibble)
}

#' @export
report.GLARMA = function(x){
  if(x$distr == 'Poi')
    distr_x = 'Poisson'
  else distr_x = 'Negative Binomial'
  x_tidy = x |> tidy()
  cat('\n')
  cat(sprintf("%s GLARMA(%i, %i)", distr_x, x$coef[1], x$coef[3]))
  cat('\n')
  x_tidy |> dplyr::filter(statistic == 'estimate' | statistic == 'std_error') |> print()
  cat('\n')
  cat(paste('log likelihood='), x$gl_model$logLik, sep = '')
  cat('\n')
  cat(paste('AIC='), x$gl_model$aic, sep = '')


}


#' Tidy a fable model
#'
#' Returns the coefficients from the model in a `tibble` format.
#'
#' @inheritParams generics::tidy
#'
#' @return The model's coefficients in a `tibble`.
#'
#' @examples
#' @export
tidy.GLARMA = function(x){
  sum_model = x$gl_model |> summary()
  out = sum_model[14 : (14 + x$gl_model$pq)] |>
    unlist() |>
    tibble::as_tibble()
  if(x$coef[1] == 0) ar_param = NULL
   else ar_param = paste('ar_', rep(1:x$coef[1], 4) |> sort(), sep = '')
  if(x$coef[3] == 0) ma_param = NULL
    else ma_param = paste('ma_', rep(1:x$coef[3], 4) |> sort(), sep = '')
  if(nrow(sum_model$coefficients1) > 1){
  coef1 =
    c(
      rep('intercept', 4),
      paste('xreg_',
            rep(1:(nrow(sum_model$coefficients1) - 1), 4) |>
              sort(), sep = '')
      )
  }
  else coef1 = rep('intercept', 4)
  out = out |>
    dplyr::mutate(term =
                    c(coef1, ar_param, ma_param),
                  statistic =
                    rep(c('estimate', 'std_error', 'z_ratio', 'p_value'),
                        nrow(sum_model$coefficients1)+x$gl_model$pq)
    ) |>
    tidyr::pivot_wider(values_from = value, names_from = term)


}


#' Extract fitted values from a fable model
#'
#' Extracts the fitted values.
#'
#' @inheritParams forecast.INGARCH
#'
#' @return A vector of fitted values.
#'
#' @examples
#' @export
fitted.GLARMA = function(object, ...){
  object$fitted
}


#' Extract residuals from a fable model
#'
#' Extracts the residuals.
#'
#' @inheritParams forecast.INGARCH
#' @param type The type of residuals to extract.
#'
#' @return A vector of fitted residuals.
#'
#' @examples
#' @export
residuals.GLARMA = function(object, ...){
  object$residuals
}

glance.GLARMA = function(x, ...){
  tibble::tibble(sigma2 = x$sigma2,
                 log_lik = x$tsmodel$logLik,
                 AIC = x$tsmodel |> AIC(),
                 BIC = x$tsmodel |> BIC(),
  )
}




#' Forecast a model from the fable package
#'
#' Produces forecasts from a trained model.
#'
#' Predict future observations based on a fitted GLM-type model for time series of counts.
#' For 1 step ahead, it returns parametric forecast, based on the Poisson distribution,
#' for multiples steps forecast, the distribution is not know analytically, so it uses a parametric bootstrap
#'
#' @inheritParams generics::forecast
#' @param new_data Tsibble, it has to contains the time points and exogenous regressors to produce forecasts for.
#' @param bootstrap Logical, if `TRUE`, then forecast distributions are computed using simulation with resampled errors.
#' @param times Numeric, the number of sample paths to use in estimating the forecast distribution when `bootstrap = TRUE`.
#'
#' @importFrom stats formula residuals
#'
#' @return A list of forecasts.
#'
#' @examples
#' @export
forecast.GLARMA = function(object, new_data,...){
  h = NROW(new_data)
  newdata = matrix(rep(1, h), ncol = 1)
  newoffset  = matrix(rep(0, h), ncol = 1)
  newm = matrix(rep(1, h), ncol = 1)
  values = glarma::forecast(object = object$gl_model,
                            n.ahead = h,
                            newdata = newdata,
                            newoffset = newoffset)
  distributional::dist_poisson(values$Y)
}

####################

xreg_teste = cbind(
  rep(1, length(tsibbledata::aus_production$Beer)),
  tsibbledata::aus_production$Gas) |> as.model.matrix()

xreg_teste = matrix(c(rep(1, 218),
                      tsibbledata::aus_production$Gas), ncol = 2
)


teste_glarma = glarma_call(y = tsibbledata::aus_production$Beer,
            xreg = tsibbledata::aus_production$Gas, ic = 'AIC', type = 'Poi',
                       p = 0, d = 1, q = 1,
                       method = 'FS', residuals = 'Pearson',
                       automatic = F, trace = F)

teste_sum = teste_glarma$gl_model |> summary()

teste_fab = tsibbledata::aus_production |>
  model(ing_automatic = GLARMA(Beer ~ pdq(1,1,0) + Gas,
                                ic = 'AIC',
                                distr = 'Poi',
                                method = 'FS',
                                residuals = 'Pearson',
                                trace = T),
        ar = ARIMA(Beer ~ pdq(1,0,1)))

teste_fab |> select(ing_automatic) |> tidy()
teste_fab |> select(ing_automatic) |> report()
teste_fab |> select(ar) |> report()

