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

report.GLARMA = function(x){
  cat(print(x$gl_model$delta[1,1]))
  cat('\n')
  cat(sprintf("%s GLARMA(%i, %i)", x$distr, x$coef[1], x$coef[2]))
  cat(print('1'))
  cat('\n')
  cat(sprint("Intercept: %i", x$gl_model$delta[1,1] |> round(0)))
  # cat('\n')
  # if(x$coef[1] != 0){
  #   cat(sprint('AR Params'))
  #   print.default(
  #     setNames(x$gl_model$delta[1:x$coef[1],1], paste0("AR", seq_len(x$coef[1]))),
  #     print.gap = 2
  #     )
  # }
  # if(x$coef[3] != 0){
  #   cat(sprint('MA Params'))
  #   print.default(
  #     setNames(x$gl_model$delta[x$coef[1]+1:x$coef[3]+1,1], paste0("AR", seq_len(x$coef[3]))),
  #     print.gap = 2
  #   )
  # }
}


report.GLARMA = function(x){
  cat('\n')
  cat(sprintf("%s GLARMA(%i, %i)", x$distr, x$coef[1], x$coef[3]))
  cat('\n')
  cat(print(x$gl_model$delta))
}

tidy.GLARMA = function(x){
  model_summary = x$gl_model |>
    summary()
  model_summary$coefficients1 |> print()

  for (i in 1:(model_summary$pq+1)) {
    coef_name = paste("model_summary$coefficients", i, sep = "")
    coefficients = rbind(coefficients, get(coef_name))
  }

}

fitted.GLARMA = function(object, ...){
  object$fitted
}

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



teste_glarma = glarma_call(y = tsibbledata::aus_production$Beer,
            xreg = NULL, ic = 'AIC', type = 'Poi',
                       p = 0, d = 1, q = 1,
                       method = 'FS', residuals = 'Pearson',
                       automatic = F, trace = F)

teste_fab = tsibbledata::aus_production |>
  mutate(Beer = (Beer/20) |> ceiling()) |>
  model(ing_automatic = GLARMA(Beer ~ pdq(0,1,1),
                                ic = 'AIC',
                                distr = 'Poi',
                                method = 'FS',
                                residuals = 'Pearson',
                                trace = T))

teste_fab |> select(ing_automatic) |> forecast(h = 11) |> autoplot(tsibbledata::aus_production |>
                                                                     mutate(Beer = (Beer/20) |> ceiling()))
