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
  if(automatic){
    warning('Automatic parameter selection is not avaliable yet')
  }
  gl_model = tryCatch(expr = glarma::glarma(y = y,
                               X = xreg,
                               thetaLags = q,
                               phiLags = p,
                               type = type,
                               method = method,
                               residuals = residuals) ,
                      error = function(err){
                        stop('An error has occurred in the model estimation, try changing the model parameters')
                        return(NA)
                      })
  params = c(p, d, q)

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
  # fit = glarma_model$gl_model |> fitted()
  # e = y - fit

  # Create S3 model object
  # It should be small, but contain everything needed for methods below
  structure(
    list(
      fit = glarma_model$gl_model
    ),
    class = "GLARMA"
  )
}
####################



glarma_call(y = tsibbledata::aus_production$Bricks,
            xreg = NULL, ic = 'AIC', type = 'Poi',
                       p = 2, d = 1, q = 1,
                       method = 'FS', residuals = 'Pearson',
                       automatic = F, trace = F)

tsibbledata::aus_production |>
  model(ing_automatic = GLARMA(Beer,
                                ic = 'AIC',
                                distr = 'Poi',
                                method = 'FS',
                                residuals = 'Pearson',
                                trace = T))
