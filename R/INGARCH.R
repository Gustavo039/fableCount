globalVariables(c("p", "q", "constant"))
library(fpp3)
library(tscount)
library(tidyverse)
library(rlang)

ingarch_tscall = function(x, p = 0, q = 0,
                          automatic = T, trace = T,
                          ic, link, distr){
  if(automatic){
    # search_matrix storage a INGARCH(p,q) model in a matrix[p,q]
    search_matrix = matrix(ncol = 3, nrow = 3)
    for(i in 0:2){
      for(j in 0:2){
        step_model =
          tryCatch(expr =
                     {tscount::tsglm(x,
                                    model = list(past_obs = 0:i,
                                                 past_mean = 0:j),
                                    link = link, distr = distr)},
                  error = function(err){
                    return(NA)
                  })
        search_matrix[i + 1, j + 1] = ifelse(is.na(get(ic)(step_model)) == F, get(ic)(step_model), NA)
        if(trace){
          cat(sprintf('INGARCH[%d, %d], %s:%.2f \n', i, j, ic, search_matrix[i+1,j+1]))
        }
      }
    }
    params = which(search_matrix == min(search_matrix), arr.ind = TRUE)
    params = params - 1
    tscount_model = tscount::tsglm(x,
                                   model =
                                     list(past_obs = params[1],
                                          past_mean = params[2]),
                                   link = link,distr = distr
                                   )
  }
  else{
    tscount_model = tscount::tsglm(x,
                                   model =
                                     list(past_obs = p$p,
                                          past_mean = p$q),
                                   link = link, distr = distr
    )
    params = c(p$p,p$q)
  }
  return(list(params = params, tscount_model = tscount_model))
}

INGARCH = function(formula,
                    ic = c('AIC', 'BIC', 'QIC'),
                    link = c('identity', 'log'),
                    distr = c('poisson', 'nbinom'),
                    trace = F) {

  ic = match.arg(ic)
  link = match.arg(link)
  distr = match.arg(distr)

  model_INGARCH = new_model_class("INGARCH",
                                 train = train_INGARCH,
                                 specials = specials_INGARCH,
                                 check = function(.data) {
                                   if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                 }
  )

  # Return a model definition which stores the user's model specification
  new_model_definition(model_INGARCH, {{formula}},
                       ic = ic,
                       link = link,
                       distr = distr,
                       trace = trace)
}

specials_INGARCH = new_specials(
  pq = function(p = 'not choosen', q = 'not choosen',
                 p_init = 2, q_init = 2,
                 fixed = list()) {
    # if (self$stage %in% c("estimate", "refit")) {
    #   p <- p[p <= floor(NROW(self$data) / 3)]
    #   q <- q[q <= floor(NROW(self$data) / 3)]
    # }
    # p_init <- p[which.min(abs(p - p_init))]
    # q_init <- q[which.min(abs(q - q_init))]
    if(!all(grepl("^(ma|ar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for pq() must begin with ar or ma, followed by a lag number.")
    }
    as.list(environment())
  },
  common_xregs,
  xreg = function(..., fixed = list()) {
    dots <- enexprs(...)
    env <- map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
    env <- if (!is_empty(env)) get_env(env[[1]]) else base_env()

    constants <- map_lgl(dots, inherits, "numeric")
    constant_forced <- any(map_lgl(dots[constants], `%in%`, 1))

    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )

    # Mask user defined lag to retain history when forecasting
    env <- env_bury(env, lag = lag)

    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    tm <- terms(xreg)
    constant <- as.logical(tm %@% "intercept")
    xreg <- model.matrix(tm, xreg)

    if (constant) {
      xreg <- xreg[, -1, drop = FALSE]
    }

    list(
      constant = if (!constant || constant_forced) constant else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg,
      fixed = fixed
    )
  },
  .required_specials = "pq",
  .xreg_specials = names(common_xregs)
)

train_INGARCH = function(.data, specials, ic,
                         link, distr,
                         trace, ...){
  # Extract a vector of response data
  mv = tsibble::measured_vars(.data)
  if(length(mv) > 1) stop("INGARCH is a univariate model.")
  y = .data[[mv]]

  # Pull out inputs from the specials
  if(length(specials$season) > 1) stop("The `season()` special of `SMEAN()` should only be used once.")
  m = specials$season[[1]]

  automatic = F
  if(specials$pq[[1]]$p |>
     is.character()) automatic = T
  tsglm_model = ingarch_tscall(x = y,
                               ic = ic,
                               p = specials$pq[[1]],
                               q = specials$pq[[1]],
                               link = link,
                               distr = distr,
                               automatic = automatic,
                               trace = trace)

  # Compute fitted values and residuals
  fit = tsglm_model$tscount_model |> fitted()
  e = y - fit

  # Create S3 model object
  # It should be small, but contain everything needed for methods below
  structure(
    list(
      coef = tsglm_model$params,
      tsmodel = tsglm_model$tscount_model,
      distr = distr,
      link = link,
      n = length(y),
      y_name = mv,
      fitted = fit,
      residuals = e,
      sigma2 = var(e, na.rm = TRUE)
    ),
    class = "INGARCH"
  )
}

model_sum.INGARCH = function(x){
  sprintf("INGARCH(%i, %i)", x$coef[1],x$coef[2])
}

report.INGARCH = function(x){
  cat('\n')
  cat(sprintf("%s INGARCH[%i, %i] w/ %s link", x$distr, x$coef[1], x$coef[2], x$link))
  cat('\n')
  cat(print(x$tsmodel$coefficients))
}

tidy.INGARCH = function(x){
  model_summary = x$tsmodel |>
    summary() |>
    {\(x)x$coefficients}() |>
    tibble::rownames_to_column() |>
    tibble::as_tibble() |>
    dplyr::rename(term = 1,
                  estimate = 2,
                  std.error = 3) |>
    dplyr::mutate(term = stringr::str_replace(term, '(Intercept)', 'constant'),
                  term = stringr::str_replace(term, 'beta', 'ar'),
                  term = stringr::str_replace(term, 'alpha', 'ma'))


  return(model_summary)
}

fitted.INGARCH = function(object, ...){
  object$fitted
}

residuals.INGARCH = function(object, ...){
  object$residuals
}

#################

teste_tscount = tscount::tsglm(tscount::campy, model = list(past_obs = 5, past_mean = 1:10),
                               link = 'identity',
                               distr = 'poisson')

teste_me = ingarch_tscall(tsibbledata::aus_production$Beer, ic = 'QIC', distr = 'poisson', link = 'identity')


teste_ingarch = tsibbledata::aus_production |>
  model(ing = INGARCH(Beer, ic = 'BIC', link = 'identity', distr = 'poisson', trace = T, automatic = T),
        ar = ARIMA(Beer ~ pdq(1,0,0) + PDQ(0,0,0)))

tsibbledata::aus_production |>
  model(ing = INGARCH(Beer ~ pq(0,2), link = 'identity', distr = 'poisson'),
        ing_automatic = INGARCH(Beer, ic = 'BIC', link = 'identity', distr = 'poisson', trace = T))

tsibbledata::aus_production |>
  model(ing_automatic = INGARCH(Beer, ic = 'BIC', link = 'identity', distr = 'poisson', trace = T))
