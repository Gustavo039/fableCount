glarma_call = function(){
  if(automatic){}
  else{
    expla_var =
      data.frame('(Intercept)' =
                   rep(1, y |>
                         length()
                       )
                 ) |>
      as.matrix()
    gl_model = glarma::glarma(y = y,
                                 X = expla_var,
                                 thetaLags = q,
                                 phiLags = p,
                                 type = distr,
                                 method = method,
                                 residuals = residuals)
    params = c(p,q)
  }
  return(list(params = params, gl_model = gl_model))
}

GLARMA = function(formula,
                   ic = c('AIC', 'BIC'),
                   distr = c('poisson', 'nbinom'),
                   automatic = T,
                   trace = F) {

  ic = match.arg(ic)
  link = match.arg(link)
  distr = match.arg(distr)

  model_INGARCH = new_model_class("GLARMA",
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
                       automatic = automatic,
                       trace = trace)
}


specials_GLARMA = new_specials(
  season = function(period = NULL) {
    # Your input handling code here.
    get_frequencies(period, self$data, .auto = "smallest")
  },
  xreg = function(...) {
    # This model doesn't support exogenous regressors, time to error.
    stop("Exogenous regressors aren't supported by `SMEAN()`")
  },
  # This model requires `season()`
  # Adding this allows `SMEAN(y)` to automatically include the `season()` special
  .required_specials = "season"
)

train_GLARMA = function(.data, specials, ic,
                         link, distr, automatic,
                         trace, ...){
  # Extract a vector of response data
  mv = tsibble::measured_vars(.data)
  if(length(mv) > 1) stop("INGARCH is a univariate model.")
  y = .data[[mv]]

  # Pull out inputs from the specials
  if(length(specials$season) > 1) stop("The `season()` special of `SMEAN()` should only be used once.")
  m = specials$season[[1]]

  glarma_model = glarma_call(y = y,
                               ic = ic,
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
    class = "GLARMA"
  )
}
