% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/INGARCH.R
\name{INGARCH}
\alias{INGARCH}
\title{Estimate a INGARCH model}
\usage{
INGARCH(
  formula,
  ic = c("aic", "bic", "qic"),
  link = c("identity", "log"),
  distr = c("poisson", "nbinom"),
  algorithm = c("naive_search", "arma_to_ingarch"),
  trace = FALSE
)
}
\arguments{
\item{formula}{Model specification (see "Specials" section).}

\item{ic}{Character, can be 'aic' 'bic' or 'qic'. The information criterion used in selecting the model.}

\item{link}{Character, can be 'identity' or 'log' The link function used for the generalized model}

\item{distr}{Character, can be 'poisson' or 'nbinom'. The probabilty distribution used for the generalized model}

\item{algorithm}{Character, specifies the automatic parameter selection algorithm. Can be 'naive_search' or 'arma_to_ingarch'.
If 'naive_search' is selected, a search in a 4x4 matrix parameter space is performed, where the model to minimize the criterion value is selected.
If 'arma_to_ingarch' is selected, uses an auto_arma as the starting point for the selection algorithm.
The ‘arma_to_ingarch’ is the only one to perform a seasonal adjustment}

\item{trace}{Logical. If the automatic parameter algorithm is runnig, print the path to the best model estimation}
}
\value{
A model specification.
}
\description{
Estimate Integer-valued Generalized Autoregressive Conditional Heteroscedasticity model
with Poisson or Negative Binomial distribution.
Also is provide a automatic parameter algorithm selection for the Autorregressive and Moving Avarege params
}
\section{Specials}{


\subsection{pq}{
pq defines the non-seasonal autoregressive and moving avarages terms,
it can be define by the user,
or if it's omited, the automatic parameter selection algorithm is trigered
The automatic parameter selection algorithm gonna fit the best model based on the information criterion
}

\subsection{PQ}{
PQ defines the seasonal autoregressive and moving avarages terms,
it can be define by the user,
or if it's omited, the automatic parameter selection algorithm is trigered (only for 'arma_to_ingarch' algorithm)
The automatic parameter selection algorithm gonna fit the best model based on the information criterion
}


\subsection{xreg}{
Exogenous regressors can be included in a INGARCH model without explicitly using the `xreg()` special.
Common exogenous regressor specials as specified in [`common_xregs`] can also be used.
These regressors are handled using [stats::model.frame()],
and so interactions and other functionality behaves similarly to [stats::lm()].

The inclusion of a constant in the model follows the similar rules to [`stats::lm()`],
where including `1` will add a constant and `0` or `-1` will remove the constant.
If left out, the inclusion of a constant will be determined by minimising `ic`.
}

If a xreg is provided, the model forecast is not avaliable

\preformatted{
xreg(..., fixed = list())
}

\tabular{ll}{
  `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)\cr
  `fixed`    \tab A named list of fixed parameters for coefficients. The names identify the coefficient, and should match the name of the regressor. For example, `fixed = list(constant = 20)`.
}
}

\examples{
\donttest{
# Manual INGARCH specification
tsibbledata::aus_production |>
  fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1)))

# Automatic INGARCH specification
 tsibbledata::aus_production |>
fabletools::model(auto_ing_naive =
                    INGARCH(Beer,
                            ic = 'aic',
                            trace = TRUE,
                           algorithm = 'naive_search'),
                  auto_ing_arm_ing =
                    INGARCH(Beer,
                            ic = 'aic',
                            trace = TRUE,
                            algorithm = 'arma_to_ingarch'))
                            }


}
