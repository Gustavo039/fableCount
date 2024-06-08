# teste = fableCount::covid_sp |>
#   tibble::as_tibble() |>
#   dplyr::select(newCases) |>
#   dplyr::mutate(lag1 = dplyr::lead(newCases),
#                 lag2 = dplyr::lead(newCases, 2)) |>
#   dplyr::mutate(
#     across(c(1:3),
#     ~dplyr::case_when(is.na(.) ~ 0, .default = .)
#     )
#   )
library(fpp3)

model = arima.sim(model = list(order = c(3, 0, 1), ar = c(0.4, 0.3, 0.2), ma = c(0.7)), n = 100)

teste_ar = model |>
  tibble::as_tibble() |>
  dplyr::select(x) |>
  dplyr::mutate(lag1 = dplyr::lead(x),
                lag2 = dplyr::lead(x, 2),
                lag3 = dplyr::lead(x, 3),
                lag4 = dplyr::lead(x, 4),
                lag5 = dplyr::lead(x, 5)) |>
  dplyr::mutate(
    across(c(1:6),
    ~dplyr::case_when(is.na(.) ~ 0, .default = .)
    )
  )


lasso_ar_select = glmnet::glmnet(y = teste_ar[,1] |> as.matrix(),
               x = teste_ar[,2:6] |> as.matrix(),
               alpha = 1, intercept = F)


lags_coef = lasso_ar_select |>
  coef(s = 0.1) |>
  as.vector()

lags_coef = which(lags_coef == 0)[2] - 2

arima_model = arima(model, order = c(lags_coef, 0, 6 - lags_coef))

teste_ma = arima_model$residuals |>
  tibble::as_tibble() |>
  dplyr::select(x) |>
  dplyr::mutate(ma1 = dplyr::lead(x),
                ma2 = dplyr::lead(x, 2),
                ma3 = dplyr::lead(x, 3),
                ma4 = dplyr::lead(x, 4),
                ma5 = dplyr::lead(x, 5)) |>
  dplyr::mutate(
    across(c(1:6),
           ~dplyr::case_when(is.na(.) ~ 0, .default = .)
    )
  )


teste_ar_ma = teste_ar |>
  dplyr::bind_cols(teste_ma |>
                     select(-x))

lasso_ar_ma_select = glmnet::glmnet(y = teste_ar_ma[,1] |>
                                      as.matrix(),
                                    x = teste_ar_ma[,2:11] |>
                                      as.matrix(),
                                    alpha = 1, intercept = F)


lasso_ar_ma_select |>
  coef(s = 0.1)

library(fpp3)
library(patchwork)


lasso_ma_select = glmnet::glmnet(y = teste_ar[,1] |> as.matrix(),
                                 x = teste[,2:6] |> as.matrix(),
                                 alpha = 1, intercept = F)


tscount::rdistr(100, 3, distr = 'nbinom', 5) |>
  {\(x) paste(x |> mean(),
              x |> var())}()


modelo_teste = tscount::tsglm(tscount::rdistr(100, 3, distr = 'nbinom', 5))

######
df_isolamento = readxl::read_xlsx('D:/TCC/dados/isolamento_sp.xlsx') |>
  janitor::row_to_names(1) |>
  filter(Município1 == 'SÃO PAULO') |>
  select(-c(2:4)) |>
  rename(mun = 1) |>
  tidyr::pivot_longer(-1) |>
  mutate(name =
           name |>
           lubridate::dmy(),
         value =
           value |>
           as.numeric())

df_isolamento |>
  tsibble::as_tsibble(index = name) |>
  tsibble::fill_gaps()


######

fig1 = tsibbledata::global_economy |>
  dplyr::filter(Country == 'Brazil') |>
  select(GDP) |>
  ggplot(aes(x = Year, y = GDP))+
  geom_line(col = 'orange', linewidth = 1.15) +
  labs(title = 'Evolução do PIB Brasileiro', x = 'Ano', y = 'PIB') +
  theme_minimal()

fig2 = df_isolamento |>
  tsibble::as_tsibble(index = name) |>
  tsibble::fill_gaps() |>
  ggplot(aes(x = name, y = value))+
  geom_line(col = 'black', linewidth = 1.15) +
  labs(title = 'Nível de Isolamento Social em SP', x = 'Ano-Mês', y = 'Proporção Isolado') +
  theme_minimal()

devtools::load_all()
fig3 = fableCount::covid_sp |>
  ggplot(aes(x = date, y = newDeaths))+
  geom_line(col = 'red', linewidth = 1.15) +
  labs(title = 'Mortes por COVID em SP', x = 'Ano-Semana' , y = 'Número de Mortos') +
  theme_minimal()

fig4 = tsibbledata::global_economy |>
  dplyr::filter(Country == 'Brazil') |>
  tibble::as_tibble() |>
  select(Imports) |>
  mutate(data = seq(from = as.POSIXct("2024-09-20 15:10:00"),
                    to = as.POSIXct("2024-09-20 15:10:57"),
                    by = "1 sec"),
         Imports = log(Imports*10000)) |>
  as_tsibble(index = data) |>
  ggplot(aes(x = data, y = Imports))+
  geom_line(col = 'blue', linewidth = 1.15) +
  labs(title = 'Variação No Preço das Ações CSNA3', x = 'Segundos', y = 'Preço da Ação') +
  theme_minimal()

(fig1 + fig2)/(fig3 + fig4)

###################


## overdiperson and underdisperson tests
library(AER)
library(tidyverse)
# Simular séries

##Simulando Poisson

simula_poisson = function(lambda = 5){
  dados_poisson = sapply(c(50, 100, 200, 500, 1000),
                         function(i){
                           tscount::rdistr(i, meanvalue = lambda, distr = 'poisson')
                         })

  # Testar modelos

  p_valores = dados_poisson |>
    lapply(function(i){
      dados_iteração = i |>
        tibble::as_tibble()
      glm(value ~ ., family = poisson(), data = dados_iteração) |>
        AER::dispersiontest() |>
        {\(x) x$p.value}()
    })

  return(p_valores |> unlist())

}



replicate(100, expr = simula_poisson()) |>
  apply(1, function(i){
    sum(i <= 0.05)/length(i)
  })


simula_nbinom = function(lambda = 5){
  dados_nbinom = sapply(c(50, 100, 200, 500, 1000),
                        function(i){
                          tscount::rdistr(i, meanvalue = lambda, distr = 'nbinom', distrcoefs = lambda*2)
                        })

  # Testar modelos

  p_valores = dados_nbinom |>
    lapply(function(i){
      dados_iteração = i |>
        tibble::as_tibble()
      glm(value ~ ., family = poisson(), data = dados_iteração) |>
        AER::dispersiontest() |>
        {\(x) x$p.value}()
    })

  return(p_valores |> unlist())

}

replicate(100, expr = simula_nbinom()) |>
  apply(1, function(i){
    sum(i <= 0.05)/length(i)
  })

##Simulando Binomial Negativa
dados_nbinom = sapply(c(50, 100, 200, 500, 1000),
                      function(i){
                        tscount::rdistr(i, meanvalue = lambda, distr = 'nbinom', distrcoefs = lambda^2)
                      })


teste = dados_poisson[[1]] |>
  tibble::as_tibble()

glm(value ~ ., family = poisson(), data = teste) |>
  AER::dispersiontest() |>
  {\(x) x$p.value}()





simula_tscount = function(n, distr, distr_coef){
  dados_poisson = tscount::tsglm.sim(n,
                                     param = list(intercept = 1, past_obs = c(0.5, 0.4),
                                                  past_mean = NULL, xreg = NULL),
                                     model = list(past_obs = c(1,2), past_mean = NULL, external = FALSE),
                                     xreg = NULL, link = 'identity', distr = distr, distrcoefs = distr_coef)

  dados_poisson |>
    tibble::as_tibble() |>
    dplyr::select(1) |>
    mutate(date = as.Date(lubridate::dmy('01-01-1923'):lubridate::dmy('10-04-1923'))) |>
    tsibble::tsibble(index = date) |>
    model(teste = ARIMA(ts))

  aic_poi = tscount::tsglm(dados_poisson$ts |>
                             as.vector(),
                           model = list(past_obs = 2, past_mean = NULL), distr = 'poisson') |>
    summary() |>
    {\(x) x$AIC}()
  aic_nbin = tscount::tsglm(dados_poisson$ts |>
                              as.vector(),
                            model = list(past_obs = 2, past_mean = NULL), distr = 'nbinom') |>
    summary() |>
    {\(x) x$AIC}()

  if(aic_poi < aic_nbin) return(1)
  else return(0)
}

replicate(100, simula_tscount(100, 'poisson')) |> mean()

