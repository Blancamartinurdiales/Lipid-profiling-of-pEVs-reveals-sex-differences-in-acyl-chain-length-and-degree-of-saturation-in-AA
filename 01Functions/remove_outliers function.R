
# FUNCION remove_outliers() PARA ELIMINAR OUTLIERS

library(stringr)
library(plotly)
library(dplyr)
library(tidyverse)
library(data.table)
library(DT)

# Cargamos la funcion para eliminar outliers:
# https://stackoverflow.com/questions/4787332/how-to-remove-outliers-from-a-dataset

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Para esta función necetitamos un data frame por cada grupo experimental y la transversa (t())
