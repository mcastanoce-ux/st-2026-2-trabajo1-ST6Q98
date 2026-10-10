
if (!file.exists(file.path("R", "00-lectura.R"))) {
  stop("Ejecute este script desde la raiz del repositorio st-2026-2-trabajo1-ST6Q98.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

inicio_ejecucion <- Sys.time()
CODIGO <- "ST6Q98"

for (archivo in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(archivo, encoding = "UTF-8")
}
dir.create("figs", showWarnings = FALSE)

guardar_fig <- function(grafico, nombre, ancho = 8, alto = 4.5) {
  ggplot2::ggsave(file.path("figs", paste0(nombre, ".png")), grafico,
                  width = ancho, height = alto, dpi = 150, bg = "white")
  invisible(nombre)
}


series_cfg <- list(
  serie1 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-1.csv")),
    nombre = "Tiempo medio diario de espera en una estacion de bus",
    corto = "Espera diaria, estacion Los Almendros",
    unidad = "segundos", frecuencia = "diaria",
    s_candidato = 7,
    H = 7    # decision de cada viernes para la semana siguiente (revisar en 2b)
  ),
  serie2 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-2.csv")),
    nombre = "Pasajeros trimestrales de un aeropuerto regional",
    corto = "Pasajeros trimestrales, aeropuerto regional",
    unidad = "miles de pasajeros", frecuencia = "trimestral",
    s_candidato = 4,
    H = 4    # plan para los cuatro trimestres siguientes
  ),
  serie3 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-3.csv")),
    nombre = "Viajes mensuales en una plataforma de movilidad",
    corto = "Viajes mensuales, plataforma de movilidad",
    unidad = "miles de viajes", frecuencia = "mensual",
    s_candidato = 12,
    H = 6    # proximos seis meses
  )
)

faltantes <- vapply(series_cfg, function(cfg) !file.exists(cfg$archivo), logical(1))
if (any(faltantes)) {
  stop("Faltan los archivos de datos: ",
       paste(vapply(series_cfg[faltantes], `[[`, "", "archivo"), collapse = ", "),
       ". Copie en datos/ los tres csv tal como se entregaron.")
}


verificacion <- local({
  set.seed(3009297)
  x <- cumsum(rnorm(150)) + 10 * sin(2 * pi * (1:150) / 12)
  u <- rnorm(200)
  e <- numeric(200)
  e[1] <- u[1]
  for (t in 2:200) e[t] <- 0.6 * e[t - 1] + u[t]     # AR(1) simulado a mano
  m <- 24
  list(
    diferenciar = c(
      d1 = max(abs(diferenciar(x, d = 1) - diff(x))),
      d2 = max(abs(diferenciar(x, d = 2) - diff(x, differences = 2))),
      D1_s12 = max(abs(diferenciar(x, D = 1, s = 12) - diff(x, lag = 12))),
      d1_D1_s12 = max(abs(diferenciar(x, d = 1, D = 1, s = 12) - diff(diff(x, lag = 12))))
    ),
    acf_muestral = max(abs(acf_muestral(e, m) - stats::acf(e, lag.max = m, plot = FALSE)$acf[-1])),
    pacf_dl = max(abs(pacf_dl(acf_muestral(e, m)) - as.numeric(stats::pacf(e, lag.max = m, plot = FALSE)$acf))),
    ljung_box = abs(ljung_box(e, m)$observado -
                      unname(stats::Box.test(e, lag = m, type = "Ljung-Box")$statistic))
  )
})
stopifnot(
  all(verificacion$diferenciar < 1e-10),
  verificacion$acf_muestral < 1e-10,
  verificacion$pacf_dl < 1e-10,
  verificacion$ljung_box < 1e-8
)
message("Verificaciones superadas: diferenciar(), acf_muestral(), pacf_dl(), ljung_box().")


procesar_serie <- function(id, cfg) {
  datos <- leer_serie(cfg$archivo)
  y <- datos$valor
  T <- length(y)
  m <- min(floor(T / 4), 24)

  # Seccion 2: ficha y particion
  ficha <- ficha_serie(datos, cfg$nombre, cfg$unidad, cfg$frecuencia)
  part <- particion(T, cfg$H, datos$fecha)

  # Seccion 3: descripcion
  g_tiempo <- graficar_serie(datos, cfg$corto, cfg$unidad)
  guardar_fig(g_tiempo, paste0(id, "-tiempo"))

  rezagos <- sort(unique(c(1, 2, cfg$s_candidato)))     # revisar eleccion en 3a
  g_rez <- grafico_rezagos(y, rezagos, paste0(cfg$corto, ": graficos de rezagos"), cfg$unidad)
  guardar_fig(g_rez, paste0(id, "-rezagos"), ancho = 9, alto = 3.6)

  g_est <- grafico_estacional(datos, cfg$s_candidato,
                              paste0(cfg$corto, ": grafico estacional (s = ", cfg$s_candidato, ")"),
                              cfg$unidad)
  guardar_fig(g_est, paste0(id, "-estacional"))

  r <- acf_muestral(y, m)
  g_acf <- correlograma(r, T, paste0(cfg$corto, ": ACF de la serie sin transformar"), "ACF")
  guardar_fig(g_acf, paste0(id, "-acf-original"))

  # TODO seccion 3b: descomposicion clasica a mano (si hay estacionalidad)
  # TODO seccion 4:  transformaciones candidatas, polinomios y raices
  # TODO seccion 5:  ACF/PACF, Ljung-Box y medias por bloques de la serie elegida
  # TODO seccion 6:  optimizar(), origen_movil(), medidas(), validar_errores()

  list(cfg = cfg, datos = datos, T = T, m = m, ficha = ficha, particion = part, acf = r)
}

resultados <- purrr::imap(series_cfg, function(cfg, id) procesar_serie(id, cfg))

fichas <- dplyr::bind_rows(lapply(resultados, `[[`, "ficha"), .id = "serie")
particiones <- dplyr::bind_rows(lapply(resultados, `[[`, "particion"), .id = "serie")
print(fichas)
print(particiones)

#CIERRE.
writeLines(capture.output(sessionInfo()), "sesion-info.txt")
tiempo_ejecucion <- difftime(Sys.time(), inicio_ejecucion, units = "secs")
message("correr-todo.R termino en ", round(as.numeric(tiempo_ejecucion), 1), " segundos.")
