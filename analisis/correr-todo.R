# =============================================================================
# analisis/correr-todo.R  --  Equipo ST6Q98
# Se ejecuta desde la raiz del repositorio:  source("analisis/correr-todo.R")
#
# Organizacion:
#   Bloque 0   Configuracion: paquetes, funciones de R/, series
#   Bloque V   Verificacion de funciones propias (secciones 4 y 5)
#   Seccion 2  Series y pregunta de pronostico   (2a y 2c hechos; 2b va en el informe)
#   Seccion 3  Descripcion y patron               (3a, 3b y 3c hechos)
#   Seccion 4  Operador de rezago                 (4a, 4b y 4c hechos; 4d va en el informe)
#   Seccion 5  Dependencia y estacionariedad      (pendiente)
#   Seccion 6  Pronostico y evaluacion            (pendiente)
#   Cierre     resultados.rds, sessionInfo() y tiempo de ejecucion
# =============================================================================


# =============================================================================
# BLOQUE 0. CONFIGURACION
# =============================================================================

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

# Carga de las funciones propias, en orden (00-, 01-, ..., 05-)
for (archivo in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(archivo, encoding = "UTF-8")
}

dir.create("figs", showWarnings = FALSE)

guardar_fig <- function(grafico, nombre, ancho = 8, alto = 4.5) {
  ggplot2::ggsave(file.path("figs", paste0(nombre, ".png")), grafico,
                  width = ancho, height = alto, dpi = 150, bg = "white")
  invisible(nombre)
}

# Configuracion de las tres series.
# `esperado` son los valores que el enunciado declara para cada serie (T y fechas);
# la seccion 2(a) los contrasta contra lo que se lee de los csv.
series_cfg <- list(
  serie1 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-1.csv")),
    nombre = "Tiempo medio diario de espera en una estacion de bus",
    corto = "Espera diaria, estacion Los Almendros",
    variable = "Tiempo medio de espera",
    unidad = "segundos", frecuencia = "diaria",
    s_candidato = 7,
    # La serie termina el viernes 2025-10-03 y el operador decide ese dia para la
    # semana siguiente (lunes a domingo = horizontes h = 3, ..., 9).
    H = 9,
    esperado = list(T = 180, inicio = as.Date("2025-04-07"), fin = as.Date("2025-10-03"))
  ),
  serie2 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-2.csv")),
    nombre = "Pasajeros trimestrales de un aeropuerto regional",
    corto = "Pasajeros trimestrales, aeropuerto regional",
    variable = "Pasajeros",
    unidad = "miles de pasajeros", frecuencia = "trimestral",
    s_candidato = 4,
    H = 4,   # plan para los cuatro trimestres de 2026
    esperado = list(T = 120, inicio = as.Date("1996-01-01"), fin = as.Date("2025-10-01"))
  ),
  serie3 = list(
    archivo = file.path("datos", paste0(CODIGO, "-serie-3.csv")),
    nombre = "Viajes mensuales en una plataforma de movilidad",
    corto = "Viajes mensuales, plataforma de movilidad",
    variable = "Viajes",
    unidad = "miles de viajes", frecuencia = "mensual",
    s_candidato = 12,
    H = 6,   # enero a junio de 2026
    esperado = list(T = 144, inicio = as.Date("2014-01-01"), fin = as.Date("2025-12-01"))
  )
)

faltantes <- vapply(series_cfg, function(cfg) !file.exists(cfg$archivo), logical(1))
if (any(faltantes)) {
  stop("Faltan los archivos de datos: ",
       paste(vapply(series_cfg[faltantes], `[[`, "", "archivo"), collapse = ", "),
       ". Copie en datos/ los tres csv tal como se entregaron.")
}


# =============================================================================
# BLOQUE V. VERIFICACION DE FUNCIONES PROPIAS (secciones 4 y 5)
# Aqui, y solo aqui, se permiten diff(), acf(), pacf() y Box.test() de stats.
# =============================================================================

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

verif_polinomio <- local({
  r1  <- polinomio_rezago(c(-1))                 # 1 - B
  r2  <- polinomio_rezago(c(-2, 1))              # (1 - B)^2 = 1 - 2B + B^2
  r4  <- polinomio_rezago(c(0, 0, 0, -1))        # 1 - B^4
  r12 <- polinomio_rezago(c(rep(0, 11), -1))     # 1 - B^12
  ar1 <- polinomio_rezago(c(-0.6))               # 1 - 0.6 B
  evalpol <- function(coef, z) vapply(z, function(zz) sum(c(1, coef) * zz^(0:length(coef))), complex(1))
  list(
    raiz_unitaria = Mod(r1$raiz[1] - 1),
    raiz_doble    = max(Mod(r2$raiz - 1)),
    modulos_s4    = max(abs(r4$modulo - 1)),
    modulos_s12   = max(abs(r12$modulo - 1)),
    n_raices_s12  = nrow(r12),
    raiz_ar1      = abs(ar1$modulo - 1 / 0.6),
    residuo_s12   = max(Mod(evalpol(c(rep(0, 11), -1), r12$raiz)))
  )
})
stopifnot(
  verif_polinomio$raiz_unitaria < 1e-8,
  verif_polinomio$raiz_doble < 1e-6,          # raiz doble: polyroot pierde precision
  verif_polinomio$modulos_s4 < 1e-8,
  verif_polinomio$modulos_s12 < 1e-8,
  verif_polinomio$n_raices_s12 == 12,
  verif_polinomio$raiz_ar1 < 1e-8,
  verif_polinomio$residuo_s12 < 1e-8
)
message("Verificacion superada: polinomio_rezago().")


# =============================================================================
# SECCION 2. SERIES Y PREGUNTA DE PRONOSTICO
# =============================================================================

# ---- 2(a) Ficha de cada serie, verificada con codigo ------------------------
# La ficha sale de leer_serie() + ficha_serie(). Se contrasta (i) contra los
# valores que declara el enunciado y (ii) contra una lectura independiente del
# csv hecha con utils::read.csv, sin pasar por las funciones propias.

hacer_ficha <- function(id, cfg) {
  datos <- leer_serie(cfg$archivo)
  ficha <- ficha_serie(datos, cfg$nombre, cfg$unidad, cfg$frecuencia)
  
  crudo <- utils::read.csv(cfg$archivo)
  stopifnot(
    "T distinto del enunciado"            = ficha$T == cfg$esperado$T,
    "fecha inicial distinta del enunciado" = ficha$fecha_inicial == cfg$esperado$inicio,
    "fecha final distinta del enunciado"   = ficha$fecha_final == cfg$esperado$fin,
    "T distinto de la lectura independiente" = nrow(crudo) == ficha$T,
    "minimo distinto de la lectura independiente" = isTRUE(all.equal(ficha$minimo, min(crudo$valor))),
    "maximo distinto de la lectura independiente" = isTRUE(all.equal(ficha$maximo, max(crudo$valor))),
    "media distinta de la lectura independiente"  = isTRUE(all.equal(ficha$media, mean(crudo$valor)))
  )
  message(id, ": ficha verificada (T = ", ficha$T, ").")
  list(datos = datos, ficha = ficha)
}

lecturas <- purrr::imap(series_cfg, function(cfg, id) hacer_ficha(id, cfg))

datos_series <- lapply(lecturas, `[[`, "datos")     # tibbles fecha, valor
fichas <- dplyr::bind_rows(lapply(lecturas, `[[`, "ficha"), .id = "serie")
print(fichas)

# ---- 2(b) Pregunta de pronostico: va redactada en informe/informe.qmd --------
# ---- 2(c) Particion en tramo de estimacion y evaluacion de origen movil ------
# Misma regla para todos los metodos de una serie: n0 = floor(0.8 T).

particiones <- dplyr::bind_rows(
  purrr::imap(series_cfg, function(cfg, id) {
    d <- datos_series[[id]]
    particion(nrow(d), cfg$H, d$fecha)
  }),
  .id = "serie"
)
print(particiones)


# =============================================================================
# SECCION 3. DESCRIPCION Y PATRON
# =============================================================================

# ---- 3(a) Grafico en el tiempo y graficos de rezagos ------------------------
# Rezagos k elegidos para y_t contra y_{t-k}. La justificacion de cada k, con
# las autocorrelaciones muestrales, se escribe en informe/informe.qmd.
#   Serie 1 (diaria):     k = 1 (un dia), 4 (el de autocorrelacion mas negativa), 7 (semana)
#   Serie 2 (trimestral): k = 1 (trimestre vecino), 2 (fase opuesta), 4 (mismo trimestre del anio anterior)
#   Serie 3 (mensual):    k = 1 (mes vecino), 6 (medio anio), 12 (mismo mes del anio anterior)
rezagos_3a <- list(serie1 = c(1, 4, 7), serie2 = c(1, 2, 4), serie3 = c(1, 6, 12))

hacer_3a <- function(id, cfg, datos, rezagos) {
  y <- datos$valor
  n <- length(y)
  m <- min(floor(n / 4), 24)                       # rezagos del correlograma (3c)
  etiqueta_y <- paste0(cfg$variable, " (", cfg$unidad, ")")
  
  guardar_fig(graficar_serie(datos, cfg$corto, etiqueta_y), paste0(id, "-tiempo"))
  guardar_fig(
    grafico_rezagos(y, rezagos, paste0(cfg$corto, ": graficos de rezagos"), cfg$unidad),
    paste0(id, "-rezagos"), ancho = 9, alto = 3.6
  )
  message(id, ": figuras de 3(a) guardadas (tiempo y rezagos k = ",
          paste(rezagos, collapse = ", "), ").")
  list(rezagos = rezagos, m = m, acf = acf_muestral(y, m))
}

res3a <- purrr::imap(series_cfg, function(cfg, id) {
  hacer_3a(id, cfg, datos_series[[id]], rezagos_3a[[id]])
})

# ACF completa de cada serie (m rezagos): sirve para ELEGIR los rezagos de rezagos_3a.
# Si prefiere otros rezagos, cambie rezagos_3a mas arriba y vuelva a correr.
for (id in names(res3a)) {
  z <- res3a[[id]]
  cat("\n", id, ": autocorrelacion muestral r(k), k = 1, ..., ", z$m,
      "   (banda de ruido blanco +/- ", round(1.96 / sqrt(nrow(datos_series[[id]])), 3), ")\n", sep = "")
  print(round(setNames(z$acf, paste0("k", seq_along(z$acf))), 2))
  cat("Rezagos usados en los graficos:", paste(z$rezagos, collapse = ", "), "\n")
}

# ---- 3(a) Cifras de apoyo para la lectura del grafico, serie 1 --------------
d1 <- datos_series$serie1

resumen_s1 <- c(media = mean(d1$valor), sd = sd(d1$valor),
                min = min(d1$valor), max = max(d1$valor))
fecha_max_s1 <- d1$fecha[which.max(d1$valor)]

# Dias sobre la meta contractual de 330 s
dias_sobre_meta <- d1[d1$valor > 330, ]
n_sobre_meta <- nrow(dias_sobre_meta)

# Efecto del dia de la semana (1 = lunes, ..., 7 = domingo)
por_dia_semana <- d1 |>
  mutate(dia_sem = as.integer(format(fecha, "%u"))) |>
  group_by(dia_sem) |>
  summarise(media = mean(valor), sd = sd(valor), n = n())

# Efecto de inicio/fin de mes
por_borde_mes <- d1 |>
  mutate(borde_mes = as.integer(format(fecha, "%d")) %in% c(1:2, 29:31)) |>
  group_by(borde_mes) |>
  summarise(media = mean(valor), n = n())

print(round(resumen_s1, 2))
cat("Fecha del maximo:", format(fecha_max_s1), "\n")
cat("Dias sobre 330 s:", n_sobre_meta, "de", nrow(d1), "\n")
print(dias_sobre_meta)
print(por_dia_semana)
print(por_borde_mes)

# ---- 3(a) Cifras de apoyo para la lectura del grafico, serie 2 --------------
d2 <- datos_series$serie2 |>
  mutate(anio = as.integer(format(fecha, "%Y")),
         trim = (as.integer(format(fecha, "%m")) - 1) %/% 3 + 1)

media_anual_s2 <- d2 |> group_by(anio) |> summarise(media = mean(valor))
por_trimestre_s2 <- d2 |> group_by(trim) |> summarise(media = mean(valor), sd = sd(valor))

apoyo_s2 <- list(
  primero = d2$valor[1], ultimo = d2$valor[nrow(d2)],
  media_anual_ini = media_anual_s2$media[1],
  media_anual_fin = media_anual_s2$media[nrow(media_anual_s2)],
  crec_anual_medio = (media_anual_s2$media[nrow(media_anual_s2)] - media_anual_s2$media[1]) /
    (nrow(media_anual_s2) - 1),
  por_trimestre = por_trimestre_s2
)
print(apoyo_s2)

# ---- 3(b) Estacionalidad y descomposicion clasica a mano ----------------------
# DECISION DEL EQUIPO (edite estos valores si su analisis de la evidencia dice otra cosa):
#   estacional: la serie tiene estacionalidad con periodo s_candidato (TRUE/FALSE).
#   tipo:       version "aditiva" o "multiplicativa". Si la serie NO es estacional, es la
#               version con que se le quita la tendencia para mostrar la evidencia.
decision_3b <- list(
  serie1 = list(estacional = FALSE, tipo = "aditiva"),
  serie2 = list(estacional = TRUE,  tipo = "aditiva"),
  serie3 = list(estacional = FALSE, tipo = "multiplicativa")
)

hacer_3b <- function(id, cfg, datos, dec) {
  y <- datos$valor
  n <- length(y)
  s <- cfg$s_candidato
  m <- min(floor(n / 4), 24)
  etiqueta_y <- paste0(cfg$variable, " (", cfg$unidad, ")")
  est <- estacion_de_fecha(datos$fecha, s)
  cic <- ciclo_de_fecha(datos$fecha, s)
  etiquetas <- etiquetas_estacion(s)
  
  desc <- descomposicion_clasica(y, s, est, dec$tipo)
  pe <- prueba_estacionalidad(desc, que = cfg$corto)
  crit <- criterio_aditivo_multiplicativo(y, s, est, cic)
  
  # Version alternativa (para comparar cual de las dos describe mejor la serie)
  tipo_alt <- if (dec$tipo == "aditiva") "multiplicativa" else "aditiva"
  pe_alt <- prueba_estacionalidad(descomposicion_clasica(y, s, est, tipo_alt), que = cfg$corto)
  
  # Figuras para todas las series
  guardar_fig(grafico_estacional(datos, s, paste0(cfg$corto, ": grafico estacional (s = ", s, ")"),
                                 etiqueta_y, etiquetas), paste0(id, "-estacional"))
  guardar_fig(grafico_amplitud(crit$ciclos, paste0(cfg$corto, ": amplitud de cada ciclo contra su nivel"),
                               cfg$unidad), paste0(id, "-amplitud"))
  
  if (dec$estacional) {
    # Serie con estacionalidad: descomposicion completa (tendencia, indices, residuo)
    guardar_fig(grafico_descomposicion(desc, datos$fecha,
                                       paste0(cfg$corto, ": descomposicion clasica ", dec$tipo), cfg$unidad),
                paste0(id, "-descomposicion"), ancho = 9, alto = 7)
    guardar_fig(grafico_indices(desc, etiquetas, paste0(cfg$corto, ": indices estacionales normalizados")),
                paste0(id, "-indices"), ancho = 7, alto = 4)
  } else {
    # Serie sin estacionalidad: evidencia = grafico estacional de la serie sin tendencia
    sin_t <- tibble::tibble(fecha = datos$fecha, valor = desc$sin_tendencia)
    sin_t <- sin_t[!is.na(sin_t$valor), ]
    et_y <- if (dec$tipo == "aditiva") paste0("y - tendencia (", cfg$unidad, ")") else "y / tendencia (razon)"
    guardar_fig(grafico_estacional(sin_t, s, paste0(cfg$corto, ": grafico estacional sin tendencia (s = ", s, ")"),
                                   et_y, etiquetas), paste0(id, "-estacional-sin-tendencia"))
  }
  
  rango <- max(desc$indices) - min(desc$indices)
  list(decision = dec, s = s, n = n, desc = desc, pe = pe, pe_alt = pe_alt, tipo_alt = tipo_alt,
       crit = crit, etiquetas = etiquetas,
       efecto_pct = if (dec$tipo == "aditiva") 100 * rango / mean(y) else 100 * rango,
       r_s = acf_muestral(y, m)[s], banda = 1.96 / sqrt(n))
}

res3b <- purrr::imap(series_cfg, function(cfg, id) {
  hacer_3b(id, cfg, datos_series[[id]], decision_3b[[id]])
})

# Resultados en consola: la prueba F con sus seis elementos, los indices y el criterio de amplitud
for (id in names(res3b)) {
  z <- res3b[[id]]
  cat("\n", strrep("=", 70), "\n", id, ": s = ", z$s, ", version ", z$decision$tipo,
      ", decision: ", if (z$decision$estacional) "CON estacionalidad" else "SIN estacionalidad", "\n", sep = "")
  print(z$pe$prueba)
  cat(sprintf("R2 estacional = %.3f (version alternativa, %s: %.3f)\n",
              z$pe$R2_estacional, z$tipo_alt, z$pe_alt$R2_estacional))
  cat(sprintf("Rango de los indices = %.1f %% del nivel medio\n", z$efecto_pct))
  cat(sprintf("r(s = %d) = %.3f   (banda de ruido blanco +/- %.3f)\n", z$s, z$r_s, z$banda))
  cat("Indices estacionales normalizados:\n")
  print(round(setNames(z$desc$indices, z$etiquetas), 3))
  cat("Criterio amplitud vs nivel (ciclos completos):\n")
  print(round(z$crit$resumen, 3))
}

resumen3b <- dplyr::bind_rows(lapply(names(res3b), function(id) {
  z <- res3b[[id]]
  tibble::tibble(serie = id, s = z$s, version = z$decision$tipo,
                 decision = if (z$decision$estacional) "con estacionalidad" else "sin estacionalidad",
                 F = z$pe$F, gl1 = z$pe$gl1, gl2 = z$pe$gl2, p = z$pe$p,
                 R2_estacional = z$pe$R2_estacional, rango_indices_pct = z$efecto_pct,
                 r_s = z$r_s, banda_rb = z$banda,
                 cor_amplitud_nivel = unname(z$crit$resumen["cor_amplitud_nivel"]))
}))
cat("\nRESUMEN 3(b)\n")
print(as.data.frame(resumen3b), digits = 3)
# ---- 3(c) Correlograma de la serie sin transformar, m = min(floor(T/4), 24) 


compactar <- function(x) {
  if (length(x) == 0) return("ninguno")
  grupos <- split(x, cumsum(c(1, diff(x) != 1)))
  paste(vapply(grupos, function(g) {
    if (length(g) == 1) as.character(g) else paste0(min(g), "-", max(g))
  }, character(1)), collapse = ", ")
}

hacer_3c <- function(id, cfg, datos) {
  y <- datos$valor
  n <- length(y)
  m <- min(floor(n / 4), 24)
  s <- cfg$s_candidato
  r <- acf_muestral(y, m)
  
  # Las dos bandas del correlograma (misma formula que usa correlograma())
  banda_rb <- 1.96 / sqrt(n)                                        # ruido blanco, recta
  banda_bartlett <- 1.96 * sqrt((1 + 2 * c(0, cumsum(r^2))[seq_len(m)]) / n)   # secuencial de Bartlett
  fuera_rb <- which(abs(r) > banda_rb)
  fuera_bartlett <- which(abs(r) > banda_bartlett)
  dentro_rb <- which(abs(r) <= banda_rb)
  
  guardar_fig(
    correlograma(r, n, paste0(cfg$corto, ": ACF de la serie sin transformar"), "ACF"),
    paste0(id, "-acf-original")
  )
  message(id, ": ACF de la serie sin transformar guardada (m = ", m, ").")
  
  list(m = m, n = n, s = s, r = r, banda_rb = banda_rb, banda_bartlett = banda_bartlett,
       fuera_rb = fuera_rb, fuera_bartlett = fuera_bartlett,
       primer_dentro = if (length(dentro_rb) > 0) dentro_rb[1] else NA_integer_,
       n_positivos = sum(r > 0))
}

res3c <- purrr::imap(series_cfg, function(cfg, id) hacer_3c(id, cfg, datos_series[[id]]))

resumen3c <- dplyr::bind_rows(lapply(names(res3c), function(id) {
  z <- res3c[[id]]
  tibble::tibble(
    serie = id, T = z$n, m = z$m, banda_rb = z$banda_rb,
    r1 = z$r[1], r_s = z$r[z$s],
    rezagos_fuera_rb = compactar(z$fuera_rb), n_fuera_rb = length(z$fuera_rb),
    rezagos_fuera_bartlett = compactar(z$fuera_bartlett), n_fuera_bartlett = length(z$fuera_bartlett),
    primer_rezago_dentro = z$primer_dentro, positivos = z$n_positivos
  )
}))
cat("\nRESUMEN 3(c): lectura numerica del correlograma de la serie sin transformar\n")
print(as.data.frame(resumen3c), digits = 3)
for (id in names(res3c)) {
  cat("\n", id, ": r(k), k = 1, ..., ", res3c[[id]]$m, "\n", sep = "")
  print(round(setNames(res3c[[id]]$r, paste0("k", seq_len(res3c[[id]]$m))), 2))
}



# =============================================================================
# SECCION 4. OPERADOR DE REZAGO                                    (4a, 4b y 4c hechos)
# =============================================================================

# ---- 4(b) Transformaciones candidatas, polinomios y raices ------------------
# Se derivan del diagnostico de la seccion 3:
#   Serie 1: nivel estable        -> ninguna; (1-B) solo como contraste de sobrediferenciacion
#   Serie 2: tendencia + s = 4    -> (1-B), (1-B^4), (1-B)(1-B^4); log solo si 3(b) dio multiplicativa
#   Serie 3: tendencia, variabilidad creciente, sin estacionalidad -> con y sin log: (1-B), (1-B)^2
esc_s2 <- if (identical(decision_3b$serie2$tipo, "multiplicativa")) "log" else "ninguna"

cand_4b <- dplyr::bind_rows(
  tibble::tibble(serie = "serie1", escala = "ninguna", d = c(0, 1), D = 0, s = 1),
  tibble::tibble(serie = "serie2", escala = esc_s2, d = c(0, 1, 0, 1), D = c(0, 0, 1, 1),
                 s = series_cfg$serie2$s_candidato),
  tibble::tibble(serie = "serie3", escala = rep(c("ninguna", "log"), each = 3),
                 d = rep(c(0, 1, 2), 2), D = 0, s = 1)
) |>
  dplyr::mutate(
    operador     = purrr::pmap_chr(list(d, D, s), etiqueta_operador),
    operador_txt = purrr::pmap_chr(list(d, D, s), etiqueta_operador, tex = FALSE),
    coef         = purrr::pmap(list(d, D, s), polinomio_diferencias),
    grado        = purrr::map_int(coef, ~ length(.x) - 1L),
    polinomio    = purrr::map_chr(coef, polinomio_texto)
  )

raices_4b <- purrr::pmap(cand_4b, function(serie, escala, d, D, s, operador, operador_txt,
                                           coef, grado, ...) {
  if (grado == 0L) return(NULL)                      # c(B) = 1 no tiene raices
  r <- polinomio_rezago(coef[-1])
  tibble::tibble(serie = serie, escala = escala, operador = operador,
                 candidata = paste0(if (escala == "log") "log, " else "", operador_txt),
                 real = r$real, imaginaria = r$imaginaria,
                 modulo = r$modulo, periodo = r$periodo)
}) |> dplyr::bind_rows()

# Verificaciones: una raiz por grado y todas sobre el circulo unidad
stopifnot(nrow(raices_4b) == sum(cand_4b$grado),
          all(abs(raices_4b$modulo - 1) < 1e-5))

for (id in names(series_cfg)) {
  guardar_fig(
    grafico_raices(dplyr::filter(raices_4b, serie == id),
                   paste0(series_cfg[[id]]$corto, ": raices de los operadores candidatos")),
    paste0(id, "-raices"), ancho = 9, alto = 4)
}
print(as.data.frame(cand_4b[, c("serie", "escala", "operador_txt", "polinomio", "grado")]))



# ---- 4(c) Varianza y r1 de cada candidata, eleccion y verificacion de la media ----
# Para cada candidata de 4(b): media, varianza muestral y r1 de la serie transformada,
# y la media que el operador debe dejar (nivel, pendiente, s * pendiente o 0).
tabla_4c <- purrr::pmap_dfr(
  dplyr::select(cand_4b, serie, escala, d, D, s, operador_txt),
  function(serie, escala, d, D, s, operador_txt) {
    dplyr::bind_cols(
      tibble::tibble(serie = serie, escala = escala, d = d, D = D, s = s,
                     candidata = paste0(if (escala == "log") "log, " else "", operador_txt)),
      resumen_transformacion(datos_series[[serie]]$valor, escala, d, D, s))
  }) |>
  dplyr::group_by(serie, escala) |>
  dplyr::mutate(razon_var = varianza / varianza[d == 0 & D == 0]) |>   # frente a la misma escala sin diferenciar
  dplyr::ungroup()
cat("\nTABLA 4(c): varianza y r1 de cada candidata\n")
print(as.data.frame(tabla_4c), digits = 3)

# Identidad Var(nabla Y) = 2 gamma(0) - 2 gamma(1), para las candidatas con (1-B)
identidad_4c <- purrr::pmap_dfr(dplyr::filter(cand_4b, d == 1, D == 0), function(serie, escala, ...) {
  y <- datos_series[[serie]]$valor
  if (escala == "log") y <- log(y)
  dplyr::bind_cols(tibble::tibble(serie = serie, escala = escala), identidad_varianza(y))
})
cat("\nIdentidad Var(nabla Y) = 2 gamma(0) - 2 gamma(1)\n")
print(as.data.frame(identidad_4c), digits = 3)

# Que le pasa a r1 al diferenciar una serie que ya era estacionaria (ruido blanco simulado)
set.seed(3009297)
ruido <- rnorm(500)
sobrediferencia <- tibble::tibble(
  r1_ruido  = acf_muestral(ruido, 1),
  r1_difer  = acf_muestral(diferenciar(ruido, d = 1), 1),
  razon_var = stats::var(diferenciar(ruido, d = 1)) / stats::var(ruido))
cat("\nSobrediferenciacion de un ruido blanco simulado\n")
print(sobrediferencia)

# ELECCION de la transformacion de cada serie (edite si cambia su lectura de tabla_4c)
elegida_4c <- list(
  serie1 = list(escala = "ninguna", d = 0, D = 0, s = 1),
  serie2 = list(escala = "ninguna", d = 1, D = 1, s = 4),
  serie3 = list(escala = "log",     d = 1, D = 0, s = 1)
)
z_elegida <- purrr::imap(elegida_4c, function(e, id)
  serie_transformada(datos_series[[id]]$valor, e$escala, e$d, e$D, e$s))   # entra a la seccion 5
media_elegida <- dplyr::bind_rows(lapply(names(elegida_4c), function(id) {
  e <- elegida_4c[[id]]
  dplyr::bind_cols(tibble::tibble(serie = id),
                   resumen_transformacion(datos_series[[id]]$valor, e$escala, e$d, e$D, e$s))
}))
cat("\nMedia de la serie transformada elegida frente a la esperada\n")
print(as.data.frame(media_elegida), digits = 3)



# =============================================================================
# SECCION 5. DEPENDENCIA Y ESTACIONARIEDAD

#
# Trabaja sobre la serie resultante de la transformacion elegida en la seccion 4
# (z_elegida): serie 1 sin transformar, serie 2 con (1-B)(1-B^4), serie 3 con log y (1-B).
# La verificacion de pacf_dl() contra pacf() ya esta en el bloque V (verificacion$pacf_dl).
# =============================================================================

# DECISION DEL EQUIPO: m de la prueba de Ljung-Box (5b). Propuesta inicial: dos ciclos
# estacionales, m = 2 * s_candidato. Justifique su eleccion en el informe y cambie aqui
# si prefiere otro valor. Como la prueba se aplica a la serie transformada y no a los
# residuos de un modelo ajustado, no se restan grados de libertad (g = 0).
m_lb <- c(serie1 = 14, serie2 = 8, serie3 = 24)

hacer_5 <- function(id, cfg, z_cruda, fechas, m_lb, transformacion) {
  z <- z_cruda[!is.na(z_cruda)]
  n <- length(z)
  fechas_z <- utils::tail(fechas, n)           # la transformacion pierde las primeras observaciones
  m <- min(floor(n / 4), 24)
  
  # ---- 5(a) ACF y PACF muestrales con sus bandas --------------------------------
  r <- acf_muestral(z, m)
  p <- pacf_dl(r)
  banda_rb <- 1.96 / sqrt(n)
  banda_bartlett <- 1.96 * sqrt((1 + 2 * c(0, cumsum(r^2))[seq_len(m)]) / n)
  g_acf <- correlograma(r, n, paste0(cfg$corto, ": ACF"), "ACF")
  g_pacf <- correlograma(p, n, paste0(cfg$corto, ": PACF"), "PACF")
  guardar_fig(
    (g_acf | g_pacf) + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom"),
    paste0(id, "-acf-pacf"), ancho = 11, alto = 4.8)
  
  # ---- 5(b) Ljung-Box con seis elementos, y sensibilidad al m -------------------
  lb <- ljung_box(z, m_lb, g = 0,
                  que = paste0("la serie transformada (", cfg$corto, ")"))
  m_sens <- sort(unique(c(5, 10, cfg$s_candidato, 2 * cfg$s_candidato, 24)))
  m_sens <- m_sens[m_sens < n]
  sensibilidad <- dplyr::bind_rows(lapply(m_sens, function(mm) {
    t <- ljung_box(z, mm, g = 0)
    tibble::tibble(m = mm, Q = t$observado, gl = mm, critico = t$critico, p = t$p_valor)
  }))
  
  # ---- 5(c) Cuatro bloques consecutivos, media y varianza por bloque ------------
  bloque <- ceiling(4 * seq_len(n) / n)
  bloques <- dplyr::bind_rows(lapply(1:4, function(b) {
    zb <- z[bloque == b]
    fb <- fechas_z[bloque == b]
    tibble::tibble(bloque = b, n = length(zb), desde = fb[1], hasta = fb[length(fb)],
                   media = mean(zb), varianza = stats::var(zb), desv = stats::sd(zb))
  }))
  g_bloques <- ggplot2::ggplot(data.frame(fecha = fechas_z, z = z), ggplot2::aes(x = fecha, y = z)) +
    ggplot2::geom_line(colour = "grey35", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = as.numeric(bloques$desde[-1]), linetype = "dashed", colour = "grey60") +
    ggplot2::geom_segment(data = bloques,
                          ggplot2::aes(x = desde, xend = hasta, y = media, yend = media),
                          colour = "firebrick", linewidth = 1.1, inherit.aes = FALSE) +
    ggplot2::labs(title = paste0(cfg$corto, ": serie transformada y media por bloque"),
                  x = "Fecha", y = paste0("Serie transformada: ", transformacion),
                  caption = paste0("n = ", n, " observaciones, cuatro bloques consecutivos; ",
                                   "linea roja = media de cada bloque")) +
    ggplot2::theme_minimal(base_size = 11)
  guardar_fig(g_bloques, paste0(id, "-bloques"))
  
  # Media con dos errores estandar: el ingenuo y el de largo plazo (nucleo de Bartlett)
  L <- floor(4 * (n / 100)^(2 / 9))
  media <- mean(z)
  gamma0 <- autocovarianzas(z, 0)
  var_lp <- varianza_largo_plazo(z, L)
  se_ing <- stats::sd(z) / sqrt(n)
  se_lp <- sqrt(var_lp / n)
  
  message(id, ": seccion 5 calculada (n = ", n, ", m = ", m, ", m de Ljung-Box = ", m_lb, ").")
  list(transformacion = transformacion, n = n, m = m, r = r, pacf = p,
       banda_rb = banda_rb, banda_bartlett = banda_bartlett,
       fuera_acf = which(abs(r) > banda_rb), fuera_acf_bartlett = which(abs(r) > banda_bartlett),
       fuera_pacf = which(abs(p) > banda_rb),
       m_lb = m_lb, lb = lb, sensibilidad = sensibilidad,
       bloques = bloques, L = L, media = media, gamma0 = gamma0, var_lp = var_lp,
       se_ing = se_ing, se_lp = se_lp,
       ic_ing = media + c(-1.96, 1.96) * se_ing, ic_lp = media + c(-1.96, 1.96) * se_lp)
}

transformacion_txt <- vapply(names(elegida_4c), function(id) {
  e <- elegida_4c[[id]]
  paste0(if (e$escala == "log") "log, " else "", etiqueta_operador(e$d, e$D, e$s, tex = FALSE))
}, character(1))

res5 <- purrr::imap(series_cfg, function(cfg, id) {
  hacer_5(id, cfg, z_elegida[[id]], datos_series[[id]]$fecha, m_lb[[id]], transformacion_txt[[id]])
})

resumen5 <- dplyr::bind_rows(lapply(names(res5), function(id) {
  z <- res5[[id]]
  tibble::tibble(
    serie = id, transformacion = z$transformacion, n = z$n, m = z$m,
    banda_rb = z$banda_rb,
    acf_fuera_rb = compactar(z$fuera_acf), acf_fuera_bartlett = compactar(z$fuera_acf_bartlett),
    pacf_fuera_rb = compactar(z$fuera_pacf),
    m_lb = z$m_lb, Q = z$lb$observado, p_lb = z$lb$p_valor,
    media = z$media, L = z$L, se_ing = z$se_ing, se_lp = z$se_lp,
    razon_se = z$se_lp / z$se_ing,
    ic_ing_inf = z$ic_ing[1], ic_ing_sup = z$ic_ing[2],
    ic_lp_inf = z$ic_lp[1], ic_lp_sup = z$ic_lp[2]
  )
}))

cat("\nRESUMEN SECCION 5\n")
print(as.data.frame(resumen5), digits = 4)
for (id in names(res5)) {
  cat("\n", strrep("=", 70), "\n", id, " (", res5[[id]]$transformacion, ")\n", sep = "")
  print(res5[[id]]$lb)
  cat("\nLjung-Box segun m:\n")
  print(as.data.frame(res5[[id]]$sensibilidad), digits = 4)
  cat("\nMedia por bloque:\n")
  print(as.data.frame(res5[[id]]$bloques), digits = 4)
}


# =============================================================================
# SECCION 6. PRONOSTICO Y EVALUACION                               (pendiente)
# =============================================================================


# =============================================================================
# CIERRE
# =============================================================================

# Resultados que lee informe/informe.qmd (el informe no recalcula el analisis)
saveRDS(
  list(series_cfg = series_cfg, datos_series = datos_series,
       fichas = fichas, particiones = particiones,
       res3a = res3a,
       apoyo_s1 = list(resumen = resumen_s1, fecha_max = fecha_max_s1,
                       n_sobre_meta = n_sobre_meta,
                       por_dia_semana = por_dia_semana,
                       por_borde_mes = por_borde_mes),
       apoyo_s2 = apoyo_s2,
       res3b = res3b, resumen3b = resumen3b,
       res3c = res3c, resumen3c = resumen3c,
       verificacion = verificacion, verif_polinomio = verif_polinomio,
       cand_4b = cand_4b, raices_4b = raices_4b,
       tabla_4c = tabla_4c, identidad_4c = identidad_4c,
       sobrediferencia = sobrediferencia, elegida_4c = elegida_4c,
       media_elegida = media_elegida,
       z_elegida = z_elegida,
       res5 = res5, resumen5 = resumen5),
  file.path("figs", "resultados.rds")
)

writeLines(capture.output(sessionInfo()), "sesion-info.txt")
tiempo_ejecucion <- difftime(Sys.time(), inicio_ejecucion, units = "secs")
message("correr-todo.R termino en ", round(as.numeric(tiempo_ejecucion), 1), " segundos.")