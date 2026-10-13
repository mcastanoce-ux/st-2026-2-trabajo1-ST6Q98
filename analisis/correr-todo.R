
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
# BLOQUE V6. VERIFICACION DE LOS METODOS DE PRONOSTICO (seccion 6a)
# Cada ajustar_*() se contrasta contra una cifra de las notas (clases 3, 4, 5 y 9)
# o contra un calculo a mano. lm() solo aparece aqui, como verificacion.
# =============================================================================

verif6 <- local({
  ok <- function(x, objetivo, tol) isTRUE(max(abs(x - objetivo)) < tol)
  err <- function(ajuste, y, desde = 1) {           # errores a un paso dentro de muestra
    e <- (y - ajuste$yhat)[desde:length(y)]
    e[!is.na(e)]
  }
  mse <- function(e) mean(e^2); mad <- function(e) mean(abs(e))
  
  # Datos de las notas
  x_ej1    <- c(4, 7, 6, 9, 12, 10)                                     # clase 9, ejemplo 1
  precio   <- c(18.43, 19.98, 19.51, 20.63, 19.78, 21.25, 21.18, 22.14)  # clase 3, diap. 45 y 37
  creditos <- c(133, 155, 165, 171, 194, 231, 274, 312, 313, 333, 343)  # clase 4, diap. 66
  ventas   <- c(1.5, 4.3, 8.9, 12.5, 1.8, 4.5, 9.0, 12.1,               # clase 4, diap. 34
                1.2, 4.6, 8.1, 12.9, 1.9, 5.1, 8.6, 13.2)
  
  filas <- list()
  agregar <- function(metodo, contra, cifra, obtenido, pasa) {
    filas[[length(filas) + 1]] <<- tibble::tibble(
      metodo = metodo, contra = contra, cifra_notas = cifra,
      obtenido = obtenido, verificado = pasa)
  }
  
  # 1. Ingenuo: su yhat es la fila B x_t del ejemplo 1 de la clase 9
  a <- ajustar_ingenuo(x_ej1)
  agregar("ingenuo", "Clase 9, ejemplo 1: fila B x_t", "NA 4 7 6 9 12",
          paste(a$yhat, collapse = " "), identical(a$yhat, c(NA, 4, 7, 6, 9, 12)))
  
  # 2. Ingenuo estacional (s = 2): su yhat es la fila B^2 x_t del mismo ejemplo
  a <- ajustar_ingenuo_estacional(x_ej1, s = 2)
  agregar("ingenuo estacional", "Clase 9, ejemplo 1: fila B^2 x_t (s = 2)", "NA NA 4 7 6 9",
          paste(a$yhat, collapse = " "), identical(a$yhat, c(NA, NA, 4, 7, 6, 9)))
  
  # 3. Media simple: media recursiva calculada a mano, 4, (4+7)/2, (4+7+6)/3, ...
  a <- ajustar_media(x_ej1)
  a_mano <- c(NA, 4, 11 / 2, 17 / 3, 26 / 4, 38 / 5)
  agregar("media simple", "A mano sobre el ejemplo 1: 4; 5,5; 5,667; 6,5; 7,6", "media de y1..y(t-1)",
          paste(round(a$yhat, 3), collapse = " "),
          ok(a$yhat[-1], a_mano[-1], 1e-12) && a$pronosticar(1) == mean(x_ej1))
  
  # 4. Media movil MM(3): precio de la accion (clase 3, diap. 45)
  a <- ajustar_mm(precio, 3); e <- err(a, precio)
  agregar("media movil", "Clase 3, diap. 45: MM(3), MSE, MAD, MAPE y pronostico de enero",
          "N = 5; 1,1621; 0,9780; 4,607 %; 21,5230",
          sprintf("N = %d; %.4f; %.4f; %.3f %%; %.4f", length(e), mse(e), mad(e),
                  100 * mean(abs(e / precio[4:8])), a$pronosticar(1)),
          length(e) == 5 && ok(c(mse(e), mad(e), 100 * mean(abs(e / precio[4:8])), a$pronosticar(1)),
                               c(1.1621, 0.9780, 4.607, 21.5230), 5e-4))
  
  # 5. SES(0,15): misma serie (clase 3, diap. 37), con L1 = Y1, de modo que Yhat_2 = Y1
  a <- ajustar_ses(precio, 0.15); e <- err(a, precio); e4 <- err(a, precio, 4)
  agregar("SES", "Clase 3, diap. 37: SES(0,15), N = 7 y MSE en origen comun t = 4..8",
          "N = 7; 2,8534; 1,5895; 7,609 %; comun 3,3706",
          sprintf("N = %d; %.4f; %.4f; %.3f %%; comun %.4f", length(e), mse(e), mad(e),
                  100 * mean(abs(e / precio[2:8])), mse(e4)),
          length(e) == 7 && ok(c(mse(e), mad(e), 100 * mean(abs(e / precio[2:8])), mse(e4)),
                               c(2.8534, 1.5895, 7.609, 3.3706), 5e-4))
  
  # 6. Doble media movil DMM(3): creditos al consumo (clase 4, diap. 66)
  a <- ajustar_dmm(creditos, 3); e <- err(a, creditos)
  agregar("doble media movil", "Clase 4, diap. 66: DMM(3), N, MSE, MAD y 3 pronosticos",
          "N = 6; 841,24; 28,20; 356,56 370,00 383,44",
          sprintf("N = %d; %.2f; %.2f; %s", length(e), mse(e), mad(e),
                  paste(sprintf("%.2f", a$pronosticar(3)), collapse = " ")),
          length(e) == 6 && ok(c(mse(e), mad(e)), c(841.24, 28.20), 5e-3) &&
            ok(a$pronosticar(3), c(356.56, 370.00, 383.44), 5e-3))
  
  # 7. Holt(0,3; 0,4): creditos al consumo, tabla completa (clase 4, diap. 66)
  a <- ajustar_holt(creditos, 0.3, 0.4); e <- err(a, creditos)
  yhat_notas <- c(NA, 133.000, 142.240, 154.439, 166.766, 185.563, 215.273,
                  256.017, 302.656, 336.845, 366.315)
  agregar("Holt", "Clase 4, diap. 66: Holt(0,3; 0,4), columna Yhat, MSE, MAD y 3 pronosticos",
          "N = 10; 1133,08; 28,62; 387,15 414,97 442,80",
          sprintf("N = %d; %.2f; %.2f; %s", length(e), mse(e), mad(e),
                  paste(sprintf("%.2f", a$pronosticar(3)), collapse = " ")),
          length(e) == 10 && ok(a$yhat[-1], yhat_notas[-1], 5e-4) &&
            ok(c(mse(e), mad(e)), c(1133.08, 28.62), 5e-3) &&
            ok(a$pronosticar(3), c(387.15, 414.97, 442.80), 5e-3))
  
  # 8. Tendencia con ficticias de suma cero: ventas trimestrales (clase 4, diap. 53 y 55)
  a <- ajustar_tendencia_estacional(ventas, s = 4, inicio = 1)
  cf <- c(a$parametros$coeficientes[c("intercepto", "t")], a$parametros$efectos_estacionales)
  agregar("tendencia con ficticias", "Clase 4, diap. 53 y 55: coeficientes y Yhat_17..20",
          "6,6644 0,0262 | -5,2481 -2,2494 1,7494 5,7481 | 1,863 4,887 8,912 12,938",
          paste(sprintf("%.4f", cf), collapse = " "),
          ok(cf, c(6.6644, 0.0262, -5.2481, -2.2494, 1.7494, 5.7481), 5e-4) &&
            ok(a$pronosticar(4), c(1.863, 4.887, 8.912, 12.938), 5e-4))
  
  # 9. Tendencias lineal, cuadratica y exponencial: contra lm() y contra Duan a mano,
  #    sobre los creditos; el yhat de t se compara con lm() ajustado en y[1..t-1].
  tt <- seq_along(creditos)
  for (tipo in c("lineal", "cuadratica", "exponencial")) {
    a <- ajustar_tendencia(creditos, tipo)
    form <- switch(tipo, lineal = y ~ t, cuadratica = y ~ t + I(t^2), exponencial = log(y) ~ t)
    d <- data.frame(y = creditos, t = tt)
    fit <- stats::lm(form, data = d)
    nuevo <- data.frame(t = 12:14)
    duan <- if (tipo == "exponencial") mean(exp(stats::residuals(fit))) else 1
    pr_lm <- as.numeric(stats::predict(fit, nuevo))
    pr_lm <- if (tipo == "exponencial") exp(pr_lm) * duan else pr_lm
    t0 <- 9                                                 # un yhat dentro de muestra
    fit9 <- stats::lm(form, data = d[1:(t0 - 1), ])
    y9 <- as.numeric(stats::predict(fit9, data.frame(t = t0)))
    if (tipo == "exponencial") y9 <- exp(y9) * mean(exp(stats::residuals(fit9)))
    agregar(paste("tendencia", tipo), "lm() de R: coeficientes, pronosticos 12..14 y yhat_9 con y1..y8",
            paste(sprintf("%.4f", stats::coef(fit)), collapse = " "),
            paste(sprintf("%.4f", a$parametros$coeficientes), collapse = " "),
            ok(a$parametros$coeficientes, stats::coef(fit), 1e-8) &&
              ok(a$pronosticar(3), pr_lm, 1e-8) && ok(a$yhat[t0], y9, 1e-8))
  }
  
  # 9b. regresion_mco() (usada en 6b): coeficientes y error estandar ordinario contra
  #     summary(lm()) sobre las ventas con ficticias de suma cero; el HAC con L = 0 debe
  #     coincidir con el error de White, (X'X)^-1 X' diag(e^2) X (X'X)^-1, calculado aparte.
  X <- disenar_tendencia_estacional(1:16, 4, 1)
  rm_ <- regresion_mco(ventas, X)
  q <- factor(rep(1:4, 4))
  fit <- stats::lm(ventas ~ t + q, data = data.frame(ventas = ventas, t = 1:16, q = q),
                   contrasts = list(q = "contr.sum"))
  sm <- summary(fit)$coefficients
  white <- sqrt(diag(solve(crossprod(X)) %*% crossprod(X * rm_$residuos) %*% solve(crossprod(X))))
  agregar("regresion_mco() (6b)", "summary(lm()) con contr.sum: coeficientes y e.e. ordinario; HAC con L = 0 = White",
          paste(sprintf("%.4f", sm[, 2]), collapse = " "),
          paste(sprintf("%.4f", rm_$ee_ordinario), collapse = " "),
          ok(rm_$coeficientes, sm[, 1], 1e-10) && ok(rm_$ee_ordinario, sm[, 2], 1e-10) &&
            ok(regresion_mco(ventas, X, L = 0)$ee_robusto, white, 1e-10))
  
  # 10. Holt-Winters aditivo: ventas, (alpha, beta, gamma) = (0,3; 0,05; 0,3), paso t = 5 (clase 5, diap. 54)
  a <- ajustar_hw(ventas, 4, 0.3, 0.05, 0.3)
  hw_obt <- c(a$estados$L[4], a$estados$b[4], a$estados$S[1], a$yhat[5],
              a$estados$L[5], a$estados$b[5], a$estados$S[5])
  agregar("Holt-Winters aditivo", "Clase 5, diap. 54: L4, b4, S1, Yhat5, L5, b5, S5",
          "6,8650 0,0326 -5,2471 1,6505 6,9425 0,0348 -5,2157",
          paste(sprintf("%.4f", hw_obt), collapse = " "),
          ok(hw_obt, c(6.8650, 0.0326, -5.2471, 1.6505, 6.9425, 0.0348, -5.2157), 5e-5) &&
            abs(sum(a$parametros$inicial$S0)) < 1e-12)
  
  # 11. Holt-Winters multiplicativo: no hay cifra en las notas; calculo a mano sobre las ventas.
  #     Promedios anuales 6,80; 6,85; 6,70; 7,20 -> b0 = 0,4/12, L0 = 6,80 - 1,5 b0 = 6,75.
  #     S0_j = promedio de Y / (L0 + b0 (t - 1)), dividido por su media para que promedie uno.
  #     Paso t = 1: Yhat1 = (L0 + b0) S0_1, L1 = a Y1/S0_1 + (1 - a)(L0 + b0), etc.
  a <- ajustar_hw_mult(ventas, 4, 0.3, 0.05, 0.3)
  b0 <- (7.2 - 6.8) / 12; L0 <- 6.8 - 1.5 * b0
  recta <- L0 + b0 * (0:15)
  S0 <- c(mean(ventas[c(1, 5, 9, 13)] / recta[c(1, 5, 9, 13)]),
          mean(ventas[c(2, 6, 10, 14)] / recta[c(2, 6, 10, 14)]),
          mean(ventas[c(3, 7, 11, 15)] / recta[c(3, 7, 11, 15)]),
          mean(ventas[c(4, 8, 12, 16)] / recta[c(4, 8, 12, 16)]))
  S0 <- S0 / mean(S0)
  P1 <- L0 + b0
  yhat1 <- P1 * S0[1]
  L1 <- 0.3 * 1.5 / S0[1] + 0.7 * P1
  b1 <- 0.05 * (L1 - L0) + 0.95 * b0
  S1 <- 0.3 * 1.5 / L1 + 0.7 * S0[1]
  yhat2 <- (L1 + b1) * S0[2]
  mano <- c(L0, b0, S0[1], yhat1, L1, b1, S1, yhat2)
  obt <- c(a$parametros$inicial$L0, a$parametros$inicial$b0, a$parametros$inicial$S0[1],
           a$yhat[1], a$estados$L[1], a$estados$b[1], a$estados$S[1], a$yhat[2])
  agregar("Holt-Winters multiplicativo", "A mano sobre las ventas: L0, b0, S0_1, Yhat1, L1, b1, S1, Yhat2",
          paste(sprintf("%.4f", mano), collapse = " "),
          paste(sprintf("%.4f", obt), collapse = " "),
          ok(obt, mano, 1e-12) && abs(mean(a$parametros$inicial$S0) - 1) < 1e-12)
  
  # 12. Con alpha = beta = gamma = 0, Holt-Winters aditivo pronostica la recta inicial mas
  #     el perfil inicial (clase 5, proposicion 8): Yhat_t = L0 + b0 t + S0_j.
  a <- ajustar_hw(ventas, 4, 0, 0, 0)
  ini <- a$parametros$inicial
  agregar("Holt-Winters aditivo", "Clase 5, prop. 8: con alpha = beta = gamma = 0, Yhat_t = L0 + b0 t + S0_j",
          "diferencia 0", sprintf("%.1e", max(abs(a$yhat - (ini$L0 + ini$b0 * 1:16 + rep(ini$S0, 4))))),
          ok(a$yhat, ini$L0 + ini$b0 * 1:16 + rep(ini$S0, 4), 1e-12))
  
  dplyr::bind_rows(filas)
})
print(as.data.frame(verif6), right = FALSE)
stopifnot("Algun metodo no reproduce la cifra de verificacion" = all(verif6$verificado))
message("Verificaciones superadas: los ", nrow(verif6), " contrastes de la seccion 6.")


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
# SECCION 6. PRONOSTICO Y EVALUACION
# =============================================================================

# ---- 6(a) Catalogo de metodos y asignacion a cada serie ---------------------
# Los metodos se aplican a la serie ORIGINAL (en sus unidades), porque cada uno trae su
# propio supuesto sobre nivel, tendencia y estacionalidad; la transformacion de la
# seccion 4 se usa como diagnostico para decidir cuales supuestos son sostenibles.
# La verificacion de cada ajustar_*() esta en el bloque V6.

# Catalogo: funcion, rejilla de rejillas() (NULL si no tiene constantes) y argumentos fijos.
catalogo_6 <- list(
  media           = list(nombre = "Media simple",               f = ajustar_media,     rejilla = NULL),
  ingenuo         = list(nombre = "Ingenuo",                    f = ajustar_ingenuo,   rejilla = NULL),
  ingenuo_est     = list(nombre = "Ingenuo estacional",         f = ajustar_ingenuo_estacional, rejilla = NULL),
  mm              = list(nombre = "Media movil",                f = ajustar_mm,        rejilla = "mm"),
  ses             = list(nombre = "Suavizamiento exponencial simple", f = ajustar_ses, rejilla = "ses"),
  dmm             = list(nombre = "Doble media movil",          f = ajustar_dmm,       rejilla = "dmm"),
  tend_lineal     = list(nombre = "Tendencia lineal",           f = ajustar_tendencia, rejilla = NULL,
                         fijos = list(tipo = "lineal")),
  tend_cuadratica = list(nombre = "Tendencia cuadratica",       f = ajustar_tendencia, rejilla = NULL,
                         fijos = list(tipo = "cuadratica")),
  tend_exponencial = list(nombre = "Tendencia exponencial",     f = ajustar_tendencia, rejilla = NULL,
                          fijos = list(tipo = "exponencial")),
  tend_ficticias  = list(nombre = "Tendencia con ficticias de suma cero", f = ajustar_tendencia_estacional,
                         rejilla = NULL),
  holt            = list(nombre = "Holt",                       f = ajustar_holt,      rejilla = "holt"),
  hw              = list(nombre = "Holt-Winters aditivo",       f = ajustar_hw,        rejilla = "hw"),
  hw_mult         = list(nombre = "Holt-Winters multiplicativo", f = ajustar_hw_mult,  rejilla = "hw_mult")
)

# Metodos que usan el periodo s; a cada serie se le pasa su s_candidato (7, 4, 12).
# inicio = 1 en las tres series: la serie 1 empieza un lunes (estacion 1 con %u),
# la 2 en el primer trimestre y la 3 en enero.
usa_s <- c("ingenuo_est", "tend_ficticias", "hw", "hw_mult")

# DECISION DEL EQUIPO: asignacion de metodos (edite aqui si cambia el diagnostico).
#   rol = "sostenido"  el diagnostico de las secciones 3 a 5 sostiene su supuesto
#         "referente"  ingenuo con que se compara todo (6c)
#         "contraste"  el diagnostico contradice su supuesto; se aplica para mostrarlo
#         "no aplica"  no se aplica a esa serie
asignacion_6a <- tibble::tribble(
  ~serie,   ~metodo,            ~rol,         ~evidencia,
  # ---- Serie 1: tipo I, nivel estable, r1 = 0,48, sin estacionalidad relevante -----------
  "serie1", "media",            "sostenido",  "Tipo I (3c): nivel estable sin tendencia; no se diferencia (4c: Var(nabla Y)/Var(Y) = 1,03 = 2(1 - r1)). La media simple estima ese nivel fijo.",
  "serie1", "ingenuo",          "referente",  "Sin estacionalidad (3b: r7 = 0,09 dentro de +/-0,146), luego el referente es el ingenuo (s = 1), no el estacional.",
  "serie1", "ingenuo_est",      "contraste",  "Supone que se repite el dia de la semana; 3b no le atribuye estacionalidad semanal (efecto pequeno, r7 dentro de la banda).",
  "serie1", "mm",               "sostenido",  "Estacionaria en media con dependencia de corto plazo (r1 = 0,48, ACF que cae rapido): un nivel local promediado sobre k dias.",
  "serie1", "ses",              "sostenido",  "Mismo supuesto de nivel local; con r1 = 0,48 el pronostico optimo vuelve a la media sin olvidar del todo el ultimo dato.",
  "serie1", "dmm",              "contraste",  "Supone tendencia lineal local; la pendiente ajustada es 0,015 s/dia (2,7 s en 180 dias frente a s = 18 s).",
  "serie1", "tend_lineal",      "contraste",  "Supone tendencia global; la pendiente de 0,015 s/dia es despreciable y diferenciar sube la varianza (4c).",
  "serie1", "tend_cuadratica",  "no aplica",  "Sin tendencia, menos aun curva; la recta ya es el contraste de tendencia global.",
  "serie1", "tend_exponencial", "no aplica",  "Sin crecimiento porcentual: el nivel no crece.",
  "serie1", "tend_ficticias",   "no aplica",  "Supone tendencia y estacionalidad; 3b no atribuye estacionalidad y no hay tendencia.",
  "serie1", "holt",             "contraste",  "Supone pendiente local; el diagnostico dice nivel estable (tipo I). Si beta optima es grande, rastrea ruido.",
  "serie1", "hw",               "contraste",  "Supone estacionalidad s = 7; la F rechaza pero el efecto es pequeno y r7 queda en la banda: muestra si esa diferencia sirve para pronosticar.",
  "serie1", "hw_mult",          "no aplica",  "Sin estacionalidad ni amplitud que crezca con el nivel.",
  # ---- Serie 2: tipo IV, tendencia lineal (15,5 miles/trim.) + estacionalidad aditiva s = 4 --
  "serie2", "media",            "contraste",  "Supone nivel constante; la media anual crece de 1996 a 2025 y la ACF decae lentamente (r24 = 0,41).",
  "serie2", "ingenuo",          "contraste",  "Repite el ultimo trimestre; ignora la estacionalidad s = 4 (F rechaza con p < 0,001).",
  "serie2", "ingenuo_est",      "referente",  "Estacional con s = 4 (3b, r4 = 0,90 > r2 = 0,87): el referente natural es el ingenuo estacional, que repite el mismo trimestre del ano anterior.",
  "serie2", "mm",               "no aplica",  "Supone nivel constante; ya cubierto por el contraste de la media simple.",
  "serie2", "ses",              "no aplica",  "Supone nivel sin tendencia ni estacionalidad.",
  "serie2", "dmm",              "no aplica",  "Supone tendencia sin estacionalidad; el contraste de tendencia sin estacion lo hacen la recta y Holt.",
  "serie2", "tend_lineal",      "contraste",  "Capta la tendencia lineal pero no la estacionalidad: sus residuos deben conservar el patron de s = 4.",
  "serie2", "tend_cuadratica",  "no aplica",  "La tendencia es lineal: (1-B^4) deja media 60,4 = 4 x 15,5 y (1-B)(1-B^4) deja media 0,56, a 0,09 e.e. de cero.",
  "serie2", "tend_exponencial", "no aplica",  "Sin crecimiento porcentual: la amplitud no crece con el nivel y se eligio la version aditiva.",
  "serie2", "tend_ficticias",   "sostenido",  "Tipo IV: tendencia lineal mas estacionalidad aditiva de periodo 4 (3b y 3c); (1-B)(1-B^4) es justamente el operador que anula ese modelo.",
  "serie2", "holt",             "contraste",  "Rastrea nivel y pendiente pero no tiene estado estacional; el diagnostico exige s = 4.",
  "serie2", "hw",               "sostenido",  "Tendencia y estacionalidad aditiva s = 4 (3b): la amplitud casi no cambia mientras el nivel sube. Version local del modelo con ficticias.",
  "serie2", "hw_mult",          "contraste",  "Supone amplitud proporcional al nivel; 3b eligio aditiva porque la amplitud no crece con el nivel.",
  # ---- Serie 3: tipo III, tendencia curva con crecimiento de 1,24 % mensual, sin estacionalidad --
  "serie3", "media",            "no aplica",  "Supone nivel constante con una serie que pasa de 83 a 485 miles; el contraste de nivel lo hace el SES.",
  "serie3", "ingenuo",          "referente",  "Sin estacionalidad (F no rechaza; r12 = 0,72 entre r11 = 0,74 y r13 = 0,70): el referente es el ingenuo, pronostico optimo del paseo aleatorio que sugiere nabla log Y.",
  "serie3", "ingenuo_est",      "contraste",  "Supone que se repite el mes del ano anterior; 3b no encuentra estacionalidad anual y la tendencia lo deja un ano atrasado.",
  "serie3", "mm",               "no aplica",  "Supone nivel constante; el contraste de nivel lo hace el SES.",
  "serie3", "ses",              "contraste",  "Supone nivel sin tendencia; la ACF de la serie decae lentamente (r1 = 0,96) y el log crece 1,24 % por mes.",
  "serie3", "dmm",              "sostenido",  "Tendencia local (3c) que puede cambiar de pendiente: la serie se acelera.",
  "serie3", "tend_lineal",      "contraste",  "Supone recta; la tendencia es curva y crece a tasa constante, no en cantidad constante.",
  "serie3", "tend_cuadratica",  "sostenido",  "Tendencia curva (3a: el crecimiento se acelera).",
  "serie3", "tend_exponencial", "sostenido",  "Crecimiento porcentual constante: media de nabla log Y = 0,0124, casi igual a la pendiente del log (0,0121); r1 = -0,42 de nabla log Y es compatible con tendencia determinista en log.",
  "serie3", "tend_ficticias",   "no aplica",  "Sin estacionalidad (3b).",
  "serie3", "holt",             "sostenido",  "Tendencia local; nivel y pendiente se rastrean y la pendiente puede crecer con el nivel.",
  "serie3", "hw",               "no aplica",  "Sin estacionalidad; el contraste estacional lo hace la version multiplicativa, que es la que corresponde a una variabilidad que crece con el nivel.",
  "serie3", "hw_mult",          "contraste",  "Supone estacionalidad s = 12 multiplicativa; la F no rechaza y el rango de los indices es pequeno: la amplitud crece por el crecimiento porcentual, no por una estacion."
)

# Comprobaciones de la asignacion: cada metodo del catalogo se aplica a alguna serie,
# cada serie tiene exactamente un referente, al menos un sostenido y al menos un contraste.
stopifnot(
  setequal(unique(asignacion_6a$metodo), names(catalogo_6)),
  all(table(asignacion_6a$serie) == length(catalogo_6)),
  all(names(catalogo_6) %in% asignacion_6a$metodo[asignacion_6a$rol != "no aplica"]),
  all(tapply(asignacion_6a$rol == "referente", asignacion_6a$serie, sum) == 1),
  all(tapply(asignacion_6a$rol == "sostenido", asignacion_6a$serie, sum) >= 1),
  all(tapply(asignacion_6a$rol == "contraste", asignacion_6a$serie, sum) >= 1),
  all(asignacion_6a$metodo[asignacion_6a$rol == "referente"] %in% c("ingenuo", "ingenuo_est"))
)

# Argumentos fijos de cada metodo en cada serie (s e inicio cuando corresponde)
fijos_6 <- function(metodo, id) {
  f <- catalogo_6[[metodo]]$fijos
  if (is.null(f)) f <- list()
  if (metodo %in% usa_s) f$s <- series_cfg[[id]]$s_candidato
  if (metodo == "tend_ficticias") f$inicio <- 1
  f
}

# Tabla larga de lo que se aplica: es la entrada de 6(b), 6(c) y 6(d)
aplicados_6 <- asignacion_6a |>
  dplyr::filter(rol != "no aplica") |>
  dplyr::mutate(nombre = vapply(metodo, function(m) catalogo_6[[m]]$nombre, character(1)),
                rejilla = vapply(metodo, function(m) {
                  r <- catalogo_6[[m]]$rejilla; if (is.null(r)) NA_character_ else r
                }, character(1)))

# Tabla de 6(a): metodo por serie, con el rol en cada una
tabla_6a <- asignacion_6a |>
  dplyr::mutate(nombre = vapply(metodo, function(m) catalogo_6[[m]]$nombre, character(1))) |>
  dplyr::select(nombre, metodo, serie, rol) |>
  tidyr::pivot_wider(names_from = serie, values_from = rol)
cat("\nTABLA 6(a): asignacion de metodos (rol en cada serie)\n")
print(as.data.frame(tabla_6a), right = FALSE)
cat("\nMetodos aplicados por serie:\n")
print(table(aplicados_6$serie, aplicados_6$rol))

# ---- 6(b) Regresiones: coeficientes, errores estandar y supuestos -----------
# Cada tendencia ajustada por MCO (con o sin ficticias) se ajusta sobre la serie completa
# (1, ..., T): no tiene constantes que elegir, y la inferencia sobre la tendencia usa todos
# los datos. Error estandar ordinario y robusto HAC con nucleo de Bartlett y
# L = floor(4 (T/100)^(2/9)) rezagos (regresion_mco()). Las pruebas t usan el robusto,
# referido a una t de Student con T - p g.l. Los residuos se revisan con el correlograma,
# Ljung-Box referido a chi^2_m (los coeficientes de regresores deterministas no restan
# grados de libertad), Jarque-Bera y Durbin-Watson, leido con d ~ 2(1 - r1).
# La exponencial se revisa en escala logaritmica, donde esta su error: ln Y = a + theta t + e.

etiqueta_coef <- function(nombre) {
  dplyr::case_when(
    nombre == "intercepto" ~ "beta0 (intercepto)",
    nombre == "t" ~ "beta1 (t)",
    nombre == "t2" ~ "beta2 (t^2)",
    TRUE ~ paste0("gamma", sub("^e", "", nombre), " (estacion ", sub("^e", "", nombre), ")"))
}

lectura_coef <- function(nombre, tipo, rechaza, b, unidad) {
  if (nombre == "intercepto") {
    return(if (rechaza) "El nivel en t = 0 es distinto de cero." else
      "No se distingue de cero el nivel en t = 0 (sin interes practico).")
  }
  if (nombre == "t") {
    que <- if (tipo == "exponencial") paste0("una tasa de crecimiento de ", formatC(100 * b, digits = 2, format = "f", decimal.mark = ","), " % por periodo")
    else paste0("un cambio de ", formatC(b, digits = 3, format = "f", decimal.mark = ","), " ", unidad, " por periodo")
    return(if (rechaza) paste0("Hay tendencia significativa: ", que, ".") else
      paste0("La pendiente no se distingue de cero con el error robusto: no hay tendencia (estimada: ", que, ")."))
  }
  if (nombre == "t2") {
    return(if (rechaza) "La curvatura es significativa: la tendencia no es una recta." else
      "La curvatura no se distingue de cero: no hay evidencia contra la recta.")
  }
  j <- sub("^e", "", nombre)
  if (rechaza) paste0("La estacion ", j, " se aparta del nivel de su ano en ", formatC(b, digits = 1, format = "f", decimal.mark = ","), " ", unidad, ".")
  else paste0("La estacion ", j, " no se distingue del nivel de su ano.")
}

hacer_6b <- function(serie, metodo, rol, nombre, ...) {
  cfg <- series_cfg[[serie]]
  datos <- datos_series[[serie]]
  y <- datos$valor
  n <- length(y)
  aj <- do.call(catalogo_6[[metodo]]$f, c(list(y = y), fijos_6(metodo, serie)))
  reg <- aj$regresion
  tipo <- aj$parametros$tipo %||% "ficticias"
  unidad <- if (tipo == "exponencial") "(log)" else cfg$unidad
  
  # Coeficientes; con ficticias se agrega el efecto recuperado gamma_s = -(gamma_1 + ... + gamma_{s-1})
  b <- reg$coeficientes; ee_o <- reg$ee_ordinario; ee_r <- reg$ee_robusto
  if (metodo == "tend_ficticias") {
    s <- fijos_6(metodo, serie)$s
    w <- as.numeric(grepl("^e", names(b)))
    b <- c(b, setNames(-sum(b * w), paste0("e", s)))
    ee_o <- c(ee_o, sqrt(drop(t(w) %*% reg$V_ordinario %*% w)))
    ee_r <- c(ee_r, sqrt(drop(t(w) %*% reg$V_robusto %*% w)))
  }
  gl <- reg$gl
  t_rob <- b / ee_r
  critico <- stats::qt(0.975, gl)
  p_rob <- 2 * stats::pt(abs(t_rob), gl, lower.tail = FALSE)
  coeficientes <- tibble::tibble(
    coeficiente = names(b), etiqueta = etiqueta_coef(names(b)), estimacion = unname(b),
    ee_ordinario = unname(ee_o), ee_robusto = unname(ee_r), razon_ee = unname(ee_r / ee_o),
    t_robusto = unname(t_rob), critico = critico, p_valor = unname(p_rob))
  
  # Pruebas t con el error robusto, con sus seis elementos
  pruebas_t <- lapply(seq_along(b), function(i) {
    nm <- names(b)[i]
    prueba_seis(
      nombre = paste0("t robusta de ", etiqueta_coef(nm)),
      h0 = paste0(etiqueta_coef(nm), " = 0"), h1 = paste0(etiqueta_coef(nm), " != 0"),
      estadistico = "t = b / ee_HAC(b)",
      distribucion = paste0("t de Student con ", gl, " g.l. (aproximada)"),
      observado = t_rob[i], critico = c(-critico, critico),
      region = paste0("|t| > ", formatC(critico, digits = 3, format = "f")),
      p = p_rob[i],
      lectura_rechazo = lectura_coef(nm, tipo, TRUE, b[i], unidad),
      lectura_no_rechazo = lectura_coef(nm, tipo, FALSE, b[i], unidad))
  })
  
  # Supuestos sobre los residuos
  e <- reg$residuos
  m_acf <- min(floor(n / 4), 24)
  r_e <- acf_muestral(e, m_acf)
  que <- paste0("los residuos de la ", tolower(nombre))
  lb <- ljung_box(e, m_lb[[serie]], g = 0, que = que)
  val <- validar_errores(e, m_lb[[serie]])
  jb <- val$jarque_bera
  dw <- val$durbin_watson
  fuera <- which(abs(r_e) > 1.96 / sqrt(n))
  
  # Figura: residuos en el tiempo y su correlograma
  df_e <- data.frame(fecha = datos$fecha, e = e)
  g_e <- ggplot2::ggplot(df_e, ggplot2::aes(x = fecha, y = e)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_line(colour = "grey25", linewidth = 0.4) +
    ggplot2::labs(title = paste0(cfg$corto, ": residuos de la ", tolower(nombre)),
                  x = "Fecha", y = paste0("Residuo ", if (tipo == "exponencial") "(escala log)" else paste0("(", cfg$unidad, ")"))) +
    ggplot2::theme_minimal(base_size = 11)
  g_r <- correlograma(r_e, n, "ACF de los residuos", "ACF")
  id_fig <- paste0(serie, "-reg-", sub("tend_", "", metodo))
  guardar_fig(g_e / g_r, id_fig, ancho = 9, alto = 7)
  
  list(serie = serie, metodo = metodo, rol = rol, nombre = nombre, tipo = tipo, n = n,
       p = length(reg$coeficientes), gl = gl, L = reg$L, sigma2 = reg$sigma2,
       coeficientes = coeficientes, pruebas_t = pruebas_t,
       r_residuos = r_e, m_acf = m_acf, fuera_banda = fuera,
       ljung_box = lb, jarque_bera = jb, durbin_watson = dw, figura = id_fig)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

regresiones_6b <- aplicados_6 |>
  dplyr::filter(metodo %in% c("tend_lineal", "tend_cuadratica", "tend_exponencial", "tend_ficticias"))
res6b <- purrr::pmap(regresiones_6b, hacer_6b)
names(res6b) <- paste(regresiones_6b$serie, regresiones_6b$metodo, sep = "_")

resumen6b <- dplyr::bind_rows(lapply(res6b, function(z) {
  tibble::tibble(serie = z$serie, regresion = z$nombre, rol = z$rol, n = z$n, p = z$p, L = z$L,
                 LB_m = as.numeric(sub(".*m = (\\d+).*", "\\1", z$ljung_box$nombre)),
                 LB_Q = z$ljung_box$observado, LB_p = z$ljung_box$p_valor,
                 JB = z$jarque_bera$observado, JB_p = z$jarque_bera$p_valor,
                 DW = z$durbin_watson[["d"]], r1 = z$durbin_watson[["r1"]],
                 DW_aprox = z$durbin_watson[["aprox_2_1_menos_r1"]],
                 rezagos_fuera = compactar(z$fuera_banda))
}))
cat("\nRESUMEN 6(b): supuestos sobre los residuos de cada regresion\n")
print(as.data.frame(resumen6b), digits = 3)
for (z in res6b) {
  cat("\n", strrep("=", 70), "\n", z$serie, ": ", z$nombre, "\n", sep = "")
  print(as.data.frame(z$coeficientes), digits = 4)
}


# ---- 6(c) Origen movil -------------------------------------------------------
# Para cada serie y cada metodo aplicado en 6(a), incluido el contraste y el referente:
#   1. las constantes se eligen con optimizar() SOLO en y[1..n0], sobre la rejilla del
#      catalogo (la curva o el mapa de MSE se guarda en figs/);
#   2. se evalua por origen movil con esas constantes fijas, con los mismos origenes para
#      todos los metodos de la serie (h = 1: o = n0..T-1; h = H: o = n0..T-H).
# El referente de cada serie es el ingenuo que asignacion_6a marca con rol "referente"
# (la justificacion esta en la asignacion y en el informe). Es una fila mas de la tabla,
# calculada sobre los mismos origenes. Un metodo que no lo supera no se descarta.
# El GANADOR de cada serie es el de menor RMSE a h = 1; tambien se guarda el mejor a h = H.

hacer_6c <- function(id, cfg, datos, aplic, part) {
  y <- datos$valor
  n0 <- part$n0[1]
  H <- part$H[1]
  stopifnot("n0 distinto de floor(0.8 T)" = n0 == floor(0.8 * length(y)),
            "H distinto del de la configuracion" = H == cfg$H)
  rej_todas <- rejillas()
  evals <- list()
  filas <- list()
  for (i in seq_len(nrow(aplic))) {
    m <- aplic$metodo[i]
    rej <- if (is.na(aplic$rejilla[i])) NULL else rej_todas[[aplic$rejilla[i]]]
    message(id, ": 6(c) ", aplic$nombre[i], " ...")
    ev <- evaluar_metodo(y, catalogo_6[[m]]$f, n0, H, rej, fijos_6(m, id))
    evals[[m]] <- ev
    if (!is.null(ev$curva)) {
      npar <- ncol(ev$curva) - 1L
      guardar_fig(
        grafico_mse(ev$curva, paste0(cfg$corto, ": MSE de ", aplic$nombre[i])),
        paste0(id, "-6c-mse-", m),
        ancho = if (npar == 1L) 6.5 else if (npar == 2L) 7 else 10,
        alto  = if (npar == 1L) 3.8 else if (npar == 2L) 5 else 7.5)
    }
    filas[[i]] <- tibble::tibble(
      metodo = m, nombre = aplic$nombre[i], rol = aplic$rol[i],
      parametros = texto_parametros(ev$parametros), mse_est = ev$mse_est,
      n1 = ev$m1$n, RMSE1 = ev$m1$RMSE, MAD1 = ev$m1$MAD, MAPE1 = ev$m1$MAPE,
      nH = ev$mH$n, RMSEH = ev$mH$RMSE, MADH = ev$mH$MAD, MAPEH = ev$mH$MAPE)
  }
  tabla <- dplyr::bind_rows(filas)
  
  # Una tabla por horizonte, ordenada de menor a mayor RMSE, con el RMSE relativo al referente
  armar <- function(h1) {
    if (h1) {
      t <- dplyr::transmute(tabla, metodo, nombre, rol, parametros, n = n1,
                            RMSE = RMSE1, MAD = MAD1, MAPE = MAPE1)
    } else {
      t <- dplyr::transmute(tabla, metodo, nombre, rol, parametros, n = nH,
                            RMSE = RMSEH, MAD = MADH, MAPE = MAPEH)
    }
    t <- t[order(t$RMSE), ]
    t$RMSE_rel <- t$RMSE / t$RMSE[t$rol == "referente"]
    t
  }
  tabla_h1 <- armar(TRUE)
  tabla_hH <- armar(FALSE)
  
  # Origenes: comprobacion de que todos los metodos usan los mismos
  stopifnot("distinto numero de errores entre metodos" =
              length(unique(tabla$n1)) == 1L && length(unique(tabla$nH)) == 1L,
            "numero de origenes h = 1 distinto de T - n0" = tabla$n1[1] == length(y) - n0,
            "numero de origenes h = H distinto de T - H - n0 + 1" = tabla$nH[1] == length(y) - H - n0 + 1)
  
  list(tabla = tabla, tabla_h1 = tabla_h1, tabla_hH = tabla_hH, evals = evals,
       ganador = tabla_h1$metodo[1], ganador_H = tabla_hH$metodo[1],
       referente = tabla$metodo[tabla$rol == "referente"],
       n0 = n0, H = H, T = length(y))
}

res6c <- purrr::imap(series_cfg, function(cfg, id) {
  hacer_6c(id, cfg, datos_series[[id]], aplicados_6[aplicados_6$serie == id, ],
           particiones[particiones$serie == id, ])
})

resumen6c <- dplyr::bind_rows(lapply(names(res6c), function(id) {
  z <- res6c[[id]]
  g <- z$tabla[z$tabla$metodo == z$ganador, ]
  r <- z$tabla[z$tabla$metodo == z$referente, ]
  gH <- z$tabla[z$tabla$metodo == z$ganador_H, ]
  tibble::tibble(serie = id, transformacion = transformacion_txt[[id]],
                 ganador = g$nombre, parametros = g$parametros, referente = r$nombre,
                 RMSE1_ganador = g$RMSE1, RMSE1_referente = r$RMSE1,
                 ganador_H = gH$nombre, parametros_H = gH$parametros,
                 RMSEH_ganador = gH$RMSEH, RMSEH_referente = r$RMSEH)
}))
cat("\nRESUMEN 6(c): ganador a h = 1 y mejor a h = H, frente al referente\n")
print(as.data.frame(resumen6c), digits = 4)
for (id in names(res6c)) {
  z <- res6c[[id]]
  cat("\n", strrep("=", 70), "\n", id, " (n0 = ", z$n0, ", H = ", z$H, ")  h = 1:\n", sep = "")
  print(as.data.frame(z$tabla_h1), digits = 4)
  cat("\nh = ", z$H, ":\n", sep = "")
  print(as.data.frame(z$tabla_hH), digits = 4)
}

# ---- 6(d) Validacion del ganador sobre sus errores de origen movil a un paso -----
# Errores en el tiempo, correlograma (banda sobre el numero de errores), prueba t de media
# cero, Ljung-Box referido a chi^2_m y Jarque-Bera con sus seis elementos, y Durbin-Watson.
# m de Ljung-Box = min(m_lb de la serie, floor(n/4)): con pocos errores no se piden mas
# rezagos de los que sostienen.

res6d <- purrr::imap(series_cfg, function(cfg, id) {
  z <- res6c[[id]]
  ev <- z$evals[[z$ganador]]
  v <- validar_ganador(ev$om, m_lb[[id]])
  fechas_obj <- datos_series[[id]]$fecha[v$objetivo]
  guardar_fig(grafico_errores_6d(fechas_obj, v$e, v$r, v$n, cfg$corto, cfg$unidad),
              paste0(id, "-6d-errores"), ancho = 11, alto = 4.5)
  v$fechas <- fechas_obj
  v$metodo <- z$ganador
  v$nombre <- z$tabla$nombre[z$tabla$metodo == z$ganador]
  v$lectura_rezagos <- lectura_rezagos(v$fuera, cfg$s_candidato)
  v
})
for (id in names(res6d)) {
  v <- res6d[[id]]
  cat("\n", strrep("=", 70), "\n", id, ": 6(d) errores a un paso de ", v$nombre,
      " (n = ", v$n, ", m del correlograma = ", v$m, ", m de Ljung-Box = ", v$m_test, ")\n", sep = "")
  print(v$validacion$t); cat("\n")
  print(v$validacion$ljung_box); cat("\n")
  print(v$validacion$jarque_bera); cat("\n")
  print(round(v$validacion$durbin_watson, 3))
  cat(v$lectura_rezagos, "\n")
}

# ---- 6(e) Pronostico final ----------------------------------------------------------
# El ganador reajustado sobre toda la serie, con las constantes elegidas en y[1..n0];
# H pronosticos junto al ultimo tramo de la serie, con una banda HEURISTICA de +/- 1,96
# veces la desviacion de los errores de origen movil a cada horizonte.
ultimo_tramo <- c(serie1 = 45, serie2 = 40, serie3 = 48)

res6e <- purrr::imap(series_cfg, function(cfg, id) {
  z <- res6c[[id]]
  ev <- z$evals[[z$ganador]]
  d <- datos_series[[id]]
  nombre_g <- z$tabla$nombre[z$tabla$metodo == z$ganador]
  pf <- pronostico_final(d$valor, d$fecha, catalogo_6[[z$ganador]]$f, ev$parametros, z$H, ev$om)
  umbral <- if (id == "serie1") 330 else NULL        # meta contractual de la serie 1
  guardar_fig(grafico_pronostico(d$fecha, d$valor, pf, ultimo_tramo[[id]], cfg$corto, cfg$unidad,
                                 nombre_g, umbral),
              paste0(id, "-6e-pronostico"))
  out <- list(metodo = z$ganador, nombre = nombre_g, parametros = texto_parametros(ev$parametros),
              pronostico = pf, errores_h = errores_por_horizonte(ev$om))
  if (id == "serie1") {
    # Insumo para la seccion 8: dias de la semana siguiente (h = 3..9) con pronostico sobre la meta
    # y probabilidad heuristica de superarla, suponiendo errores normales con la desviacion de cada h.
    sub <- pf[pf$h >= 3, ]
    out$meta <- 330
    out$sobre_meta <- tibble::tibble(h = sub$h, fecha = sub$fecha, pronostico = sub$pronostico,
                                     sd_error = sub$sd_error,
                                     prob_supera = 1 - stats::pnorm((330 - sub$pronostico) / sub$sd_error))
    out$dias_esperados <- sum(out$sobre_meta$prob_supera)
  }
  out
})
for (id in names(res6e)) {
  cat("\n", strrep("=", 70), "\n", id, ": 6(e) pronostico final con ", res6e[[id]]$nombre,
      " (", res6e[[id]]$parametros, ")\n", sep = "")
  print(as.data.frame(res6e[[id]]$pronostico), digits = 4)
}
cat("\nSerie 1, insumo para 8: probabilidad heuristica de superar 330 s en h = 3..9\n")
print(as.data.frame(res6e$serie1$sobre_meta), digits = 3)
cat("Dias esperados sobre la meta (suma de probabilidades):", round(res6e$serie1$dias_esperados, 2), "\n")
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
       res5 = res5, resumen5 = resumen5,
       verif6 = verif6, asignacion_6a = asignacion_6a, tabla_6a = tabla_6a,
       aplicados_6 = aplicados_6, res6b = res6b, resumen6b = resumen6b,
       res6c = res6c, resumen6c = resumen6c, res6d = res6d, res6e = res6e,
       transformacion_txt = transformacion_txt, m_lb = m_lb),
  file.path("figs", "resultados.rds")
)

writeLines(capture.output(sessionInfo()), "sesion-info.txt")
tiempo_ejecucion <- difftime(Sys.time(), inicio_ejecucion, units = "secs")
message("correr-todo.R termino en ", round(as.numeric(tiempo_ejecucion), 1), " segundos.")
