
particion <- function(T, H, fechas = NULL) {
  if (!is.numeric(T) || T < 10) stop("particion(): `T` debe ser un entero >= 10.")
  if (!is.numeric(H) || H < 4) stop("particion(): el horizonte `H` debe ser >= 4.")
  n0 <- floor(0.8 * T)
  out <- tibble::tibble(
    T = T, n0 = n0, H = H,
    origenes_h1 = (T - 1) - n0 + 1,
    origenes_H = (T - H) - n0 + 1
  )
  if (!is.null(fechas)) {
    out$estimacion_desde <- fechas[1]
    out$estimacion_hasta <- fechas[n0]
    out$evaluacion_desde <- fechas[n0 + 1]
    out$evaluacion_hasta <- fechas[T]
    out$ultimo_origen_h1 <- fechas[T - 1]
    out$ultimo_origen_H <- fechas[T - H]
  }
  out
}

#' Rejillas del catalogo (seccion 6). Cambiar aqui si el curso indica otra.
rejillas <- function() {
  p05 <- seq(0.05, 0.95, by = 0.05)
  p10 <- seq(0.1, 0.9, by = 0.1)
  list(
    mm = data.frame(k = 2:12),
    ses = data.frame(alpha = seq(0.02, 0.98, by = 0.02)),
    dmm = data.frame(k = 3:12),
    holt = expand.grid(alpha = p05, beta = p05),
    hw = expand.grid(alpha = p10, beta = p10, gamma = p10),
    hw_mult = expand.grid(alpha = p10, beta = p10, gamma = p10)
  )
}


optimizar <- function(y, metodo, rejilla, n0, fijos = list()) {
  if (!is.function(metodo)) stop("optimizar(): `metodo` debe ser una funcion ajustar_*().")
  if (!is.data.frame(rejilla) || nrow(rejilla) < 1L) stop("optimizar(): `rejilla` debe ser un data.frame no vacio.")
  if (n0 > length(y)) stop("optimizar(): `n0` es mayor que la longitud de y.")
  y_est <- y[seq_len(n0)]
  yhats <- lapply(seq_len(nrow(rejilla)), function(i) {
    args <- c(list(y = y_est), fijos, as.list(rejilla[i, , drop = FALSE]))
    do.call(metodo, args)$yhat
  })
  comunes <- Reduce(`&`, lapply(yhats, function(v) !is.na(v)))
  if (!any(comunes)) stop("optimizar(): ningun punto evaluable comun en el tramo de estimacion.")
  mse <- vapply(yhats, function(v) mean((y_est[comunes] - v[comunes])^2), numeric(1))
  curva <- rejilla
  curva$mse <- mse
  mejor <- which.min(mse)
  list(
    parametros = c(fijos, as.list(rejilla[mejor, , drop = FALSE])),
    mse = mse[mejor],
    curva = curva,
    n_evaluados = sum(comunes)
  )
}


origen_movil <- function(y, metodo, parametros = list(), n0, H) {
  if (!is.function(metodo)) stop("origen_movil(): `metodo` debe ser una funcion ajustar_*().")
  T <- length(y)
  if (n0 >= T) stop("origen_movil(): `n0` debe ser menor que T.")
  filas <- lapply(n0:(T - 1), function(o) {
    hmax <- min(H, T - o)
    ajuste <- do.call(metodo, c(list(y = y[seq_len(o)]), parametros))
    pron <- ajuste$pronosticar(hmax)
    data.frame(origen = o, h = seq_len(hmax), objetivo = o + seq_len(hmax),
               real = y[o + seq_len(hmax)], pronostico = pron)
  })
  out <- do.call(rbind, filas)
  out$error <- out$real - out$pronostico
  tibble::as_tibble(out)
}


medidas <- function(e, real) {
  if (length(e) != length(real)) stop("medidas(): `e` y `real` deben tener la misma longitud.")
  ok <- !is.na(e)
  e <- e[ok]; real <- real[ok]
  tibble::tibble(
    n = length(e),
    RMSE = sqrt(mean(e^2)),
    MAD = mean(abs(e)),
    MAPE = 100 * mean(abs(e / real))
  )
}


validar_errores <- function(e, m, alfa = 0.05) {
  if (!is.numeric(e)) stop("validar_errores(): `e` debe ser numerico.")
  e <- e[!is.na(e)]
  n <- length(e)
  if (n < 8) stop("validar_errores(): muy pocos errores.")

  # 1) t de media cero
  t_obs <- mean(e) / (stats::sd(e) / sqrt(n))
  t_crit <- stats::qt(1 - alfa / 2, n - 1)
  p_t <- 2 * stats::pt(abs(t_obs), n - 1, lower.tail = FALSE)
  prueba_t <- prueba_seis(
    nombre = "t de media cero de los errores",
    h0 = "mu_e = 0", h1 = "mu_e != 0",
    estadistico = "t = ebar / (s_e / sqrt(n))",
    distribucion = paste0("t de Student con ", n - 1, " g.l."),
    observado = t_obs, critico = c(-t_crit, t_crit),
    region = paste0("|t| > ", formatC(t_crit, digits = 3, format = "f")),
    p = p_t, alfa = alfa,
    lectura_rechazo = "El metodo tiene sesgo: sus errores no estan centrados en cero.",
    lectura_no_rechazo = "No hay evidencia de sesgo en los pronosticos a un paso."
  )

  # 2) Ljung-Box, chi^2_m (no se restan grados de libertad)
  prueba_lb <- ljung_box(e, m, g = 0, alfa = alfa, que = "los errores a un paso")

  # 3) Jarque-Bera
  ec <- e - mean(e)
  s2 <- mean(ec^2)
  asim <- mean(ec^3) / s2^1.5
  curt <- mean(ec^4) / s2^2
  jb <- n / 6 * (asim^2 + (curt - 3)^2 / 4)
  jb_crit <- stats::qchisq(1 - alfa, 2)
  prueba_jb <- prueba_seis(
    nombre = "Jarque-Bera de los errores",
    h0 = "asimetria = 0 y curtosis = 3 (normalidad)",
    h1 = "asimetria != 0 o curtosis != 3",
    estadistico = "JB = n/6 [A^2 + (K - 3)^2/4]",
    distribucion = "chi-cuadrado con 2 g.l. (asintotica)",
    observado = jb, critico = jb_crit,
    region = paste0("JB > ", formatC(jb_crit, digits = 3, format = "f")),
    p = stats::pchisq(jb, 2, lower.tail = FALSE), alfa = alfa,
    lectura_rechazo = "Los errores no son normales: la banda +/-1,96 sd es solo heuristica.",
    lectura_no_rechazo = "No se rechaza normalidad de los errores."
  )

  # 4) Durbin-Watson (diagnostico)
  dw <- sum(diff_simple(e)^2) / sum(e^2)
  r1 <- acf_muestral(e, 1)

  list(
    n = n,
    t = prueba_t,
    ljung_box = prueba_lb,
    jarque_bera = prueba_jb,
    durbin_watson = c(d = dw, r1 = r1, aprox_2_1_menos_r1 = 2 * (1 - r1)),
    acf = acf_muestral(e, m)
  )
}

# diferencia simple sin diff(), para Durbin-Watson
diff_simple <- function(x) x[-1] - x[-length(x)]
