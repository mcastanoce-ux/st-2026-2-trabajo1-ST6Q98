
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

#ejercicio6.c,d,e

texto_parametros <- function(par) {
  nom <- intersect(c("k", "alpha", "beta", "gamma"), names(par))
  if (length(nom) == 0L) return("sin constantes")
  paste0(nom, " = ", vapply(nom, function(n) format(round(par[[n]], 3), trim = TRUE), character(1)),
         collapse = ", ")
}


evaluar_metodo <- function(y, metodo, n0, H, rejilla = NULL, fijos = list()) {
  if (!is.function(metodo)) stop("evaluar_metodo(): `metodo` debe ser una funcion ajustar_*().")
  if (!is.null(rejilla)) {
    opt <- optimizar(y, metodo, rejilla, n0, fijos = fijos)
    par <- opt$parametros; mse_est <- opt$mse; curva <- opt$curva
  } else {
    par <- fijos; mse_est <- NA_real_; curva <- NULL
  }
  om <- origen_movil(y, metodo, par, n0, H)
  e1 <- om[om$h == 1, ]
  eH <- om[om$h == H, ]
  list(parametros = par, mse_est = mse_est, curva = curva, om = om,
       m1 = medidas(e1$error, e1$real), mH = medidas(eH$error, eH$real))
}

#' Resumen de los errores de origen movil en cada horizonte h = 1..H:
#' numero de errores, media, desviacion estandar y RMSE.
errores_por_horizonte <- function(om) {
  hs <- sort(unique(om$h))
  tibble::tibble(
    h = hs,
    n = vapply(hs, function(k) sum(om$h == k), numeric(1)),
    media = vapply(hs, function(k) mean(om$error[om$h == k]), numeric(1)),
    sd = vapply(hs, function(k) stats::sd(om$error[om$h == k]), numeric(1)),
    RMSE = vapply(hs, function(k) sqrt(mean(om$error[om$h == k]^2)), numeric(1))
  )
}

#' Curva (1 constante) o mapa (2 o 3 constantes) del MSE de un paso en el tramo de
#' estimacion, con el optimo marcado.
#'
#' @param curva data.frame devuelto por optimizar(): una columna por constante y `mse`.
#' @param titulo titulo del grafico.
grafico_mse <- function(curva, titulo) {
  pars <- setdiff(names(curva), "mse")
  if (length(pars) < 1L || length(pars) > 3L) stop("grafico_mse(): se esperan 1 a 3 constantes.")
  pie <- "MSE de un paso en el tramo de estimacion y[1..n0]; el marcador senala el minimo."
  if (length(pars) == 1L) {
    d <- data.frame(x = curva[[pars[1]]], mse = curva$mse)
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = x, y = mse)) +
        ggplot2::geom_line(colour = "grey45") +
        ggplot2::geom_point(size = 1.4, colour = "grey25") +
        ggplot2::geom_point(data = d[which.min(d$mse), ], colour = "firebrick", size = 3.2) +
        ggplot2::labs(title = titulo, x = pars[1], y = "MSE de un paso", caption = pie) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }
  d <- data.frame(a = curva[[pars[1]]], b = curva[[pars[2]]], mse = curva$mse)
  if (length(pars) == 3L) d$cc <- curva[[pars[3]]]
  mejor <- d[which.min(d$mse), ]
  g <- ggplot2::ggplot(d, ggplot2::aes(x = a, y = b, fill = mse)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "MSE") +
    ggplot2::geom_point(data = mejor, ggplot2::aes(x = a, y = b), inherit.aes = FALSE,
                        shape = 4, size = 3, stroke = 1.5, colour = "white") +
    ggplot2::labs(title = titulo, x = pars[1], y = pars[2],
                  caption = paste0(pie, if (length(pars) == 3L) paste0(" Un panel por ", pars[3], ".") else "")) +
    ggplot2::theme_minimal(base_size = 10)
  if (length(pars) == 3L) {
    etq <- paste0(pars[3], " = ", format(sort(unique(d$cc))))
    names(etq) <- as.character(sort(unique(d$cc)))
    g <- g + ggplot2::facet_wrap(~ cc, labeller = ggplot2::as_labeller(etq))
  }
  g
}

#' Validacion del ganador sobre sus errores de origen movil a UN paso (6d).
#'
#' El correlograma usa m = min(floor(n/4), 24) rezagos y su banda se calcula sobre el
#' numero n de errores. Ljung-Box se refiere a chi^2_m (no se restan grados de libertad)
#' con m = min(m_lb, m), para no pedir mas rezagos de los que los n errores sostienen.
#'
#' @param om tibble de origen_movil() del ganador.
#' @param m_lb m de Ljung-Box de la serie en la seccion 5.
validar_ganador <- function(om, m_lb) {
  e1 <- om[om$h == 1, ]
  e <- e1$error
  n <- length(e)
  m <- min(floor(n / 4), 24)
  m_test <- max(2, min(m_lb, m))
  r <- acf_muestral(e, m)
  val <- validar_errores(e, m_test)
  fuera <- which(abs(r) > 1.96 / sqrt(n))
  m_sens <- sort(unique(c(max(2, floor(m_test / 2)), m_test, min(n - 2, 2 * m_test))))
  sens <- do.call(rbind, lapply(m_sens, function(mm) {
    t <- ljung_box(e, mm, g = 0)
    data.frame(m = mm, Q = t$observado, gl = mm, critico = t$critico[1], p = t$p_valor)
  }))
  list(e = e, n = n, objetivo = e1$objetivo, origen = e1$origen, m = m, m_test = m_test,
       r = r, banda = 1.96 / sqrt(n), fuera = fuera, validacion = val,
       sensibilidad = tibble::as_tibble(sens), media = mean(e), sd = stats::sd(e))
}

#' Lectura en palabras de los rezagos que quedan fuera de la banda (6d).
#' @param fuera rezagos fuera de la banda; s periodo estacional candidato.
lectura_rezagos <- function(fuera, s) {
  if (length(fuera) == 0L) {
    return("Ningun rezago de la ACF de los errores sale de la banda: no queda autocorrelacion que el metodo no aproveche.")
  }
  partes <- character(0)
  if (1 %in% fuera) partes <- c(partes, "el rezago 1 indica persistencia de corto plazo: el error de hoy anticipa el de manana, es decir, el metodo no sigue el nivel con la rapidez suficiente")
  mult <- fuera[fuera %% s == 0 & fuera >= s]
  if (length(mult) > 0) partes <- c(partes, paste0("los rezagos ", paste(mult, collapse = ", "), " son multiplos de s = ", s, ": queda estructura estacional sin capturar"))
  otros <- setdiff(fuera, c(1, mult))
  if (length(otros) > 0) partes <- c(partes, paste0("los rezagos ", paste(otros, collapse = ", "), " son sueltos; con ", "una banda al 95 % algunos pueden salir por azar, y solo se tratan como estructura si forman un patron"))
  paste0("Salen de la banda los rezagos ", paste(fuera, collapse = ", "), ": ", paste(partes, collapse = "; "), ".")
}

#' Paso entre fechas consecutivas ("day", "month" o "3 months").
paso_fecha <- function(fechas) {
  d <- as.numeric(fechas[-1] - fechas[-length(fechas)])
  md <- stats::median(d)
  if (md <= 1.5) return("day")
  if (md >= 80 && md <= 100) return("3 months")
  if (md >= 27 && md <= 32) return("month")
  stop("paso_fecha(): frecuencia de fechas no reconocida.")
}

#' Pronostico final (6e): el ganador reajustado sobre toda la serie, H pronosticos y
#' banda HEURISTICA pronostico +/- z * desviacion de los errores de origen movil a cada
#' horizonte. No es un intervalo de confianza.
pronostico_final <- function(y, fechas, metodo, parametros, H, om, z = 1.96) {
  ajuste <- do.call(metodo, c(list(y = y), parametros))
  fc <- ajuste$pronosticar(H)
  eh <- errores_por_horizonte(om)
  eh <- eh[match(seq_len(H), eh$h), ]
  fut <- seq(fechas[length(fechas)], by = paso_fecha(fechas), length.out = H + 1)[-1]
  tibble::tibble(h = seq_len(H), fecha = fut, pronostico = fc, sd_error = eh$sd,
                 n_errores = eh$n, inferior = fc - z * eh$sd, superior = fc + z * eh$sd)
}

#' Errores a un paso en el tiempo (izquierda) y su correlograma (derecha), 6d.
grafico_errores_6d <- function(fechas_obj, e, r, n, titulo, unidad) {
  df <- data.frame(fecha = fechas_obj, error = e)
  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = fecha, y = error)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_line(linewidth = 0.4, colour = "grey25") +
    ggplot2::geom_point(size = 1, colour = "grey25") +
    ggplot2::labs(title = paste0(titulo, ": errores de origen movil a un paso"),
                  x = "Fecha del dato evaluado", y = paste0("Error (", unidad, ")"),
                  caption = paste0("n = ", n, " errores (y real - pronostico)")) +
    ggplot2::theme_minimal(base_size = 10)
  g2 <- correlograma(r, n, "ACF de los errores a un paso", "ACF")
  patchwork::wrap_plots(g1, g2, nrow = 1)
}

#' Ultimo tramo de la serie, pronostico final y banda heuristica (6e).
grafico_pronostico <- function(fechas, y, pf, n_ultimo, titulo, unidad, nombre_metodo, umbral = NULL) {
  Tn <- length(y)
  ini <- max(1, Tn - n_ultimo + 1)
  hist <- data.frame(fecha = fechas[ini:Tn], valor = y[ini:Tn])
  g <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pf, ggplot2::aes(x = fecha, ymin = inferior, ymax = superior),
                         fill = "firebrick", alpha = 0.2) +
    ggplot2::geom_line(data = hist, ggplot2::aes(x = fecha, y = valor), colour = "grey25") +
    ggplot2::geom_line(data = pf, ggplot2::aes(x = fecha, y = pronostico), colour = "firebrick") +
    ggplot2::geom_point(data = pf, ggplot2::aes(x = fecha, y = pronostico), colour = "firebrick", size = 1.3)
  if (!is.null(umbral)) g <- g + ggplot2::geom_hline(yintercept = umbral, linetype = "dashed", colour = "steelblue")
  g + ggplot2::labs(
    title = paste0(titulo, ": pronostico final (", nombre_metodo, ")"),
    x = "Fecha", y = unidad,
    caption = paste0("Ultimas ", nrow(hist), " observaciones y ", nrow(pf), " pronosticos. Banda HEURISTICA: ",
                     "pronostico +/- 1,96 veces la desviacion de los errores de origen movil a cada horizonte; ",
                     "no es un intervalo de confianza.")) +
    ggplot2::theme_minimal(base_size = 10)
}
