
graficar_serie <- function(datos, titulo, unidad,
                           cortes_fecha = NULL, formato_fecha = NULL) {
  if (!all(c("fecha", "valor") %in% names(datos))) {
    stop("graficar_serie(): `datos` debe tener columnas fecha y valor.")
  }
  rango_dias <- as.numeric(diff(range(datos$fecha)))
  if (is.null(cortes_fecha)) {
    cortes_fecha <- if (rango_dias <= 400) "1 month" else if (rango_dias <= 6000) "1 year" else "5 years"
  }
  if (is.null(formato_fecha)) {
    formato_fecha <- if (rango_dias <= 400) "%b\n%Y" else "%Y"
  }
  ggplot2::ggplot(datos, ggplot2::aes(x = fecha, y = valor)) +
    ggplot2::geom_line(colour = "grey25", linewidth = 0.5) +
    ggplot2::geom_point(colour = "grey25", size = 0.6) +
    ggplot2::scale_x_date(date_breaks = cortes_fecha, date_labels = formato_fecha) +
    ggplot2::labs(
      title = titulo,
      x = "Fecha",
      y = unidad,
      caption = paste0("n = ", nrow(datos), " observaciones (",
                       format(min(datos$fecha)), " a ", format(max(datos$fecha)), ")")
    ) +
    ggplot2::theme_minimal(base_size = 11)
}


grafico_estacional <- function(datos, s, titulo, unidad, etiquetas = NULL) {
  if (!all(c("fecha", "valor") %in% names(datos))) {
    stop("grafico_estacional(): `datos` debe tener columnas fecha y valor.")
  }
  if (!is.numeric(s) || length(s) != 1L || s < 2) stop("grafico_estacional(): `s` debe ser un entero >= 2.")
  df <- tibble::tibble(
    estacion = estacion_de_fecha(datos$fecha, s),
    ciclo = factor(ciclo_de_fecha(datos$fecha, s)),
    valor = datos$valor
  )
  if (is.null(etiquetas)) {
    etiquetas <- switch(as.character(s),
      "12" = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"),
      "4"  = c("T1", "T2", "T3", "T4"),
      "7"  = c("Lun", "Mar", "Mie", "Jue", "Vie", "Sab", "Dom"),
      as.character(seq_len(s)))
  }
  n_ciclos <- nlevels(df$ciclo)
  ggplot2::ggplot(df, ggplot2::aes(x = estacion, y = valor, group = ciclo, colour = as.integer(ciclo))) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_x_continuous(breaks = seq_len(s), labels = etiquetas) +
    ggplot2::scale_colour_gradient(low = "grey80", high = "navy", name = "Ciclo",
                                   breaks = range(seq_len(n_ciclos)),
                                   labels = levels(df$ciclo)[range(seq_len(n_ciclos))]) +
    ggplot2::labs(title = titulo, x = "Estacion", y = unidad,
                  caption = paste0(n_ciclos, " ciclos (completos o parciales), s = ", s,
                                   ", n = ", nrow(df), " observaciones")) +
    ggplot2::theme_minimal(base_size = 11)
}


grafico_rezagos <- function(y, rezagos, titulo, unidad) {
  if (!is.numeric(y) || anyNA(y)) stop("grafico_rezagos(): `y` debe ser numerico y sin NA.")
  n <- length(y)
  if (any(rezagos < 1) || any(rezagos >= n)) stop("grafico_rezagos(): rezagos fuera de rango.")
  r <- acf_muestral(y, max(rezagos))
  df <- do.call(rbind, lapply(rezagos, function(k) {
    data.frame(
      rezago = paste0("k = ", k, "  (r = ", formatC(r[k], digits = 2, format = "f"), ")"),
      orden = k,
      y_t = y[(k + 1):n],
      y_tk = y[1:(n - k)]
    )
  }))
  df$rezago <- factor(df$rezago, levels = unique(df$rezago[order(df$orden)]))
  ggplot2::ggplot(df, ggplot2::aes(x = y_tk, y = y_t)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(size = 0.8, alpha = 0.7, colour = "grey20") +
    ggplot2::facet_wrap(~ rezago) +
    ggplot2::labs(title = titulo,
                  x = paste0("y(t - k) [", unidad, "]"),
                  y = paste0("y(t) [", unidad, "]"),
                  caption = paste0("n = ", n, " observaciones; r = autocorrelacion muestral al rezago k")) +
    ggplot2::theme_minimal(base_size = 11)
}


correlograma <- function(r, n, titulo, tipo = c("ACF", "PACF"), bartlett = TRUE) {
  tipo <- match.arg(tipo)
  if (!is.numeric(r) || anyNA(r)) stop("correlograma(): `r` debe ser numerico y sin NA.")
  if (!is.numeric(n) || n < 2) stop("correlograma(): `n` debe ser el numero de observaciones.")
  m <- length(r)
  df <- data.frame(rezago = seq_len(m), r = as.numeric(r))
  banda_rb <- 1.96 / sqrt(n)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = rezago, y = r)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_segment(ggplot2::aes(xend = rezago, yend = 0), linewidth = 0.8) +
    ggplot2::geom_hline(data = data.frame(b = c(-banda_rb, banda_rb)),
                        ggplot2::aes(yintercept = b, linetype = "Ruido blanco"), colour = "steelblue")
  if (tipo == "ACF" && bartlett) {
    suma_previa <- c(0, cumsum(df$r^2))[seq_len(m)]          # sum_{i<h} r_i^2
    df$bartlett <- 1.96 * sqrt((1 + 2 * suma_previa) / n)
    g <- g +
      ggplot2::geom_line(data = df, ggplot2::aes(y = bartlett, linetype = "Bartlett"), colour = "firebrick") +
      ggplot2::geom_line(data = df, ggplot2::aes(y = -bartlett, linetype = "Bartlett"), colour = "firebrick")
  }
  g +
    ggplot2::scale_linetype_manual(name = "Banda 95 %",
                                   values = c("Ruido blanco" = "dashed", "Bartlett" = "dotted")) +
    ggplot2::scale_x_continuous(breaks = seq_len(m)) +
    ggplot2::labs(title = titulo, x = "Rezago k",
                  y = if (tipo == "ACF") "Autocorrelacion r(k)" else "Autocorrelacion parcial",
                  caption = paste0("n = ", n, "; banda de ruido blanco +/- 1,96/sqrt(n) = +/- ",
                                   formatC(banda_rb, digits = 3, format = "f"))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}

# =============================================================================
# SECCION 3(b): descomposicion clasica a mano y evidencia de estacionalidad
# (pegar al final de R/01-graficos.R)
# No usa decompose(), stl(), filter() ni ninguna funcion de series de tiempo.
# =============================================================================

etiquetas_estacion <- function(s) {
  switch(as.character(s),
         "12" = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"),
         "4"  = c("T1", "T2", "T3", "T4"),
         "7"  = c("Lun", "Mar", "Mie", "Jue", "Vie", "Sab", "Dom"),
         as.character(seq_len(s)))
}


# Media movil centrada de orden s.
#  - s impar: promedio simple de s valores, centrado en t.
#  - s par:   media movil 2 x s, con peso 1/(2s) en los extremos y 1/s en el resto.
# Deja NA en los floor(s/2) primeros y ultimos periodos.
media_movil_centrada <- function(y, s) {
  if (!is.numeric(y) || anyNA(y)) stop("media_movil_centrada(): `y` debe ser numerico y sin NA.")
  if (!is.numeric(s) || length(s) != 1L || s != round(s) || s < 2) {
    stop("media_movil_centrada(): `s` debe ser un entero >= 2.")
  }
  n <- length(y)
  if (n < 2 * s) stop("media_movil_centrada(): se necesitan al menos dos ciclos completos (n >= 2s).")
  h <- s %/% 2
  w <- if (s %% 2 == 1) rep(1 / s, s) else c(0.5, rep(1, s - 1), 0.5) / s
  tendencia <- rep(NA_real_, n)
  for (t in (h + 1):(n - h)) tendencia[t] <- sum(w * y[(t - h):(t + h)])
  tendencia
}


# Descomposicion clasica. `estacion` es la estacion (1, ..., s) de cada observacion.
#  aditiva:         y = T + S + R   indices normalizados con suma cero
#  multiplicativa:  y = T * S * R   indices normalizados con promedio uno
descomposicion_clasica <- function(y, s, estacion, tipo = c("aditiva", "multiplicativa")) {
  tipo <- match.arg(tipo)
  if (!is.numeric(y) || anyNA(y)) stop("descomposicion_clasica(): `y` debe ser numerico y sin NA.")
  if (length(estacion) != length(y)) stop("descomposicion_clasica(): `estacion` debe tener la longitud de `y`.")
  if (tipo == "multiplicativa" && any(y <= 0)) stop("descomposicion_clasica(): la version multiplicativa exige y > 0.")
  
  tendencia <- media_movil_centrada(y, s)
  sin_tendencia <- if (tipo == "aditiva") y - tendencia else y / tendencia
  est_f <- factor(estacion, levels = seq_len(s))
  crudos <- as.numeric(tapply(sin_tendencia, est_f, mean, na.rm = TRUE))   # promedio por estacion
  indices <- if (tipo == "aditiva") crudos - mean(crudos) else crudos / mean(crudos)
  estacional <- indices[as.integer(est_f)]
  residuo <- if (tipo == "aditiva") y - tendencia - estacional else y / (tendencia * estacional)
  
  list(y = y, s = s, tipo = tipo, estacion = as.integer(est_f),
       tendencia = tendencia, sin_tendencia = sin_tendencia,
       indices_crudos = crudos, indices = indices,
       estacional = estacional, residuo = residuo)
}


# Prueba F de igualdad de medias estacionales sobre la serie sin tendencia
# (analisis de varianza de un factor, calculado a mano).
# Nota: supone errores independientes; con autocorrelacion el valor p es solo una guia.
prueba_estacionalidad <- function(desc, que = "la serie", alfa = 0.05) {
  z <- desc$sin_tendencia
  ok <- !is.na(z)
  z <- z[ok]
  e <- factor(desc$estacion[ok], levels = seq_len(desc$s))
  s <- desc$s
  n <- length(z)
  medias <- as.numeric(tapply(z, e, mean))
  nj <- as.numeric(table(e))
  sce <- sum(nj * (medias - mean(z))^2)                 # entre estaciones
  scd <- sum((z - medias[as.integer(e)])^2)             # dentro de estaciones
  gl1 <- s - 1
  gl2 <- n - s
  Fobs <- (sce / gl1) / (scd / gl2)
  critico <- stats::qf(1 - alfa, gl1, gl2)
  p <- stats::pf(Fobs, gl1, gl2, lower.tail = FALSE)
  prueba <- prueba_seis(
    nombre = paste0("F de estacionalidad sobre ", que, " (s = ", s, ", version ", desc$tipo, ")"),
    h0 = paste0("mu_1 = ... = mu_", s, " (la serie sin tendencia tiene la misma media en las ", s, " estaciones)"),
    h1 = "al menos una media estacional es distinta",
    estadistico = "F = [SCE/(s-1)] / [SCD/(n-s)]",
    distribucion = paste0("F con ", gl1, " y ", gl2, " g.l."),
    observado = Fobs, critico = critico,
    region = paste0("F > ", formatC(critico, digits = 3, format = "f")),
    p = p, alfa = alfa,
    lectura_rechazo = paste0("Las medias por estacion de ", que, " difieren de forma significativa: hay efecto estacional de periodo ", s, "."),
    lectura_no_rechazo = paste0("No hay evidencia de que las medias por estacion de ", que, " difieran.")
  )
  list(prueba = prueba, F = Fobs, gl1 = gl1, gl2 = gl2, p = p,
       R2_estacional = sce / (sce + scd), n = n)
}


# Evidencia para elegir entre aditiva y multiplicativa: amplitud y nivel de cada ciclo completo.
# Si la amplitud crece con el nivel, la estacionalidad es multiplicativa; si no, aditiva.
criterio_aditivo_multiplicativo <- function(y, s, estacion, ciclo) {
  if (length(ciclo) != length(y)) stop("criterio_aditivo_multiplicativo(): `ciclo` debe tener la longitud de `y`.")
  filas <- lapply(sort(unique(ciclo)), function(c) {
    idx <- which(ciclo == c)
    if (length(idx) < s) return(NULL)                  # solo ciclos completos
    data.frame(ciclo = c, nivel = mean(y[idx]), amplitud = max(y[idx]) - min(y[idx]))
  })
  tab <- do.call(rbind, filas)
  if (is.null(tab) || nrow(tab) < 4) stop("criterio_aditivo_multiplicativo(): se necesitan al menos 4 ciclos completos.")
  tab$amp_relativa <- tab$amplitud / tab$nivel
  k <- nrow(tab)
  ini <- 1:3
  fin <- (k - 2):k
  resumen <- c(
    ciclos_completos = k,
    cor_amplitud_nivel = stats::cor(tab$nivel, tab$amplitud),
    nivel_inicial = mean(tab$nivel[ini]), nivel_final = mean(tab$nivel[fin]),
    amplitud_inicial = mean(tab$amplitud[ini]), amplitud_final = mean(tab$amplitud[fin]),
    amp_relativa_inicial = mean(tab$amp_relativa[ini]), amp_relativa_final = mean(tab$amp_relativa[fin])
  )
  list(ciclos = tab, resumen = resumen)
}


grafico_amplitud <- function(ciclos, titulo, unidad) {
  x <- ciclos$nivel
  v <- ciclos$amplitud
  b <- sum((x - mean(x)) * (v - mean(v))) / sum((x - mean(x))^2)     # recta MCO a mano
  a <- mean(v) - b * mean(x)
  ggplot2::ggplot(ciclos, ggplot2::aes(x = nivel, y = amplitud)) +
    ggplot2::geom_abline(intercept = a, slope = b, colour = "firebrick", linetype = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = ciclo), vjust = -0.8, size = 3) +
    ggplot2::labs(title = titulo,
                  x = paste0("Nivel del ciclo: promedio [", unidad, "]"),
                  y = paste0("Amplitud del ciclo: maximo - minimo [", unidad, "]"),
                  caption = paste0(nrow(ciclos), " ciclos completos; r(amplitud, nivel) = ",
                                   formatC(stats::cor(x, v), digits = 2, format = "f"))) +
    ggplot2::theme_minimal(base_size = 11)
}


grafico_descomposicion <- function(desc, fechas, titulo, unidad) {
  n <- length(desc$y)
  niveles <- c("Serie", "Tendencia (media movil centrada)", "Componente estacional", "Residuo")
  df <- data.frame(
    fecha = rep(fechas, 4),
    componente = factor(rep(niveles, each = n), levels = niveles),
    valor = c(desc$y, desc$tendencia, desc$estacional, desc$residuo)
  )
  ref <- if (desc$tipo == "aditiva") 0 else 1
  refs <- data.frame(componente = factor(niveles[3:4], levels = niveles), ref = ref)
  ggplot2::ggplot(df, ggplot2::aes(x = fecha, y = valor)) +
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = ref), colour = "grey60", linetype = "dashed") +
    ggplot2::geom_line(colour = "grey25", linewidth = 0.5, na.rm = TRUE) +
    ggplot2::facet_grid(componente ~ ., scales = "free_y") +
    ggplot2::labs(title = titulo, x = "Fecha",
                  y = if (desc$tipo == "aditiva") unidad else paste0(unidad, " (serie y tendencia); razon (estacional y residuo)"),
                  caption = paste0("Descomposicion clasica ", desc$tipo, ", s = ", desc$s, ", n = ", n, " observaciones")) +
    ggplot2::theme_minimal(base_size = 11)
}


grafico_indices <- function(desc, etiquetas, titulo) {
  df <- data.frame(estacion = factor(etiquetas, levels = etiquetas), indice = desc$indices)
  ref <- if (desc$tipo == "aditiva") 0 else 1
  ggplot2::ggplot(df, ggplot2::aes(x = estacion, y = indice)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_hline(yintercept = ref, colour = "grey30") +
    ggplot2::labs(title = titulo, x = "Estacion",
                  y = if (desc$tipo == "aditiva") "Indice estacional (unidades de la serie)" else "Indice estacional (razon)",
                  caption = paste0("Indices normalizados (", desc$tipo, "): ",
                                   if (desc$tipo == "aditiva") "suman 0" else "promedian 1")) +
    ggplot2::theme_minimal(base_size = 11)
}


# Convierte una prueba de seis elementos en una tabla para el informe.
prueba_a_tabla <- function(x) {
  data.frame(
    Elemento = c("1. Hipotesis", "2. Estadistico y distribucion", "3. Region de rechazo (5 %)",
                 "4. Observado y valor p", "5. Decision", "6. Lectura"),
    Contenido = c(
      paste0("H0: ", x$hipotesis[["H0"]], ";  H1: ", x$hipotesis[["H1"]]),
      paste0(x$estadistico, "  ~  ", x$distribucion, " bajo H0"),
      x$region_rechazo,
      paste0("observado = ", formatC(x$observado, digits = 4, format = "f"),
             "; valor p = ", format.pval(x$p_valor, digits = 4)),
      x$decision,
      x$lectura),
    stringsAsFactors = FALSE
  )
}