
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
