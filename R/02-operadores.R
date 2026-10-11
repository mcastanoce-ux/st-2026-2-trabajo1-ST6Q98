# =============================================================================
# R/02-operadores.R  --  Equipo ST6Q98
# Operador de rezago: diferenciar(), polinomio_rezago() y auxiliares de 4(b)
# =============================================================================

diferenciar <- function(y, d = 0, D = 0, s = 1) { 
  if (!is.numeric(y) || anyNA(y)) stop("diferenciar(): `y` debe ser numerico y sin NA.")
  es_entero <- function(x, minimo) is.numeric(x) && length(x) == 1L && x == round(x) && x >= minimo
  if (!es_entero(d, 0)) stop("diferenciar(): `d` debe ser un entero >= 0.")
  if (!es_entero(D, 0)) stop("diferenciar(): `D` debe ser un entero >= 0.")
  if (!es_entero(s, 1)) stop("diferenciar(): `s` debe ser un entero >= 1.")
  if (length(y) - d - D * s < 1) stop("diferenciar(): la serie es demasiado corta para esas diferencias.")
  z <- as.numeric(y)
  for (i in seq_len(D)) {
    n <- length(z)
    z <- z[(s + 1):n] - z[1:(n - s)]          # (1 - B^s) z_t = z_t - z_{t-s}
  }
  for (i in seq_len(d)) {
    n <- length(z)
    z <- z[2:n] - z[1:(n - 1)]                 # (1 - B) z_t = z_t - z_{t-1}
  }
  z
}

#Raices del polinomio
polinomio_rezago <- function(coef) {
  if (!is.numeric(coef) || length(coef) < 1L || anyNA(coef)) {
    stop("polinomio_rezago(): `coef` debe ser un vector numerico c(c1, ..., cp).")
  }
  
  while (length(coef) > 0L && coef[length(coef)] == 0) coef <- coef[-length(coef)]
  if (length(coef) == 0L) stop("polinomio_rezago(): el polinomio es constante (sin raices).")
  raices <- polyroot(c(1, coef))
  # limpieza numerica de partes ~0
  re <- Re(raices); im <- Im(raices)
  re[abs(re) < 1e-10] <- 0
  im[abs(im) < 1e-10] <- 0
  raices <- complex(real = re, imaginary = im)
  argumento <- Arg(raices)
  periodo <- ifelse(abs(argumento) < 1e-10, Inf, 2 * pi / abs(argumento))
  out <- data.frame(
    raiz = raices,
    real = re,
    imaginaria = im,
    modulo = Mod(raices),
    argumento = argumento,
    periodo = periodo
  )
  out[order(out$argumento), , drop = FALSE]
}


multiplicar_polinomios <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) stop("multiplicar_polinomios(): entradas numericas.")
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    idx <- i:(i + length(b) - 1L)
    out[idx] <- out[idx] + a[i] * b
  }
  out
}


# ---- Auxiliares de 4(b) ------------------------------------------------------

# Coeficientes completos c(1, c1, ..., cp) de (1 - B)^d (1 - B^s)^D.
# Para pasarlos a polinomio_rezago() hay que quitar el 1 inicial: coef[-1].
polinomio_diferencias <- function(d = 0, D = 0, s = 1) {
  p <- 1
  for (i in seq_len(d)) p <- multiplicar_polinomios(p, c(1, -1))
  for (i in seq_len(D)) p <- multiplicar_polinomios(p, c(1, rep(0, s - 1), -1))
  p
}

# Nombre del operador: tex = TRUE para el informe, FALSE para graficos.
etiqueta_operador <- function(d = 0, D = 0, s = 1, tex = TRUE) {
  pot <- function(k) if (tex) paste0("^{", k, "}") else paste0("^", k)
  partes <- character(0)
  if (d == 1) partes <- c(partes, "(1-B)")
  if (d > 1)  partes <- c(partes, paste0("(1-B)", pot(d)))
  if (D == 1) partes <- c(partes, paste0("(1-B", pot(s), ")"))
  if (D > 1)  partes <- c(partes, paste0("(1-B", pot(s), ")", pot(D)))
  if (length(partes) == 0) "1" else paste(partes, collapse = "")
}

# Polinomio expandido como texto LaTeX, por ejemplo: 1 - B - B^{4} + B^{5}.
polinomio_texto <- function(p) {
  out <- ""
  for (k in seq_along(p)) {
    cf <- p[k]
    if (cf == 0) next
    pot <- k - 1L
    mon <- if (pot == 0L) "" else if (pot == 1L) "B" else paste0("B^{", pot, "}")
    num <- if (abs(cf) == 1 && pot > 0L) "" else format(abs(cf))
    if (out == "") {
      out <- paste0(if (cf < 0) "-" else "", num, mon)
    } else {
      out <- paste(out, if (cf < 0) "-" else "+", paste0(num, mon))
    }
  }
  out
}

# Raices en el plano complejo con el circulo unidad.
# `raices` necesita las columnas candidata, real e imaginaria.
grafico_raices <- function(raices, titulo) {
  theta <- seq(0, 2 * pi, length.out = 200)
  circulo <- data.frame(x = cos(theta), y = sin(theta))
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_path(data = circulo, ggplot2::aes(x = x, y = y), colour = "grey55") +
    ggplot2::geom_point(data = raices, ggplot2::aes(x = real, y = imaginaria), size = 2.6) +
    ggplot2::facet_wrap(~candidata) +
    ggplot2::coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    ggplot2::labs(title = titulo, x = "Parte real", y = "Parte imaginaria",
                  caption = "Circulo = circulo unidad; una raiz doble se dibuja como un solo punto") +
    ggplot2::theme_minimal(base_size = 10)
}

# ---- Auxiliares de 4(c) ------------------------------------------------------

# Pendiente de la recta de minimos cuadrados de y contra t = 1, ..., n (formula cerrada).
pendiente_mco <- function(y) {
  if (!is.numeric(y) || anyNA(y)) stop("pendiente_mco(): `y` debe ser numerico y sin NA.")
  t <- seq_along(y)
  sum((t - mean(t)) * (y - mean(y))) / sum((t - mean(t))^2)
}

# Serie transformada: log opcional y luego (1-B)^d (1-B^s)^D.
serie_transformada <- function(y, escala = "ninguna", d = 0, D = 0, s = 1) {
  if (!is.numeric(y) || anyNA(y)) stop("serie_transformada(): `y` debe ser numerico y sin NA.")
  if (!escala %in% c("ninguna", "log")) stop("serie_transformada(): escala debe ser 'ninguna' o 'log'.")
  if (escala == "log") {
    if (any(y <= 0)) stop("serie_transformada(): log requiere y > 0.")
    y <- log(y)
  }
  diferenciar(y, d, D, s)
}

# Una fila de la tabla 4(c).
# media_esperada: lo que el operador debe dejar si la serie tiene tendencia lineal de pendiente b:
#   d + D = 0 -> el nivel (media de la serie);  (1-B) -> b;  (1-B^s) -> s*b;  d + D >= 2 -> 0.
resumen_transformacion <- function(y, escala = "ninguna", d = 0, D = 0, s = 1) {
  yb <- serie_transformada(y, escala)                 # solo la escala
  z  <- serie_transformada(y, escala, d, D, s)
  b  <- pendiente_mco(yb)
  esperada <- if (d + D == 0) mean(yb) else if (d + D >= 2) 0 else if (d == 1) b else s * b
  tibble::tibble(n = length(z), media = mean(z), varianza = stats::var(z),
                 r1 = acf_muestral(z, 1), media_esperada = esperada)
}

# Identidad Var(nabla Y) = 2 gamma(0) - 2 gamma(1) con las estimaciones muestrales (divisor n).
identidad_varianza <- function(y) {
  n  <- length(y)
  g0 <- mean((y - mean(y))^2)
  g1 <- sum((y[-1] - mean(y)) * (y[-n] - mean(y))) / n
  z  <- diferenciar(y, d = 1)
  izq <- mean((z - mean(z))^2)
  der <- 2 * g0 - 2 * g1
  tibble::tibble(gamma0 = g0, gamma1 = g1, lado_izq = izq, lado_der = der,
                 dif_relativa = abs(izq - der) / izq)
}