
prueba_seis <- function(nombre, h0, h1, estadistico, distribucion,
                        observado, critico, region, p, alfa = 0.05,
                        lectura_rechazo = "", lectura_no_rechazo = "") {
  rechaza <- p < alfa
  structure(list(
    nombre = nombre,
    hipotesis = c(H0 = h0, H1 = h1),
    estadistico = estadistico,
    distribucion = distribucion,
    region_rechazo = region,
    critico = critico,
    observado = observado,
    p_valor = p,
    decision = if (rechaza) "Se rechaza H0 al 5 %" else "No se rechaza H0 al 5 %",
    lectura = if (rechaza) lectura_rechazo else lectura_no_rechazo
  ), class = "prueba_seis")
}


print.prueba_seis <- function(x, ...) {
  cat("Prueba:", x$nombre, "\n")
  cat("1. Hipotesis:  H0:", x$hipotesis[["H0"]], "   vs   H1:", x$hipotesis[["H1"]], "\n")
  cat("2. Estadistico:", x$estadistico, " ~ ", x$distribucion, " bajo H0\n")
  cat("3. Region de rechazo (5 %):", x$region_rechazo, "\n")
  cat("4. Observado =", formatC(x$observado, digits = 4, format = "f"),
      "; valor p =", format.pval(x$p_valor, digits = 4), "\n")
  cat("5. Decision:", x$decision, "\n")
  cat("6. Lectura:", x$lectura, "\n")
  invisible(x)
}


autocovarianzas <- function(y, m) {
  if (!is.numeric(y) || anyNA(y)) stop("autocovarianzas(): `y` debe ser numerico y sin NA.")
  n <- length(y)
  if (m < 0 || m >= n) stop("autocovarianzas(): `m` debe estar entre 0 y n - 1.")
  yc <- y - mean(y)
  vapply(0:m, function(k) sum(yc[(k + 1):n] * yc[1:(n - k)]) / n, numeric(1))
}
#autocorrelacion muestral
acf_muestral <- function(y, m) {
  if (!is.numeric(m) || length(m) != 1L || m < 1) stop("acf_muestral(): `m` debe ser un entero >= 1.")
  g <- autocovarianzas(y, m)
  g[-1] / g[1]
}

#autocorrelacion parcial
pacf_dl <- function(r) {
  if (!is.numeric(r) || anyNA(r) || length(r) < 1L) stop("pacf_dl(): `r` debe ser numerico, sin NA.")
  m <- length(r)
  pacf <- numeric(m)
  phi <- numeric(0)                       # phi_{k-1, 1..k-1}
  for (k in seq_len(m)) {
    if (k == 1L) {
      phikk <- r[1]
      phi <- phikk
    } else {
      num <- r[k] - sum(phi * r[(k - 1):1])
      den <- 1 - sum(phi * r[1:(k - 1)])
      phikk <- num / den
      phi <- c(phi - phikk * rev(phi), phikk)
    }
    pacf[k] <- phikk
  }
  pacf
}


ljung_box <- function(y, m, g = 0, alfa = 0.05, que = "la serie") {
  if (!is.numeric(y) || anyNA(y)) stop("ljung_box(): `y` debe ser numerico y sin NA.")
  n <- length(y)
  if (m < 1 || m >= n) stop("ljung_box(): `m` debe estar entre 1 y n - 1.")
  if (m - g < 1) stop("ljung_box(): m - g debe ser >= 1.")
  r <- acf_muestral(y, m)
  Q <- n * (n + 2) * sum(r^2 / (n - seq_len(m)))
  gl <- m - g
  critico <- stats::qchisq(1 - alfa, gl)
  p <- stats::pchisq(Q, gl, lower.tail = FALSE)
  prueba_seis(
    nombre = paste0("Ljung-Box sobre ", que, " (m = ", m, ")"),
    h0 = paste0("rho(1) = ... = rho(", m, ") = 0"),
    h1 = paste0("rho(k) != 0 para algun k en 1, ..., ", m),
    estadistico = "Q = n(n+2) sum_{k=1}^{m} r_k^2/(n-k)",
    distribucion = paste0("chi-cuadrado con ", gl, " g.l."),
    observado = Q,
    critico = critico,
    region = paste0("Q > ", formatC(critico, digits = 3, format = "f")),
    p = p, alfa = alfa,
    lectura_rechazo = paste0("Hay autocorrelacion conjunta significativa en los primeros ", m,
                             " rezagos de ", que, "."),
    lectura_no_rechazo = paste0("No hay evidencia de autocorrelacion conjunta en los primeros ", m,
                                " rezagos de ", que, "; es compatible con ruido blanco.")
  )
}

#nucleo dde barlett
varianza_largo_plazo <- function(y, L = NULL) {
  if (!is.numeric(y) || anyNA(y)) stop("varianza_largo_plazo(): `y` debe ser numerico y sin NA.")
  n <- length(y)
  if (is.null(L)) L <- floor(4 * (n / 100)^(2 / 9))
  if (L < 0 || L >= n) stop("varianza_largo_plazo(): `L` debe estar entre 0 y n - 1.")
  g <- autocovarianzas(y, L)
  if (L == 0) return(g[1])
  pesos <- 1 - seq_len(L) / (L + 1)
  g[1] + 2 * sum(pesos * g[-1])
}
