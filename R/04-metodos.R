
validar_y <- function(y, min_n = 2, positiva = FALSE, quien = "metodo") {
  if (!is.numeric(y) || !is.null(dim(y))) stop(quien, "(): `y` debe ser un vector numerico.")
  if (anyNA(y)) stop(quien, "(): `y` no puede tener NA.")
  if (length(y) < min_n) stop(quien, "(): `y` necesita al menos ", min_n, " observaciones.")
  if (positiva && any(y <= 0)) stop(quien, "(): `y` debe ser estrictamente positiva.")
  invisible(TRUE)
}

validar_param <- function(x, nombre, quien, minimo = 0, maximo = 1, entero = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x >= minimo && x <= maximo
  if (entero) ok <- ok && x == round(x)
  if (!ok) stop(quien, "(): `", nombre, "` debe estar en [", minimo, ", ", maximo, "]",
                if (entero) " y ser entero" else "", ".")
  invisible(TRUE)
}

validar_h <- function(h) {
  if (!is.numeric(h) || length(h) != 1L || h < 1 || h != round(h)) {
    stop("pronosticar(): `h` debe ser un entero >= 1.")
  }
  invisible(TRUE)
}

#ajuste-minimos cuadrados
regresion_mco <- function(y, X, L = NULL) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (length(y) != n) stop("regresion_mco(): `y` y `X` tienen distinto numero de filas.")
  qrX <- qr(X)
  if (qrX$rank < p) stop("regresion_mco(): la matriz de diseno no tiene rango completo.")
  beta <- as.numeric(qr.coef(qrX, y))
  names(beta) <- colnames(X)
  ajustados <- as.numeric(X %*% beta)
  e <- y - ajustados
  XtX_inv <- chol2inv(qr.R(qrX))
  sigma2 <- sum(e^2) / (n - p)
  ee_ord <- sqrt(diag(XtX_inv) * sigma2)
  if (is.null(L)) L <- floor(4 * (n / 100)^(2 / 9))
  Xe <- X * e
  S <- crossprod(Xe)
  for (k in seq_len(L)) {
    w <- 1 - k / (L + 1)
    G <- crossprod(Xe[(k + 1):n, , drop = FALSE], Xe[1:(n - k), , drop = FALSE])
    S <- S + w * (G + t(G))
  }
  V_hac <- XtX_inv %*% S %*% XtX_inv
  ee_hac <- sqrt(diag(V_hac))
  list(coeficientes = beta, ee_ordinario = ee_ord, ee_robusto = ee_hac,
       residuos = e, ajustados = ajustados, sigma2 = sigma2, L = L, gl = n - p)
}


#' Media simple: 
ajustar_media <- function(y) {
  validar_y(y, 2, quien = "ajustar_media")
  n <- length(y)
  yhat <- c(NA_real_, cumsum(y)[-n] / seq_len(n - 1))
  media <- mean(y)
  pronosticar <- function(h) { validar_h(h); rep(media, h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "media", media = media))
}

#' Ingenuo
ajustar_ingenuo <- function(y) {
  validar_y(y, 2, quien = "ajustar_ingenuo")
  n <- length(y)
  yhat <- c(NA_real_, y[-n])
  ultimo <- y[n]
  pronosticar <- function(h) { validar_h(h); rep(ultimo, h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "ingenuo"))
}

#' Ingenuo 
ajustar_ingenuo_estacional <- function(y, s) {
  validar_param(s, "s", "ajustar_ingenuo_estacional", 2, Inf, entero = TRUE)
  validar_y(y, s + 1, quien = "ajustar_ingenuo_estacional")
  n <- length(y)
  yhat <- c(rep(NA_real_, s), y[1:(n - s)])
  ultimo_ciclo <- y[(n - s + 1):n]
  pronosticar <- function(h) {
    validar_h(h)
    ultimo_ciclo[((seq_len(h) - 1) %% s) + 1]
  }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "ingenuo estacional", s = s))
}

# ---- Promedios y suavizamientos -------------------------------------------------

#' Media movil de orden k
ajustar_mm <- function(y, k) {
  validar_param(k, "k", "ajustar_mm", 1, Inf, entero = TRUE)
  validar_y(y, k + 1, quien = "ajustar_mm")
  n <- length(y)
  cs <- c(0, cumsum(y))
  M <- rep(NA_real_, n)                         # M[t] = media de y[(t-k+1):t]
  M[k:n] <- (cs[(k + 1):(n + 1)] - cs[1:(n - k + 1)]) / k
  yhat <- c(NA_real_, M[-n])
  ultimo <- M[n]
  pronosticar <- function(h) { validar_h(h); rep(ultimo, h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "media movil", k = k))
}

#SES:Suavizamiento exponencial simple

ajustar_ses <- function(y, alpha) {
  validar_param(alpha, "alpha", "ajustar_ses", 0, 1)
  validar_y(y, 2, quien = "ajustar_ses")
  n <- length(y)
  L <- numeric(n)
  L[1] <- y[1]
  for (t in 2:n) L[t] <- alpha * y[t] + (1 - alpha) * L[t - 1]
  yhat <- c(NA_real_, L[-n])
  nivel <- L[n]
  pronosticar <- function(h) { validar_h(h); rep(nivel, h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "SES", alpha = alpha, nivel_final = nivel))
}

#Doble media movil de orden k.

ajustar_dmm <- function(y, k) {
  validar_param(k, "k", "ajustar_dmm", 2, Inf, entero = TRUE)
  validar_y(y, 2 * k, quien = "ajustar_dmm")
  n <- length(y)
  media_movil <- function(x, k) {
    out <- rep(NA_real_, length(x))
    for (t in k:length(x)) out[t] <- mean(x[(t - k + 1):t])
    out
  }
  M <- media_movil(y, k)
  M2 <- rep(NA_real_, n)
  for (t in (2 * k - 1):n) M2[t] <- mean(M[(t - k + 1):t])
  a <- 2 * M - M2
  b <- 2 * (M - M2) / (k - 1)
  yhat <- c(NA_real_, (a + b)[-n])            # NA hasta la posicion 2k - 1
  a_n <- a[n]; b_n <- b[n]
  pronosticar <- function(h) { validar_h(h); a_n + b_n * seq_len(h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "doble media movil", k = k, a = a_n, b = b_n))
}

# ---- Regresiones de tendencia --------------------------------------------------

disenar_tendencia <- function(t_idx, tipo) {
  switch(tipo,
    lineal = cbind(intercepto = 1, t = t_idx),
    cuadratica = cbind(intercepto = 1, t = t_idx, t2 = t_idx^2),
    exponencial = cbind(intercepto = 1, t = t_idx))
}

#Tendencia determinista lineal, cuadratica o exponencial por MCO.

ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial")) {
  tipo <- match.arg(tipo)
  validar_y(y, 5, positiva = (tipo == "exponencial"), quien = "ajustar_tendencia")
  n <- length(y)
  z <- if (tipo == "exponencial") log(y) else y
  p <- ncol(disenar_tendencia(1, tipo))
  yhat <- rep(NA_real_, n)
  for (t in (p + 2):n) {                       # ajuste con y[1..t-1], al menos p+1 datos
    X <- disenar_tendencia(seq_len(t - 1), tipo)
    b <- qr.coef(qr(X), z[1:(t - 1)])
    pred <- sum(disenar_tendencia(t, tipo) * b)
    yhat[t] <- if (tipo == "exponencial") exp(pred) else pred
  }
  X <- disenar_tendencia(seq_len(n), tipo)
  ajuste <- regresion_mco(z, X)
  b <- ajuste$coeficientes
  pronosticar <- function(h) {
    validar_h(h)
    pred <- as.numeric(disenar_tendencia(n + seq_len(h), tipo) %*% b)
    if (tipo == "exponencial") exp(pred) else pred
  }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = paste("tendencia", tipo), tipo = tipo, coeficientes = b))
}

disenar_tendencia_estacional <- function(t_idx, s, inicio) {
  est <- ((inicio - 1 + t_idx - 1) %% s) + 1           # estacion de cada t
  F <- matrix(0, length(t_idx), s - 1,
              dimnames = list(NULL, paste0("e", seq_len(s - 1))))
  for (j in seq_len(s - 1)) F[, j] <- (est == j) - (est == s)   # suma cero
  cbind(intercepto = 1, t = t_idx, F)
}


ajustar_tendencia_estacional <- function(y, s, inicio = 1) {
  validar_param(s, "s", "ajustar_tendencia_estacional", 2, Inf, entero = TRUE)
  validar_param(inicio, "inicio", "ajustar_tendencia_estacional", 1, s, entero = TRUE)
  p <- s + 1
  validar_y(y, p + 2, quien = "ajustar_tendencia_estacional")
  n <- length(y)
  yhat <- rep(NA_real_, n)
  primera <- max(p + 2, s + 2)                  # todas las estaciones vistas al menos una vez
  for (t in primera:n) {
    X <- disenar_tendencia_estacional(seq_len(t - 1), s, inicio)
    q <- qr(X)
    if (q$rank < ncol(X)) next
    b <- qr.coef(q, y[1:(t - 1)])
    yhat[t] <- sum(disenar_tendencia_estacional(t, s, inicio) * b)
  }
  ajuste <- regresion_mco(y, disenar_tendencia_estacional(seq_len(n), s, inicio))
  b <- ajuste$coeficientes
  efectos <- c(b[paste0("e", seq_len(s - 1))], -sum(b[paste0("e", seq_len(s - 1))]))
  names(efectos) <- paste0("estacion", seq_len(s))
  pronosticar <- function(h) {
    validar_h(h)
    as.numeric(disenar_tendencia_estacional(n + seq_len(h), s, inicio) %*% b)
  }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "tendencia con ficticias de suma cero", s = s, inicio = inicio,
                         coeficientes = b, efectos_estacionales = efectos))
}

# ---- Holt y Holt-Winters -------------------------------------------------------

#' Holt (tendencia lineal local).
ajustar_holt <- function(y, alpha, beta) {
  validar_param(alpha, "alpha", "ajustar_holt")
  validar_param(beta, "beta", "ajustar_holt")
  validar_y(y, 3, quien = "ajustar_holt")
  n <- length(y)
  L <- Tt <- numeric(n)
  L[1] <- y[1]; Tt[1] <- y[2] - y[1]
  yhat <- rep(NA_real_, n)
  for (t in 2:n) {
    if (t >= 3) yhat[t] <- L[t - 1] + Tt[t - 1]
    L[t] <- alpha * y[t] + (1 - alpha) * (L[t - 1] + Tt[t - 1])
    Tt[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Tt[t - 1]
  }
  L_n <- L[n]; T_n <- Tt[n]
  pronosticar <- function(h) { validar_h(h); L_n + T_n * seq_len(h) }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "Holt", alpha = alpha, beta = beta,
                         nivel_final = L_n, pendiente_final = T_n))
}

#' Holt-Winters aditivo con periodo s.

ajustar_hw <- function(y, s, alpha, beta, gamma) {
  validar_param(s, "s", "ajustar_hw", 2, Inf, entero = TRUE)
  validar_param(alpha, "alpha", "ajustar_hw")
  validar_param(beta, "beta", "ajustar_hw")
  validar_param(gamma, "gamma", "ajustar_hw")
  validar_y(y, 2 * s + 1, quien = "ajustar_hw")
  n <- length(y)
  L <- Tt <- S <- rep(NA_real_, n)
  L[s] <- mean(y[1:s])
  Tt[s] <- (mean(y[(s + 1):(2 * s)]) - mean(y[1:s])) / s
  S[1:s] <- y[1:s] - L[s]
  yhat <- rep(NA_real_, n)
  for (t in (s + 1):n) {
    pred <- L[t - 1] + Tt[t - 1] + S[t - s]
    if (t > 2 * s) yhat[t] <- pred
    L[t] <- alpha * (y[t] - S[t - s]) + (1 - alpha) * (L[t - 1] + Tt[t - 1])
    Tt[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Tt[t - 1]
    S[t] <- gamma * (y[t] - L[t]) + (1 - gamma) * S[t - s]
  }
  L_n <- L[n]; T_n <- Tt[n]; S_fin <- S[(n - s + 1):n]
  pronosticar <- function(h) {
    validar_h(h)
    k <- seq_len(h)
    L_n + T_n * k + S_fin[k - s * ceiling(k / s) + s]
  }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "Holt-Winters aditivo", s = s, alpha = alpha, beta = beta,
                         gamma = gamma, nivel_final = L_n, pendiente_final = T_n,
                         estacionales_finales = S_fin))
}

#' Holt-Winters multiplicativo con periodo s.

ajustar_hw_mult <- function(y, s, alpha, beta, gamma) {
  validar_param(s, "s", "ajustar_hw_mult", 2, Inf, entero = TRUE)
  validar_param(alpha, "alpha", "ajustar_hw_mult")
  validar_param(beta, "beta", "ajustar_hw_mult")
  validar_param(gamma, "gamma", "ajustar_hw_mult")
  validar_y(y, 2 * s + 1, positiva = TRUE, quien = "ajustar_hw_mult")
  n <- length(y)
  L <- Tt <- S <- rep(NA_real_, n)
  L[s] <- mean(y[1:s])
  Tt[s] <- (mean(y[(s + 1):(2 * s)]) - mean(y[1:s])) / s
  S[1:s] <- y[1:s] / L[s]
  yhat <- rep(NA_real_, n)
  for (t in (s + 1):n) {
    pred <- (L[t - 1] + Tt[t - 1]) * S[t - s]
    if (t > 2 * s) yhat[t] <- pred
    L[t] <- alpha * (y[t] / S[t - s]) + (1 - alpha) * (L[t - 1] + Tt[t - 1])
    Tt[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Tt[t - 1]
    S[t] <- gamma * (y[t] / L[t]) + (1 - gamma) * S[t - s]
  }
  L_n <- L[n]; T_n <- Tt[n]; S_fin <- S[(n - s + 1):n]
  pronosticar <- function(h) {
    validar_h(h)
    k <- seq_len(h)
    (L_n + T_n * k) * S_fin[k - s * ceiling(k / s) + s]
  }
  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(metodo = "Holt-Winters multiplicativo", s = s, alpha = alpha, beta = beta,
                         gamma = gamma, nivel_final = L_n, pendiente_final = T_n,
                         estacionales_finales = S_fin))
}
