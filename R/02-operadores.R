
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
