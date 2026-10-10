
leer_serie <- function(ruta) {
  if (!is.character(ruta) || length(ruta) != 1L) {
    stop("leer_serie(): `ruta` debe ser una cadena de texto.")
  }
  if (!file.exists(ruta)) {
    stop("leer_serie(): no existe el archivo ", ruta)
  }
  datos <- readr::read_csv(
    ruta,
    col_types = readr::cols(fecha = readr::col_date(), valor = readr::col_double()),
    progress = FALSE
  )
  faltan <- setdiff(c("fecha", "valor"), names(datos))
  if (length(faltan) > 0L) {
    stop("leer_serie(): faltan columnas en ", ruta, ": ", paste(faltan, collapse = ", "))
  }
  if (anyNA(datos$fecha)) stop("leer_serie(): hay fechas que no se pudieron leer en ", ruta)
  if (anyNA(datos$valor)) stop("leer_serie(): hay valores faltantes en ", ruta)
  if (is.unsorted(datos$fecha, strictly = TRUE)) {
    stop("leer_serie(): las fechas de ", ruta, " no estan en orden estrictamente creciente.")
  }
  tibble::tibble(fecha = datos$fecha, valor = datos$valor)
}


estacion_de_fecha <- function(fechas, s) {
  if (!inherits(fechas, "Date")) stop("estacion_de_fecha(): `fechas` debe ser Date.")
  if (s == 12) return(as.integer(format(fechas, "%m")))
  if (s == 4) return((as.integer(format(fechas, "%m")) - 1L) %/% 3L + 1L)
  if (s == 7) return(as.integer(format(fechas, "%u")))
  ((seq_along(fechas) - 1L) %% s) + 1L
}


ciclo_de_fecha <- function(fechas, s) {
  if (s %in% c(4, 12)) return(as.integer(format(fechas, "%Y")))
  est <- estacion_de_fecha(fechas, s)
  # nuevo ciclo cada vez que la estacion vuelve a empezar
  cumsum(c(TRUE, diff(est) < 0))

ficha_serie <- function(datos, nombre, unidad, frecuencia) {
  if (!all(c("fecha", "valor") %in% names(datos))) {
    stop("ficha_serie(): `datos` debe tener columnas fecha y valor.")
  }
  tibble::tibble(
    nombre = nombre,
    unidad = unidad,
    frecuencia = frecuencia,
    fecha_inicial = min(datos$fecha),
    fecha_final = max(datos$fecha),
    T = nrow(datos),
    minimo = min(datos$valor),
    maximo = max(datos$valor),
    media = mean(datos$valor)
  )
}
