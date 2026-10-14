https://github.com/mcastanoce-ux/st-2026-2-trabajo1-ST6Q98 · Equipo **ST6Q98**

# Trabajo 1 — Tres series: descripción, operador de rezago, dependencia, pronóstico y argumentación

Series de tiempo univariadas (3009297), Universidad Nacional de Colombia, sede Medellín. Semestre 2026-II. Prof. Juan Pablo Valencia Arango.

## Integrantes

| Integrante | Usuario de GitHub | Qué hizo |
|---|---|---|
| Mariana Castaño Celis | mcastanoce-ux | Escribió `R/00-lectura.R` y `R/05-evaluacion.R`. Terminó la sección 6 (pronóstico y evaluación) y escribió la sección 7 (conclusión técnica). Participó en `R/04-metodos.R`, `correr-todo.R` y el informe. |
| Sarah Catalina Criollo García | Sarah-Cr1ollo | Desarrolló los puntos 6(a) (métodos y asignación) y 6(b) (regresiones). Hizo el README y la exposición. Participó en `R/04-metodos.R`, `correr-todo.R` y el informe. |
| Sofia Cadavid Rueda | Sofia-Cadavid-Rueda | Escribió `R/01-graficos.R`, `R/02-operadores.R` y `R/03-dependencia.R`. Desarrolló las secciones 2 a 5 (series, descripción, operador de rezago y dependencia) e hizo la nota para los interesados (sección 8). Participó en `R/04-metodos.R`, `correr-todo.R` y el informe. |

## Qué contiene el repositorio

| Archivo | Contenido | Autor |
|---|---|---|
| `README.md` | Este documento | Sarah Catalina Criollo García |
| `.gitignore` | Excluye `.Rhistory`, `.RData` y `.Rproj.user` | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `st-2026-2-trabajo1-ST6Q98.Rproj` | Proyecto de RStudio | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `datos/ST6Q98-serie-1.csv` | Serie 1: tiempo medio diario de espera en la estación Los Almendros (segundos, diaria, T = 180), tal como se entregó | Profesor |
| `datos/ST6Q98-serie-2.csv` | Serie 2: pasajeros trimestrales del aeropuerto regional (miles, trimestral, T = 120), tal como se entregó | Profesor |
| `datos/ST6Q98-serie-3.csv` | Serie 3: viajes mensuales en la plataforma de movilidad (miles, mensual, T = 144), tal como se entregó | Profesor |
| `R/00-lectura.R` | `leer_serie()` | Mariana Castaño Celis |
| `R/01-graficos.R` | `graficar_serie()`, `grafico_estacional()`, `grafico_rezagos()`, `correlograma()` | Sofia Cadavid Rueda |
| `R/02-operadores.R` | `diferenciar()`, `polinomio_rezago()` | Sofia Cadavid Rueda |
| `R/03-dependencia.R` | `acf_muestral()`, `pacf_dl()`, `ljung_box()`, `varianza_largo_plazo()` | Sofia Cadavid Rueda |
| `R/04-metodos.R` | Los métodos del catálogo, uno por función (`ajustar_media()`, `ajustar_ingenuo()`, `ajustar_ingenuo_estacional()`, `ajustar_mm()`, `ajustar_ses()`, `ajustar_dmm()`, `ajustar_tendencia()`, `ajustar_tendencia_estacional()`, `ajustar_holt()`, `ajustar_hw()`, `ajustar_hw_mult()`) | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `R/05-evaluacion.R` | `optimizar()`, `origen_movil()`, `medidas()`, `validar_errores()` | Mariana Castaño Celis |
| `analisis/correr-todo.R` | Carga `R/`, verifica los datos y las funciones propias contra R, y corre las tres series en orden (secciones 2 a 6) | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `informe/informe.qmd` | Informe en Quarto | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `informe/informe.html` | Informe renderizado | Mariana Castaño Celis, Sarah Catalina Criollo García y Sofia Cadavid Rueda |
| `informe/nota-interesados.pdf` | Nota de una página para los interesados (sección 8) | Sofia Cadavid Rueda |
| `exposicion/ST6Q98.pdf` | Diapositivas de la exposición (sección 10) | Sarah Catalina Criollo García |
| `figs/` | Figuras generadas por `correr-todo.R`, con prefijo `serie1-`, `serie2-` y `serie3-`, y `resultados.rds` | Generado por código |
| `sesion-info.txt` | Salida de `sessionInfo()` de la máquina en la que se cerró la entrega | Generado por código |

### Funciones reutilizadas de la Tarea 1

| Archivo o función | De quién es | Qué se le cambió |
|---|---|---|
| No aplica | No aplica | No aplica |

No se reutilizaron funciones de la Tarea 1. Todas las funciones de `R/` se escribieron para este trabajo.

## Cómo se corre

Requisitos: R 4.3 o posterior, Quarto y los paquetes `dplyr`, `tidyr`, `purrr`, `tibble`, `readr`, `ggplot2` y `patchwork`. La entrega se cerró con R 4.5.2 en Windows 10 (ver `sesion-info.txt`).

```r
# 1. Instalar los paquetes (una sola vez)
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2", "patchwork"))

# 2. Desde la raíz del repositorio: correr todo el análisis
source("analisis/correr-todo.R")
```

```bash
# 3. Desde la raíz del repositorio: renderizar el informe
quarto render informe/informe.qmd
```

`correr-todo.R` usa rutas relativas, no pide intervención y guarda las figuras en `figs/`.

Tiempo de ejecución observado de `correr-todo.R`: 37,40 segundos (medido con `system.time(source("analisis/correr-todo.R"))`, valor `elapsed`).

## Resumen de resultados

RMSE de origen móvil, con los mismos orígenes para todos los métodos de cada serie. Constantes elegidas por MSE de un paso en el tramo de estimación (n0 = ⌊0,8 T⌋).

| Serie | Transformación elegida | Referente | Ganador a h = 1 (parámetros) | RMSE h = 1: ganador / referente | Mejor a h = H (parámetros) | RMSE h = H: mejor / referente |
|---|---|---|---|---|---|---|
| 1. Espera en estación (s), H = 9 | Ninguna, c(B) = 1 | Ingenuo | Holt–Winters aditivo (α = 0,9, β = 0,1, γ = 0,1) | 16,88 / 17,00 (−0,7 %) | Media simple (sin constantes) | 19,16 / 27,59 (−30,6 %) |
| 2. Pasajeros aeropuerto (miles), H = 4 | (1 − B)(1 − B⁴) | Ingenuo estacional (s = 4) | Holt–Winters aditivo (α = 0,7, β = 0,1, γ = 0,1) | 38,21 / 92,18 (−58,6 %) | Tendencia con ficticias de suma cero (sin constantes) | 41,30 / 90,35 (−54,3 %) |
| 3. Viajes plataforma (miles), H = 6 | log Yₜ y luego (1 − B) | Ingenuo | Tendencia exponencial (sin constantes) | 18,31 / 23,03 (−20,5 %) | Tendencia exponencial (sin constantes) | 18,99 / 37,37 (−49,2 %) |

## Declaración de uso de IA

El equipo usó Claude (Anthropic, vía claude.ai) en octubre de 2026 como apoyo para orientarse en la construcción de las estructuras del trabajo: cómo organizar los bloques de `correr-todo.R`, las funciones de `R/`, las secciones del informe y este README. La IA no tuvo acceso de escritura al repositorio. Se le compartieron los archivos, los datos, capturas de pantalla y las diapositivas de clase, y para el README leyó el repositorio público desde la web.

**Qué se pidió.** Orientación sobre la estructura de cada parte del trabajo. También se le pidió revisar las funciones propias (`diferenciar()`, `polinomio_rezago()` y los métodos de `R/04-metodos.R`), explicar cómo leer los gráficos de rezagos y los correlogramas, y resolver errores de renderizado de Quarto. Además, se le pidió un esqueleto del README con las seis secciones de 9(a).

**Qué se recibió.** Propuestas de estructura y de orden de los bloques, observaciones sobre el código y correcciones para seguir las convenciones de inicialización de las notas de clase (Holt, Holt–Winters y tendencia exponencial). Se recibieron también funciones auxiliares, como `polinomio_diferencias()` y `grafico_raices()`, el código de las regresiones de 6(b), borradores de texto para el informe y el esqueleto del README.

**Qué se verificó por cuenta propia.** Todo el código y el texto recibidos se copiaron, se corrieron y se revisaron a mano:

- `correr-todo.R` corre completo y pasa las verificaciones de las funciones propias contra R (`diff()`, `acf()`, `pacf()`, `Box.test()`, `lm()`).
- Las cifras de los borradores se contrastaron con la salida en consola.
- El informe se renderizó y se revisó que no quedaran fórmulas sin mostrar.
- Las cifras del resumen del README se cotejaron con el informe.

Las decisiones sobre las series las tomó el equipo con base en sus propios resultados: transformaciones, estacionalidad, elección de $m$, clasificación del patrón y asignación de métodos. Los autores y las tareas de cada integrante que aparecen en este README también los aportó el equipo.