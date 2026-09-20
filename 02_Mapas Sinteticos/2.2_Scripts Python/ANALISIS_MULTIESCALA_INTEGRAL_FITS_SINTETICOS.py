# ============================================================
# ANÁLISIS MULTIESCALA DE MAPAS eFBM
# ============================================================
#
# Autor: Antonio Nombela Luengo
# Trabajo Fin de Máster en Astrofísica
#
# OBJETIVO CIENTÍFICO
# ------------------------------------------------------------
# Este script implementa tres herramientas estadísticas
# complementarias para caracterizar mapas sintéticos eFBM
# (enhanced Fractional Brownian Motion):
#
#   1) Power Spectrum
#   2) Probability Density Function (PDF)
#   3) Delta-Variance
#
# El objetivo es estudiar la influencia simultánea de:
#
#   • H  : exponente de Hurst
#   • Mach : número de Mach turbulento
#
# sobre la morfología y la estructura multiescala de los mapas.
#
# JUSTIFICACIÓN
# ------------------------------------------------------------
# Los campos tipo fBM constituyen modelos ampliamente utilizados
# para representar estructuras fractales autosimilares.
#
# La inclusión de un parámetro turbulento adicional (Mach)
# permite estudiar hasta qué punto las técnicas observacionales
# recuperan correctamente la información física originalmente
# impuesta durante la generación de los mapas.
#
# METODOLOGÍA GENERAL
# ------------------------------------------------------------
# Para cada mapa FITS:
#
# 1. Se calcula el espectro de potencia.
# 2. Se analiza la PDF de densidades.
# 3. Se calcula la Delta-Variance.
# 4. Se generan figuras comparativas.
# 5. Se exportan tablas compatibles con Excel.
#
# REFERENCIAS
# ------------------------------------------------------------
# Mandelbrot (1983)
# Stutzki et al. (1998)
# Federrath et al. (2010)
# Burkhart et al. (2013)
# Koch et al. (2017)
# ============================================================

# ============================================================
# IMPORTACIÓN DE LIBRERÍAS
# ============================================================
#
# Se utilizan librerías estándar de Python para:
#
#   - Gestión de archivos y directorios (os, glob, re)
#   - Cálculo numérico (numpy)
#   - Manipulación tabular de resultados (pandas)
#   - Visualización científica (matplotlib)
#
# Además, se emplean herramientas especializadas:
#
#   - Astropy:
#         Lectura y manipulación de archivos FITS.
#
#   - TurbuStat:
#         Implementación de estadísticas utilizadas en
#         estudios de estructura interestelar y turbulencia.
#
# Las tres técnicas principales utilizadas en este trabajo
# (Power Spectrum, PDF y Delta-Variance) se apoyan en estas
# librerías.
# ============================================================

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from turbustat.statistics import PowerSpectrum
from astropy.io import fits

# ============================================================
# LOCALIZACIÓN DE LOS MAPAS DE ENTRADA
# ============================================================
#
# Todos los mapas eFBM se encuentran almacenados en formato
# FITS dentro de un directorio común.
#
# Se seleccionan automáticamente aquellos archivos que siguen
# la nomenclatura:
#
#     efbm_Hx.xMachy.y_n.fits
#
# donde:
#
#     H     = exponente de Hurst
#     Mach  = número de Mach
#
# Esta estrategia permite automatizar completamente el análisis
# para grandes conjuntos de simulaciones.
# ============================================================

input_dir = r"E:\Astronomía\Master Astrofísica\ASGARD\Set 06 Jun 26"

fits_files = sorted(
    glob.glob(os.path.join(input_dir, "efbm*.fits"))
)

if len(fits_files) == 0:
    raise RuntimeError("No se encontraron archivos efbm*.fits")

# ============================================================
# EXTRACCIÓN AUTOMÁTICA DE PARÁMETROS FÍSICOS
# ============================================================
#
# El valor de H y Mach se recupera directamente del nombre
# del archivo mediante expresiones regulares.
#
# De esta forma se evita introducir manualmente los parámetros
# físicos y se minimiza la posibilidad de errores humanos.
#
# Ejemplo:
#
#     efbm_H0.5Mach4.0_3.fits
#
# devuelve:
#
#     H = 0.5
#     Mach = 4.0
# ============================================================

pattern = re.compile(
    r"efbm_H(?P<H>[0-9]*\.?[0-9]+)Mach(?P<Mach>[0-9]*\.?[0-9]+)",
    re.IGNORECASE
)

# ============================================================
# ANÁLISIS MEDIANTE POWER SPECTRUM
# ============================================================
#
# El espectro de potencia permite estudiar cómo se distribuye
# la energía o potencia estadística entre diferentes escalas
# espaciales.
#
# En representación log-log:
#
#     P(k) ∝ k^β
#
# donde:
#
#     k = número de onda espacial
#     β = pendiente espectral
#
# Para un campo fractal tipo fBM la teoría predice:
#
#     β = -(2H + 2)
#
# Por tanto, la comparación entre β teórico y β medido permite
# validar la recuperación del exponente de Hurst.
# ============================================================

rows = []

# ------------------------------------------------------------
# Se abre cada archivo FITS de manera individual.
#
# Para cada simulación se recuperan:
#
#     - Datos del mapa
#     - Cabecera FITS
#     - Parámetros H y Mach
#
# Los resultados obtenidos se almacenarán posteriormente en
# una tabla resumen.
# ------------------------------------------------------------

for ruta in fits_files:

    # ------------------------------------------------------------
    # Eliminación de la componente media.
    #
    # La frecuencia espacial k = 0 contiene únicamente información
    # sobre el valor medio global del mapa.
    #
    # Dado que el interés se centra en las fluctuaciones
    # estructurales, se resta la media antes de calcular el
    # espectro de potencia.
    # ------------------------------------------------------------

data = data - np.nanmean(data)