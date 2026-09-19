# ==========================================================
# DELTA VARIANCE DE NUBES MOLECULARES
#
# Objetivo:
# Analizar la estructura multiescala de las nubes
# moleculares mediante el método Delta Variance.
#
# Se compara la pendiente β obtenida a partir de los
# mapas observados con los valores publicados por
# Schneider et al. (Tabla 4).
#
# Salidas:
#   - Gráfico JPG para cada región
#   - Excel multipestaña con resultados detallados
#   - Hoja resumen final
#
# Fundamentación teórica:
#
# La Delta Variance cuantifica cómo cambia la
# cantidad de estructura presente al variar la escala
# espacial de observación.
#
# Si:
#
#       σ²Δ(l) ∝ l^α
#
# entonces:
#
#       β = α + 2
#
# donde β es el índice del espectro de potencias.
#
# Referencias:
# Stutzki et al. (1998)
# Ossenkopf et al. (2008)
# Schneider et al. (2011, 2013)
#
# ==========================================================


# ==========================================================
# LIBRERÍAS
# ==========================================================

import os
import numpy as np

# Compatibilidad temporal con algunas versiones de
# TurbuStat que aún utilizan la nomenclatura antigua.
np.NaN = np.nan

import pandas as pd

from astropy.io import fits

# Lectura de información astrométrica del FITS.
from astropy.wcs import WCS

# Obtención de la escala angular real de cada píxel.
from astropy.wcs.utils import proj_plane_pixel_scales

# Implementación de Delta Variance.
from turbustat.statistics import DeltaVariance

# Ajuste lineal en espacio log-log.
from scipy.stats import linregress

import matplotlib.pyplot as plt


# ==========================================================
# DIRECTORIOS
# ==========================================================

# Carpeta que contiene los mapas FITS analizados.

INPUT_DIR = (
    r"E:\FITS\enderezados\recortes_Schneider"
)

# Carpeta de salida.

OUTPUT_DIR = (
    r"E:\FITS\enderezados\recortes_Schneider\estudio_Deltav"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# ==========================================================
# DISTANCIAS DE SCHNEIDER
# ==========================================================
#
# Distancia de cada región en kpc.
#
# Necesaria para convertir escalas físicas (pc)
# a escalas angulares observadas (arcmin).
#
# ==========================================================

DIST_KPC = {

    "dr15"   : 1.40,
    "m16"    : 2.00,
    "m17"    : 2.20,
    "monob1" : 0.80,
    "monr2"  : 0.862,
    "ngc2264": 0.719,
    "ngc6334": 1.35,
    "ngc6357": 1.75,
    "ngc7538": 2.80,
    "rosette": 1.46,
    "velac"  : 0.70

}


# ==========================================================
# ESCALA P1 DE SCHNEIDER
# ==========================================================
#
# P1 representa la escala característica donde la
# Delta Variance cambia de comportamiento.
#
# Schneider interpreta esta escala como una posible
# transición física dominante:
#
#   - filamentos
#   - regiones HII
#   - inyección de energía
#   - escalas gravitatorias
#
# ==========================================================

P1_PC = {

    "dr15"   : 1.42,
    "m16"    : 2.62,
    "m17"    : 2.57,
    "monob1" : 0.76,
    "monr2"  : 0.37,
    "ngc2264": 0.98,
    "ngc6334": 1.10,
    "ngc6357": 2.53,
    "ngc7538": 2.10,
    "rosette": 4.76,
    "velac"  : 1.83

}


# ==========================================================
# β PUBLICADOS POR SCHNEIDER
# ==========================================================
#
# Se utilizarán para comparar directamente con
# nuestros resultados.
#
# ==========================================================

BETA_SCHNEIDER = {

    "dr15"   : 2.17,
    "m16"    : 2.17,
    "m17"    : 2.22,
    "monob1" : 3.38,
    "monr2"  : 2.20,
    "ngc2264": 2.80,
    "ngc6334": 2.41,
    "ngc6357": 2.02,
    "ngc7538": 2.93,
    "rosette": 2.42,
    "velac"  : 2.30

}


# ==========================================================
# RESOLUCIÓN ANGULAR DE HERSCHEL
# ==========================================================
#
# Escala mínima físicamente fiable.
#
# Por debajo de ella domina el efecto instrumental.
#
# ==========================================================

BEAM_ARCSEC = 18.2

BEAM_ARCMIN = (
    BEAM_ARCSEC / 60.
)


# ==========================================================
# FUNCIONES AUXILIARES
# ==========================================================

def pixel_scale_arcmin(header):

    """
    Obtiene el tamaño angular de un píxel
    utilizando la información WCS del FITS.
    """

    w = WCS(header)

    scale_deg = (
        proj_plane_pixel_scales(w)[0]
    )

    return scale_deg * 60.


def pc_to_arcmin(pc, distance_kpc):

    """
    Convierte una escala física en pc
    a tamaño angular observado.

    Aproximación de ángulo pequeño:

        θ ≈ L / D
    """

    theta_rad = (
        pc /
        (distance_kpc * 1000.)
    )

    theta_arcmin = (
        np.rad2deg(theta_rad) * 60.
    )

    return theta_arcmin


# ==========================================================
# EXCEL GLOBAL
# ==========================================================

excel_file = os.path.join(
    OUTPUT_DIR,
    "DeltaVariance_Schneider.xlsx"
)

writer = pd.ExcelWriter(
    excel_file,
    engine="openpyxl"
)

summary_rows = []


# ==========================================================
# LISTA DE MAPAS FITS
# ==========================================================

fits_files = [

    f for f in os.listdir(INPUT_DIR)

    if f.lower().endswith(".fits")

]


# ==========================================================
# BUCLE PRINCIPAL
# ==========================================================

for file in fits_files:

    root = file.split("_")[0].lower()

    # Solo se procesan regiones con parámetros
    # definidos en Schneider.

    if root not in DIST_KPC:
        continue

    print()
    print("Procesando:", file)

    # ======================================================
    # CARGA DEL MAPA
    # ======================================================

    path = os.path.join(
        INPUT_DIR,
        file
    )

    hdu = fits.open(path)[0]

    data = np.squeeze(
        hdu.data
    ).astype(float)

    header = hdu.header


    # ======================================================
    # NORMALIZACIÓN
    # ======================================================
    #
    # Divide por la media para eliminar diferencias
    # de escala absoluta entre regiones.
    #
    # Δ-Variance solo estudia la distribución relativa
    # de estructuras.
    #
    # ======================================================

    data = data / np.nanmean(data)


    # ======================================================
    # CÁLCULO DE LA DELTA VARIANCE
    # ======================================================
    #
    # Se aplica un filtro ondícula ("Mexican Hat")
    # a distintas escalas espaciales.
    #
    # La varianza resultante indica cuánta estructura
    # existe en cada escala.
    #
    # ======================================================

    dv = DeltaVariance(
        (data, header)
    )

    dv.run(
        verbose=False
    )

    lags_pix = dv.lags.value

    delta_var = dv.delta_var

    delta_err = dv.delta_var_error


    # ======================================================
    # ESCALAS EN ARCMIN
    # ======================================================
    #
    # Conversión de píxel a tamaño angular real.
    #
    # ======================================================

    pixscale_arcmin = (
        pixel_scale_arcmin(header)
    )

    lag_arcmin = (
        lags_pix *
        pixscale_arcmin
    )


    # ======================================================
    # ESCALA P1 DE SCHNEIDER
    # ======================================================

    turnover_arcmin = pc_to_arcmin(

        P1_PC[root],

        DIST_KPC[root]

    )


    # ======================================================
    # SELECCIÓN DEL INTERVALO DE AJUSTE
    # ======================================================
    #
    # Límite inferior:
    #       resolución de Herschel
    #
    # Límite superior:
    #       escala P1
    #
    # Con ello reproducimos la metodología
    # de Schneider.
    #
    # ======================================================

    mask = (

        (lag_arcmin >= BEAM_ARCMIN)

        &

        (lag_arcmin <= turnover_arcmin)

    )


    # ======================================================
    # AJUSTE LINEAL EN LOG-LOG
    # ======================================================
    #
    # log(σ²Δ) = α log(l) + b
    #
    # ======================================================

    xfit = np.log10(
        lag_arcmin[mask]
    )

    yfit = np.log10(
        delta_var[mask]
    )

    fit = linregress(
        xfit,
        yfit
    )

    alpha = fit.slope

    beta_calc = alpha + 2

    r2 = fit.rvalue**2


    # ======================================================
    # EXPORTACIÓN A EXCEL
    # ======================================================

    df = pd.DataFrame({

        "Lag_arcmin": lag_arcmin,

        "DeltaVar": delta_var,

        "Error": delta_err

    })

    stats = pd.DataFrame({

        "Parametro": [

            "beta_schneider",
            "beta_calculado",
            "alpha",
            "R2",
            "P1_pc",
            "P1_arcmin"

        ],

        "Valor": [

            BETA_SCHNEIDER[root],
            beta_calc,
            alpha,
            r2,
            P1_PC[root],
            turnover_arcmin

        ]

    })


    # Cada nube se almacena en una hoja propia.

    sheet = root[:31]

    df.to_excel(
        writer,
        sheet_name=sheet,
        index=False
    )

    stats.to_excel(
        writer,
        sheet_name=sheet,
        startcol=5,
        index=False
    )


    # ======================================================
    # TABLA RESUMEN GLOBAL
    # ======================================================

    summary_rows.append({

        "Region": root,

        "Beta_Schneider":
            BETA_SCHNEIDER[root],

        "Beta_Calculado":
            beta_calc,

        "Alpha":
            alpha,

        "R2":
            r2,

        "P1_pc":
            P1_PC[root]

    })


    # ======================================================
    # REPRESENTACIÓN GRÁFICA
    # ======================================================
    #
    # Negro:
    #      Δ-Variance observada
    #
    # Gris:
    #      resolución instrumental
    #
    # Rojo:
    #      escala P1
    #
    # Azul:
    #      ajuste utilizado para obtener β
    #
    # ======================================================

    plt.figure(figsize=(8, 6))

    