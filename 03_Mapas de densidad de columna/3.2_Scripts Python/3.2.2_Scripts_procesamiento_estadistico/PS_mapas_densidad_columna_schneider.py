# =====================================================================
# ANÁLISIS DEL ESPECTRO DE POTENCIAS DE MAPAS DE DENSIDAD DE COLUMNA
# =====================================================================
#
#
# Objetivo:
# Calcular el espectro de potencias 1D de mapas FITS de densidad de
# columna obtenidos con Herschel para estudiar la distribución de
# estructuras a diferentes escalas espaciales.
#
# El espectro de potencias permite analizar cómo se reparte la energía
# estructural de una nube molecular en función de la escala.
#
# Si:
#
#           P(l) ∝ l^β
#
# la pendiente β caracteriza la organización jerárquica de la nube.
#
# Se representan dos regímenes:
#
# β1 → escalas pequeñas/intermedias
# β2 → escalas grandes
#
# Los resultados se guardan en:
#
# 1. Excel con las tablas numéricas.
# 2. Imagen JPG con el espectro de potencias.
#
# =====================================================================

# =====================================================================
# IMPORTACIÓN DE LIBRERÍAS
# =====================================================================

import os                    # Manejo de archivos y directorios
import numpy as np           # Cálculo numérico
import pandas as pd          # Tablas y exportación a Excel

# ---------------------------------------------------------------------
# Compatibilidad TurbuStat / NumPy 2.x
# ---------------------------------------------------------------------
# Algunas versiones de TurbuStat utilizan np.NaN.
# NumPy 2.x elimina este alias.
# Se redefine para evitar errores.

np.NaN = np.nan

from astropy.io import fits
from turbustat.statistics import PowerSpectrum
import matplotlib.pyplot as plt


# =====================================================================
# DIRECTORIOS DE TRABAJO
# =====================================================================

# Carpeta que contiene los mapas FITS

input_dir = r"E:\FITS\enderezados\recortes_Schneider"

# Carpeta donde se guardarán los resultados

output_dir = (
    r"E:\FITS\enderezados\recortes_Schneider"
    r"\PowerSpectrum_Schneider"
)

# Crear carpeta automáticamente si no existe

os.makedirs(output_dir, exist_ok=True)

# Nombre del archivo Excel global

excel_file = os.path.join(
    output_dir,
    "PowerSpectrum_Todos.xlsx"
)


# =====================================================================
# RESOLUCIÓN ANGULAR DE HERSCHEL
# =====================================================================
#
# Beam efectivo de los mapas.
#
# La línea vertical correspondiente al beam marca la escala mínima
# confiable. Por debajo de ella las estructuras pueden estar afectadas
# por la respuesta instrumental.
#
# =====================================================================

beam_arcsec = 18.2
beam_arcmin = beam_arcsec / 60.0


# =====================================================================
# PARÁMETROS DE LOS AJUSTES
# =====================================================================
#
# Para cada región se almacenan:
#
# β1     Pendiente del régimen de pequeñas escalas
# lmin1  Escala mínima del ajuste β1
# lmax1  Escala máxima del ajuste β1
#
# β2     Pendiente del régimen de grandes escalas
# lmin2  Escala mínima del ajuste β2
# lmax2  Escala máxima del ajuste β2
#
# Estos valores proceden del análisis previo realizado para cada nube.
#
# =====================================================================

parametros = {

    "dr15_coldens_cf_high250_medsmo3": {
        "beta1": 4.80, "lmin1": 2, "lmax1": 20,
        "beta2": 2.56, "lmin2": 20, "lmax2": 100
    },

    "m16_coldens_cf_high250_medsmo3": {
        "beta1": 5.37, "lmin1": 8, "lmax1": 20,
        "beta2": 2.53, "lmin2": 20, "lmax2": 100
    },

    "m17_coldens_cf_high250_medsmo3": {
        "beta1": 3.54, "lmin1": 6, "lmax1": 20,
        "beta2": 2.37, "lmin2": 20, "lmax2": 200
    },

    "monob1_coldens_cf_high250_medsmo3": {
        "beta1": 5.43, "lmin1": 4, "lmax1": 20,
        "beta2": 2.98, "lmin2": 20, "lmax2": 100
    },

    "monr2_coldens_cf_high250_medsmo3": {
        "beta1": 5.94, "lmin1": 2, "lmax1": 20,
        "beta2": 2.33, "lmin2": 20, "lmax2": 100
    },

    "ngc2264_coldens_cf_high250_medsmo3": {
        "beta1": 6.11, "lmin1": 4, "lmax1": 20,
        "beta2": 2.79, "lmin2": 20, "lmax2": 100
    },

    "ngc6334_coldens_cf_high250_medsmo3": {
        "beta1": 2.68, "lmin1": 2, "lmax1": 6,
        "beta2": 0.09, "lmin2": 10, "lmax2": 100
    },

    "ngc6357_coldens_cf_high250_medsmo3": {
        "beta1": 4.95, "lmin1": 3, "lmax1": 20,
        "beta2": 2.19, "lmin2": 20, "lmax2": 100
    },

    "ngc7538_coldens_cf_high250_medsmo3": {
        "beta1": 5.29, "lmin1": 2, "lmax1": 20,
        "beta2": 2.28, "lmin2": 20, "lmax2": 200
    },

    "rosette_coldens_cf_high250_medsmo3": {
        "beta1": 4.04, "lmin1": 2, "lmax1": 20,
        "beta2": 2.51, "lmin2": 20, "lmax2": 200
    },

    "vela_coldens_cf_high250_medsmo3": {
        "beta1": 4.41, "lmin1": 2, "lmax1": 8,
        "beta2": 2.62, "lmin2": 20, "lmax2": 200
    }

}


# =====================================================================
# BÚSQUEDA DE ARCHIVOS FITS
# =====================================================================

fits_files = sorted([
    f for f in os.listdir(input_dir)
    if f.lower().endswith(".fits")
])

print()
print(f"Encontrados {len(fits_files)} archivos FITS")


# =====================================================================
# LISTA PARA ALMACENAR EL RESUMEN
# =====================================================================

resumen = []


# =====================================================================
# APERTURA DEL EXCEL GLOBAL
# =====================================================================

with pd.ExcelWriter(
    excel_file,
    engine="openpyxl"
) as writer:

    # ==============================================================
    # BUCLE PRINCIPAL
    # ==============================================================
    # Cada iteración procesa una nube molecular.
    # ==============================================================

    for fits_name in fits_files:

        fits_path = os.path.join(
            input_dir,
            fits_name
        )

        base_name = os.path.splitext(
            fits_name
        )[0]

        print()
        print("Procesando:", fits_name)

        try:

            # ======================================================
            # CARGA DEL MAPA FITS
            # ======================================================

            hdu = fits.open(fits_path)[0]

            data = np.squeeze(hdu.data).astype(float)

            header = hdu.header

            # ======================================================
            # NORMALIZACIÓN
            # ======================================================
            #
            # Divide el mapa por su valor medio:
            #
            # I_norm = I / <I>
            #
            # Esto elimina diferencias absolutas entre regiones y
            # permite comparar únicamente la estructura espacial.
            #
            # ======================================================

            data = data / np.nanmean(data)

            # ======================================================
            # CÁLCULO DEL ESPECTRO DE POTENCIAS
            # ======================================================
            #
            # TurbuStat calcula:
            #
            # P(k) = |F(k)|²
            #
            # donde F(k) es la transformada de Fourier del mapa.
            #
            # ======================================================

            pspec = PowerSpectrum((data, header))

            pspec.run(verbose=False)

            # ======================================================
            # INFORMACIÓN DIAGNÓSTICA
            # ======================================================

            print()
            print(base_name)
            print("Unidad:", pspec.freqs.unit)
            print("Frecuencia mínima:", np.nanmin(pspec.freqs))
            print("Frecuencia máxima:", np.nanmax(pspec.freqs))
            print("Primeras frecuencias:", pspec.freqs[:5])

            # ======================================================
            # EXTRACCIÓN DE RESULTADOS
            # ======================================================

            freqs = pspec.freqs.value
            power = pspec.ps1D

            # ======================================================
            # ELIMINACIÓN DE VALORES NO VÁLIDOS
            # ======================================================

            mask = (
                (freqs > 0)
                &
                (power > 0)
            )

            freqs = freqs[mask]
            power = power[mask]

            # ======================================================
            # CONVERSIÓN FRECUENCIA → ESCALA
            # ======================================================
            #
            # k = frecuencia espacial
            #
            # l = 1/k
            #
            # Es más intuitivo analizar escalas que frecuencias.
            #
            # ======================================================

            scale_arcmin = 1.0 / freqs

            # ======================================================
            # EXPORTACIÓN DE LA TABLA
            # ======================================================

            df = pd.DataFrame({
                "Scale_arcmin": scale_arcmin,
                "Frequency": freqs,
                "Power": power
            })

            sheet_name = base_name[:31]

            df.to_excel(
                writer,
                sheet_name=sheet_name,
                index=False
            )

            resumen.append({
                "Mapa": sheet_name,
                "Beam_arcmin": beam_arcmin
            })

            # ======================================================
            # REPRESENTACIÓN GRÁFICA
            # ======================================================

            plt.figure(figsize=(8, 6))

            plt.loglog(
                scale_arcmin,
                power,
                "ko",
                ms=4,
                label="Datos"
            )

            # ======================================================
            # AJUSTES β1 Y β2
            # ======================================================

            if base_name in parametros:

                p = parametros[base_name]

                # --------------------------------------------------
                # β1
                # --------------------------------------------------

                mask1 = (
                    (scale_arcmin >= p["lmin1"])
                    &
                    (scale_arcmin <= p["lmax1"])
                )

                xfit1 = scale_arcmin[mask1]
                yfit1 = power[mask1]

                if len(xfit1) > 3:

                    logA1 = np.mean(
                        np.log10(yfit1)
                        - p["beta1"] * np.log10(xfit1)
                    )

                    A1 = 10**logA1

                    xx1 = np.logspace(
                        np.log10(p["lmin1"]),
                        np.log10(p["lmax1"]),
                        300
                    )

                    yy1 = A1 * xx1**p["beta1"]

                    plt.loglog(
                        xx1,
                        yy1,
                        color="red",
                        lw=2,
                        label=f"β₁ = {p['beta1']:.2f}"
                    )

                # --------------------------------------------------
                # β2
                # --------------------------------------------------

                mask2 = (
                    (scale_arcmin >= p["lmin2"])
                    &
                    (scale_arcmin <= p["lmax2"])
                )

                xfit2 = scale_arcmin[mask2]
                yfit2 = power[mask2]

                if len(xfit2) > 3:

                    logA2 = np.mean(
                        np.log10(yfit2)
                        - p["beta2"] * np.log10(xfit2)
                    )

                    A2 = 10**logA2

                    xx2 = np.logspace(
                        np.log10(p["lmin2"]),
                        np.log10(p["lmax2"]),
                        300
                    )

                    yy2 = A2 * xx2**p["beta2"]

                    plt.loglog(
                        xx2,
                        yy2,
                        color="blue",
                        lw=2,
                        label=f"β₂ = {p['beta2']:.2f}"
                    )

                # --------------------------------------------------
                # Líneas delimitadoras de los ajustes
                # --------------------------------------------------

                plt.axvline(p["lmin1"],
                            color="red",
                            linestyle="--",
                            lw=1.8)

                plt.axvline(p["lmax1"],
                            color="green",
                            linestyle="--",
                            lw=1.8)

                plt.axvline(p["lmin2"],
                            color="blue",
                            linestyle="--",
                            lw=1.8)

                plt.axvline(p["lmax2"],
                            color="orange",
                            linestyle="--",
                            lw=1.8)

            # ======================================================
            # RESOLUCIÓN INSTRUMENTAL
            # ======================================================

            plt.axvline(
                beam_arcmin,
                color="black",
                ls=":",
                lw=2,
                label="Beam Herschel"
            )

            # ======================================================
            # CONFIGURACIÓN DEL GRÁFICO
            # ======================================================

            plt.xlabel("Escala espacial (arcmin)")
            plt.ylabel("Potencia")

            plt.title(base_name)

            plt.grid(alpha=0.3)

            plt.legend()

            # Escalas grandes a la izquierda
            plt.gca().invert_xaxis()

            plt.tight_layout()

            # ======================================================
            # GUARDADO DE LA FIGURA
            # ======================================================

            jpg_file = os.path.join(
                output_dir,
                f"{base_name}_PS.jpg"
            )

            plt.savefig(
                jpg_file,
                dpi=300
            )

            plt.close()

        except Exception as e:

            print(
                f"ERROR en {fits_name}: {e}"
            )

    # ==============================================================
    # HOJA RESUMEN FINAL
    # ==============================================================

    pd.DataFrame(resumen).to_excel(
        writer,
        sheet_name="Resumen",
        index=False
    )


# =====================================================================
# MENSAJE FINAL
# =====================================================================

print()
print("===================================================")
print("PROCESO FINALIZADO")
print("Archivo Excel:", excel_file)
print("===================================================")