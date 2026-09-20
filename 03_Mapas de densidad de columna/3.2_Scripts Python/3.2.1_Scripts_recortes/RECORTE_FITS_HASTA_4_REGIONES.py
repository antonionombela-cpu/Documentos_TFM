"""
Recorte de subregiones físicas de un mapa de densidad columnar.

Objetivo científico:
--------------------
Este script extrae regiones de interés de un mapa FITS previamente
procesado y enderezado. Cada recorte representa una subestructura
del complejo molecular que posteriormente puede analizarse de forma
independiente mediante PDF, espectro de potencias, Delta Variance
o dendrogramas.

La selección se realiza mediante coordenadas de píxel definidas
manualmente a partir de la inspección visual del mapa completo.

Convención de coordenadas:
--------------------------
Las coordenadas se expresan en píxeles FITS estándar:

    origen = esquina inferior izquierda
    x crece hacia la derecha
    y crece hacia arriba

Cuando Astropy carga la imagen:

    data[y, x]

utiliza ya esta misma convención, por lo que no es necesario invertir
el eje vertical.

Además de los recortes FITS, el script genera una imagen JPG de control
de calidad donde aparecen marcadas las cajas seleccionadas. Esto permite
verificar visualmente que las regiones extraídas coinciden con las
estructuras deseadas.

Requiere:
    pip install astropy matplotlib --break-system-packages
"""

# ============================================================
# IMPORTACIÓN DE LIBRERÍAS
# ============================================================

# Gestión de rutas y directorios
import os

# Operaciones numéricas
import numpy as np

# Lectura y escritura de archivos FITS
from astropy.io import fits

# Visualización científica
import matplotlib.pyplot as plt

# Dibujado de rectángulos sobre imágenes
from matplotlib.patches import Rectangle

# ============================================================
# CONFIGURACIÓN
# ============================================================

# FITS de entrada previamente enderezado.
# Este mapa contiene la distribución bidimensional de densidad columnar.
input_path = r"E:\FITS\enderezados\w48_coldens_cf_high250_medsmo3.fits"

# Directorio donde se almacenarán los recortes individuales.
output_dir = r"E:\FITS\enderezados\crops_w48"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------
# Definición de subregiones
# ------------------------------------------------------------
#
# Cada región queda definida mediante dos esquinas opuestas:
#
#   (x_a, y_a, x_b, y_b)
#
# No importa el orden de las coordenadas porque posteriormente
# se calculan los valores mínimo y máximo de forma automática.
#
# Estas regiones pueden corresponder a:
#   - filamentos,
#   - cúmulos,
#   - núcleos densos,
#   - zonas activas de formación estelar,
#   - o regiones de comparación.
#
# El objetivo es aislar estructuras con propiedades físicas
# diferenciadas para analizarlas individualmente.

cajas_px = {

    # Región principal seleccionada.
    # Coordenadas expresadas en píxeles.
    1: (184, 2347, 2443, 148),

    # Ejemplos de regiones alternativas.
    # Se mantienen comentadas para futuros análisis.
    #
    # 2: (3646, 2230, 3943, 2516),
    # 3: (661, 1621, 1306, 848),
    # 4: (1562, 1434, 2792, 696),
}

# ============================================================
# CARGA DEL MAPA FITS
# ============================================================

with fits.open(input_path) as hdul:

    # Se selecciona automáticamente la primera extensión
    # que contenga una matriz de datos.
    hdu = next(h for h in hdul if h.data is not None)

    # Mapa científico.
    data = hdu.data

    # Cabecera FITS original.
    header = hdu.header

# Dimensiones del mapa.
ny, nx = data.shape

print(f"Imagen original: {data.shape} (filas x columnas)")

# ============================================================
# RECORTE DE SUBREGIONES
# ============================================================

for num, (xa, ya, xb, yb) in cajas_px.items():

    # --------------------------------------------------------
    # Ordenación automática de coordenadas
    # --------------------------------------------------------
    #
    # Permite introducir esquinas opuestas en cualquier orden.
    #
    x0, x1 = sorted((xa, xb))
    y0, y1 = sorted((ya, yb))

    # --------------------------------------------------------
    # Protección frente a coordenadas fuera de la imagen
    # --------------------------------------------------------
    #
    # Si alguna coordenada excede los límites del array,
    # se recorta automáticamente al rango válido.
    #
    x0c, x1c = max(0, x0), min(nx, x1)
    y0c, y1c = max(0, y0), min(ny, y1)

    if (x0c, x1c) != (x0, x1) or (y0c, y1c) != (y0, y1):

        print(
            f"  AVISO caja {num}: "
            f"rango solicitado fuera de la imagen; "
            f"ajustado a límites válidos."
        )

    print(
        f"Caja {num}: "
        f"x[{x0c}:{x1c}] "
        f"y[{y0c}:{y1c}] "
        f"→ tamaño {x1c - x0c} × {y1c - y0c}"
    )

    # --------------------------------------------------------
    # Extracción de la subimagen
    # --------------------------------------------------------
    #
    # Se obtiene únicamente la región rectangular definida.
    #
    recorte = data[y0c:y1c, x0c:x1c]

    # ========================================================
    # ACTUALIZACIÓN DEL HEADER FITS
    # ========================================================

    header_out = header.copy()

    # Nuevas dimensiones.
    header_out["NAXIS1"] = recorte.shape[1]
    header_out["NAXIS2"] = recorte.shape[0]

    # --------------------------------------------------------
    # Ajuste del sistema de referencia
    # --------------------------------------------------------
    #
    # Si existe un punto de referencia WCS (CRPIX),
    # debe desplazarse para mantener la coherencia
    # dentro del nuevo sistema de coordenadas.
    #
    # Esto preserva correctamente la localización
    # astronómica de los píxeles recortados.
    #
    if "CRPIX1" in header_out:
        header_out["CRPIX1"] = header_out["CRPIX1"] - x0c

    if "CRPIX2" in header_out:
        header_out["CRPIX2"] = header_out["CRPIX2"] - y0c

    # --------------------------------------------------------
    # Escritura del FITS recortado
    # --------------------------------------------------------

    out_path = os.path.join(
        output_dir,
        f"w48_caja{num}.fits"
    )

    fits.writeto(
        out_path,
        recorte,
        header_out,
        overwrite=True
    )

    print(f"  guardado: {out_path}")

# ============================================================
# GENERACIÓN DE IMAGEN DE VERIFICACIÓN
# ============================================================

# Para visualizar correctamente la imagen se utiliza un
# estiramiento robusto basado en percentiles.
#
# De esta forma los valores extremos no dominan el contraste.
#
vmin, vmax = np.nanpercentile(data, [1, 99])

fig, ax = plt.subplots(figsize=(10, 10))

# origin="lower" asegura que la representación visual emplee
# la misma convención de coordenadas utilizada en los recortes.
ax.imshow(
    data,
    origin="lower",
    cmap="gray",
    vmin=vmin,
    vmax=vmax
)

# ------------------------------------------------------------
# Dibujo de las cajas seleccionadas
# ------------------------------------------------------------
#
# Cada rectángulo representa una región recortada.
# Esto permite validar visualmente la selección antes
# de realizar cualquier análisis físico.
#
for num, (xa, ya, xb, yb) in cajas_px.items():

    x0, x1 = sorted((xa, xb))
    y0, y1 = sorted((ya, yb))

    rect = Rectangle(
        (x0, y0),
        x1 - x0,
        y1 - y0,
        linewidth=1.5,
        linestyle="--",
        edgecolor="yellow",
        facecolor="none"
    )

    ax.add_patch(rect)

    # Etiqueta identificativa.
    ax.text(
        x1,
        y1,
        f" {num}",
        color="yellow",
        fontsize=12,
        va="bottom",
        ha="left",
        weight="bold"
    )

# ------------------------------------------------------------
# Elementos gráficos
# ------------------------------------------------------------

ax.set_xlabel("x (píxeles)")
ax.set_ylabel("y (píxeles)")
ax.set_title("W48 — regiones seleccionadas para el análisis")

fig.tight_layout()

# ============================================================
# GUARDADO DE LA FIGURA DE CONTROL
# ============================================================

jpg_path = os.path.join(
    output_dir,
    "w48_cajas_verificacion.jpg"
)

fig.savefig(
    jpg_path,
    dpi=150,
    pil_kwargs={"quality": 90}
)

plt.close(fig)

print(f"\nImagen de verificación guardada: {jpg_path}")

# ============================================================
# FIN DEL PROCESO
# ============================================================

print("\nListo. Recortes almacenados en:")
print(output_dir)