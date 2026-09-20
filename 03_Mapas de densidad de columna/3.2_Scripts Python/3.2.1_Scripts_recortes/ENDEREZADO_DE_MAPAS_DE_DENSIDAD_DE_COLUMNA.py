"""
===============================================================================
ENDEREZADO GEOMÉTRICO DE MAPAS FITS DE DENSIDAD COLUMNAR (HERSCHEL/HOBYS)
===============================================================================

Objetivo
--------
Este script corrige la orientación geométrica de mosaicos FITS cuya región
observada aparece inclinada respecto a los ejes del array de píxeles.

Es importante destacar que NO se trata de una corrección astrométrica del WCS.
En estos mapas, el sistema de coordenadas celestes es correcto y no presenta
errores de orientación. La aparente rotación corresponde únicamente a la forma
geométrica del mosaico generado durante el escaneo observacional.

Por tanto, el procedimiento consiste en:

1. Identificar los píxeles válidos del mapa.
2. Determinar la orientación dominante del contorno observado.
3. Rotar la matriz de píxeles para alinear dicho contorno con los ejes.
4. Recortar automáticamente las zonas vacías introducidas por la rotación.
5. Guardar un nuevo mapa listo para análisis posteriores.

Aplicación científica
---------------------
Este preprocesado resulta especialmente útil antes de aplicar técnicas de
análisis multiescala como:

    - Funciones de distribución de probabilidad (PDF)
    - Espectros de potencias (Power Spectrum)
    - Delta Variance
    - Dendrogramas

ya que evita sesgos asociados a bordes inclinados y regiones vacías de gran
tamaño.

Requisitos
----------
pip install astropy scipy --break-system-packages
===============================================================================
"""

import os
import glob
import numpy as np

from astropy.io import fits
from scipy.ndimage import rotate as nd_rotate
from scipy.spatial import ConvexHull


# =============================================================================
# CONFIGURACIÓN GENERAL
# =============================================================================

# Directorio donde se encuentran los archivos FITS originales.
input_dir = r"E:\FITS"

# Directorio de salida para los mapas corregidos.
output_dir = os.path.join(input_dir, "enderezados")

# Creación automática de la carpeta de salida.
os.makedirs(output_dir, exist_ok=True)

# Umbral mínimo para considerar un píxel como válido.
#
# Se exige:
#   - que sea un número finito
#   - que supere este valor mínimo
#
# Si existen valores negativos físicamente significativos en el mapa,
# este parámetro debe ajustarse.
VALID_MIN = 0.0


# =============================================================================
# CONSTRUCCIÓN DE LA MÁSCARA DE DATOS VÁLIDOS
# =============================================================================

def build_valid_mask(data):
    """
    Genera una máscara booleana de píxeles observados.

    True  -> píxel con datos válidos
    False -> NaN o valor por debajo del umbral

    Esta máscara constituye la base geométrica utilizada durante todo
    el proceso de detección de orientación.
    """
    return np.isfinite(data) & (data > VALID_MIN)


# =============================================================================
# CÁLCULO DEL ÁNGULO MEDIANTE RECTÁNGULO DE ÁREA MÍNIMA
# =============================================================================

def min_area_rect_angle(mask):
    """
    Determina la orientación dominante de la región observada.

    Procedimiento:
    --------------
    1. Se extraen las coordenadas de todos los píxeles válidos.
    2. Se calcula la envolvente convexa (Convex Hull).
    3. Para cada lado de dicha envolvente se calcula un sistema de
       referencia rotado.
    4. Se evalúa el área del rectángulo envolvente.
    5. Se selecciona el ángulo que minimiza el área.

    Este método es robusto porque depende exclusivamente de la geometría
    global del mosaico y no de las estructuras físicas internas
    (filamentos, núcleos densos, etc.).

    Devuelve:
        Ángulo en grados.
    """

    ys, xs = np.where(mask)

    pts = np.column_stack([xs, ys]).astype(float)

    hull = ConvexHull(pts)

    hull_pts = pts[hull.vertices]

    best_angle = 0.0
    best_area = np.inf

    n = len(hull_pts)

    for i in range(n):

        p1 = hull_pts[i]
        p2 = hull_pts[(i + 1) % n]

        edge = p2 - p1

        ang = np.degrees(
            np.arctan2(edge[1], edge[0])
        )

        c = np.cos(np.radians(-ang))
        s = np.sin(np.radians(-ang))

        R = np.array([
            [c, -s],
            [s,  c]
        ])

        rot_pts = hull_pts @ R.T

        area = (
            np.ptp(rot_pts[:, 0]) *
            np.ptp(rot_pts[:, 1])
        )

        if area < best_area:
            best_area = area
            best_angle = ang

    return best_angle


# =============================================================================
# SELECCIÓN ÓPTIMA DEL SIGNO DE ROTACIÓN
# =============================================================================

def best_rotation_angle(mask, candidate):
    """
    Determina cuál de las soluciones geométricamente equivalentes produce
    la mejor alineación del mosaico.

    Se prueban varias combinaciones equivalentes:

        θ
        -θ
        θ + 90°
        θ - 90°
        ...

    El criterio utilizado consiste en maximizar el porcentaje de ocupación
    del rectángulo resultante.

    Cuanto mayor sea el llenado del bounding box, mejor alineado estará
    el paralelogramo observado respecto a los ejes del array.

    Devuelve:
        mejor_angulo, factor_de_llenado
    """

    candidates = [
        candidate,
        -candidate,
        candidate + 90,
        -candidate + 90,
        candidate - 90,
        -candidate - 90
    ]

    best_fill = -1
    best_ang = 0.0

    mask_f = mask.astype(np.float64)

    for ang in candidates:

        rot = nd_rotate(
            mask_f,
            angle=ang,
            reshape=True,
            order=1,
            cval=0.0
        )

        vm = rot > 0.5

        if not vm.any():
            continue

        rows = np.any(vm, axis=1)
        cols = np.any(vm, axis=0)

        r0, r1 = np.where(rows)[0][[0, -1]]
        c0, c1 = np.where(cols)[0][[0, -1]]

        fill = vm[r0:r1 + 1, c0:c1 + 1].mean()

        if fill > best_fill:
            best_fill = fill
            best_ang = ang

    return best_ang, best_fill


# =============================================================================
# ROTACIÓN CONSERVANDO LAS REGIONES NaN
# =============================================================================

def rotate_with_nan(data, mask, angle):
    """
    Rota el mapa preservando correctamente las zonas sin observación.

    La estrategia consiste en rotar de forma independiente:

        1. Los datos científicos.
        2. La máscara de validez.

    Posteriormente la máscara rotada se utiliza para restaurar los NaN.

    Esto evita que la interpolación introduzca valores artificiales en
    regiones originalmente no observadas.
    """

    data_filled = np.where(mask, data, 0.0)

    valid = mask.astype(np.float64)

    rot_data = nd_rotate(
        data_filled,
        angle=angle,
        reshape=True,
        order=1,
        cval=0.0
    )

    rot_valid = nd_rotate(
        valid,
        angle=angle,
        reshape=True,
        order=1,
        cval=0.0
    )

    vm = rot_valid > 0.5

    out = np.where(
        vm,
        rot_data,
        np.nan
    )

    return out, vm


# =============================================================================
# RECORTE AL ÁREA ÚTIL
# =============================================================================

def crop_to_valid(data, mask):
    """
    Elimina filas y columnas completamente vacías.

    Tras la rotación aparecen márgenes externos sin información.
    Estos márgenes no aportan contenido físico y aumentan
    innecesariamente el tamaño del archivo.

    Se recorta el mapa al rectángulo mínimo que contiene todos los
    píxeles válidos.
    """

    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)

    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]

    return data[r0:r1 + 1, c0:c1 + 1]


# =============================================================================
# BÚSQUEDA DE ARCHIVOS FITS
# =============================================================================

patterns = [
    "*.fits",
    "*.fit",
    "*.FITS"
]

fits_files = []

for p in patterns:
    fits_files.extend(
        glob.glob(os.path.join(input_dir, p))
    )

fits_files = sorted(set(fits_files))

print(
    f"Encontrados {len(fits_files)} archivos FITS en {input_dir}"
)


# =============================================================================
# PROCESAMIENTO PRINCIPAL
# =============================================================================

for path in fits_files:

    fname = os.path.basename(path)

    print(f"\nProcesando: {fname}")

    try:

        # ---------------------------------------------------------------------
        # Lectura del primer HDU con datos.
        # ---------------------------------------------------------------------

        with fits.open(path) as hdul:

            hdu = next(
                h for h in hdul
                if h.data is not None
            )

            data = hdu.data.astype(np.float64)

            header = hdu.header.copy()

        # ---------------------------------------------------------------------
        # Construcción de máscara.
        # ---------------------------------------------------------------------

        mask = build_valid_mask(data)

        if not mask.any():

            print(
                "  AVISO: no se encontraron píxeles válidos."
            )

            continue

        # ---------------------------------------------------------------------
        # Detección de orientación geométrica.
        # ---------------------------------------------------------------------

        candidate = min_area_rect_angle(mask)

        angle, fill = best_rotation_angle(
            mask,
            candidate
        )

        print(
            f"  ángulo candidato={candidate:.2f}° "
            f"-> ángulo aplicado={angle:.2f}° "
            f"(llenado={fill:.1%})"
        )

        # ---------------------------------------------------------------------
        # Rotación y recorte.
        # ---------------------------------------------------------------------

        rotated, rot_mask = rotate_with_nan(
            data,
            mask,
            angle
        )

        cropped = crop_to_valid(
            rotated,
            rot_mask
        )

        print(
            f"  forma original: {data.shape}"
            f" -> forma final: {cropped.shape}"
        )

        # ---------------------------------------------------------------------
        # Actualización del encabezado FITS.
        # ---------------------------------------------------------------------

        header_out = header.copy()

        header_out["NAXIS1"] = cropped.shape[1]
        header_out["NAXIS2"] = cropped.shape[0]

        header_out["HISTORY"] = (
            f"Rotación geométrica aplicada: "
            f"{angle:.3f} grados"
        )

        header_out["HISTORY"] = (
            "Escala de píxel conservada."
        )

        # ---------------------------------------------------------------------
        # Eliminación de claves WCS cuya precisión deja de estar garantizada.
        #
        # La escala física permanece inalterada, pero la correspondencia
        # exacta píxel-cielo deja de ser estrictamente válida tras una
        # rotación interpolada.
        #
        # Esto no afecta a PDF, espectros de potencias ni Delta Variance.
        # ---------------------------------------------------------------------

        for key in (
            "CRPIX1",
            "CRPIX2",
            "CRVAL1",
            "CRVAL2",
            "CROTA1",
            "CROTA2"
        ):

            if key in header_out:
                del header_out[key]

        # ---------------------------------------------------------------------
        # Escritura del mapa final.
        # ---------------------------------------------------------------------

        out_path = os.path.join(
            output_dir,
            fname
        )

        fits.writeto(
            out_path,
            cropped.astype(np.float32),
            header_out,
            overwrite=True
        )

        print(
            f"  guardado: {out_path}"
        )

    except Exception as e:

        print(
            f"  ERROR procesando {fname}: {e}"
        )

        continue


# =============================================================================
# FIN DEL PROCESO
# =============================================================================

print(
    "\nProceso completado."
)

print(
    "Todos los archivos enderezados se encuentran en:"
)

print(output_dir)