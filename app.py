
"""
HidroSed · Generador de Cuenca Consolidada v1
---------------------------------------------

Objetivo:
Desde PC-HIDRO + PC-DESCARGA + eje del cauce + nombre de cuenca:
1) descarga o carga un DEM común;
2) delimita cuenca hidrológica y cuenca de descarga;
3) valida si las cuencas están truncadas por borde de DEM;
4) genera curvas de nivel desde DEM;
5) exporta una "Cuenca consolidada <nombre>.kmz".

Diseño:
- PC-HIDRO: punto de control hidrológico.
- PC-DESCARGA: punto de descarga/cierre de cuenca soporte.
- Eje: línea del cauce que conecta ambos puntos.
- Cuenca consolidada: contiene cuenca descarga, cuenca hidrológica, intercuenca, eje, puntos y curvas.
"""

from __future__ import annotations

import io
import json
import math
import re
import heapq
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

try:
    import requests
except Exception:
    requests = None

import rasterio
from rasterio.features import shapes, geometry_mask
from rasterio.transform import Affine
from rasterio.crs import CRS
from rasterio.io import MemoryFile
from rasterio.warp import calculate_default_transform, reproject, Resampling

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, MultiLineString, shape, mapping
from shapely.ops import unary_union, transform as shp_transform
from pyproj import CRS as PyCRS, Transformer

try:
    from skimage import measure
except Exception:
    measure = None


APP_TITLE = "HidroSed · Generador de Cuenca Consolidada v1.2 · Motor DEM COP30 probado"
OPENTOPO_URL = "https://portal.opentopography.org/API/globaldem"

st.set_page_config(page_title=APP_TITLE, page_icon="🌊", layout="wide")

D8_OFFSETS = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, -1),           (0, 1),
    (1, -1),  (1, 0),  (1, 1),
]


@dataclass
class ControlPoint:
    name: str
    lon: float
    lat: float
    x: float
    y: float


@dataclass
class BasinResult:
    name: str
    original_rc: Tuple[int, int]
    snapped_rc: Tuple[int, int]
    mask: np.ndarray
    geom: object
    area_km2: float
    snapped_x: float
    snapped_y: float
    snapped_lon: float
    snapped_lat: float
    accumulation_cells: float


# =============================================================================
# KML / KMZ
# =============================================================================

def extract_kml_from_upload(uploaded_file) -> str:
    raw = uploaded_file.getvalue()
    name = uploaded_file.name.lower()
    if name.endswith(".kmz"):
        with zipfile.ZipFile(io.BytesIO(raw), "r") as zf:
            kml_names = [n for n in zf.namelist() if n.lower().endswith(".kml")]
            if not kml_names:
                raise ValueError("El KMZ no contiene un KML interno.")
            preferred = next((n for n in kml_names if n.lower().endswith("doc.kml")), kml_names[0])
            return zf.read(preferred).decode("utf-8", errors="ignore")
    if name.endswith(".kml"):
        return raw.decode("utf-8", errors="ignore")
    raise ValueError("El archivo debe ser KML o KMZ.")


def _coord_tokens(text: str) -> List[Tuple[float, float]]:
    coords = []
    if not text:
        return coords
    raw = text.strip().replace("\n", " ").replace("\t", " ")
    for token in raw.split():
        parts = [p for p in token.split(",") if p != ""]
        if len(parts) >= 2:
            coords.append((float(parts[0]), float(parts[1])))
    return coords


def parse_first_point(uploaded_file) -> Tuple[float, float]:
    kml = extract_kml_from_upload(uploaded_file)
    root = ET.fromstring(kml.encode("utf-8"))
    ns = {"k": "http://www.opengis.net/kml/2.2"}
    for pnode in root.findall(".//k:Point", ns):
        cnode = pnode.find(".//k:coordinates", ns)
        if cnode is not None and cnode.text:
            coords = _coord_tokens(cnode.text)
            if coords:
                return coords[0]
    # Respaldo sin namespace
    for elem in root.iter():
        if elem.tag.lower().endswith("point"):
            for child in elem.iter():
                if child.tag.lower().endswith("coordinates") and child.text:
                    coords = _coord_tokens(child.text)
                    if coords:
                        return coords[0]
    raise ValueError("No se encontró un punto válido en el KML/KMZ.")


def parse_first_line(uploaded_file) -> List[Tuple[float, float]]:
    kml = extract_kml_from_upload(uploaded_file)
    root = ET.fromstring(kml.encode("utf-8"))
    ns = {"k": "http://www.opengis.net/kml/2.2"}
    best = []
    for lnode in root.findall(".//k:LineString", ns):
        cnode = lnode.find(".//k:coordinates", ns)
        if cnode is not None and cnode.text:
            coords = _coord_tokens(cnode.text)
            if len(coords) > len(best):
                best = coords
    if len(best) >= 2:
        return best

    for elem in root.iter():
        if elem.tag.lower().endswith("linestring"):
            for child in elem.iter():
                if child.tag.lower().endswith("coordinates") and child.text:
                    coords = _coord_tokens(child.text)
                    if len(coords) > len(best):
                        best = coords
    if len(best) >= 2:
        return best
    raise ValueError("No se encontró un eje/LineString válido en el KML/KMZ.")


def kml_escape(v: Any) -> str:
    return escape(str(v), quote=True)


# =============================================================================
# DEM download/read/reproject
# =============================================================================

def km_per_degree_lon(lat: float) -> float:
    return max(1e-6, 111.320 * math.cos(math.radians(lat)))


def bbox_from_lonlat_points(points: List[Tuple[float, float]], buffer_km: float) -> Tuple[float, float, float, float]:
    if not points:
        raise ValueError("No hay coordenadas para construir el área de descarga DEM.")
    lons = [p[0] for p in points]
    lats = [p[1] for p in points]
    midlat = (min(lats) + max(lats)) / 2.0
    dlat = buffer_km / 111.320
    dlon = buffer_km / km_per_degree_lon(midlat)
    south = min(lats) - dlat
    north = max(lats) + dlat
    west = min(lons) - dlon
    east = max(lons) + dlon
    return south, north, west, east


def bbox_area_km2(south: float, north: float, west: float, east: float) -> float:
    midlat = (south + north) / 2.0
    return abs((north - south) * 111.320) * abs((east - west) * km_per_degree_lon(midlat))


def mask_api_key(api_key: str) -> str:
    api_key = (api_key or "").strip()
    if not api_key:
        return ""
    if len(api_key) <= 8:
        return "****"
    return f"{api_key[:4]}...{api_key[-4:]}"


def margin_to_degrees(lat: float, margin_value: float, margin_unit: str) -> Tuple[float, float]:
    if margin_value <= 0:
        raise ValueError("El margen debe ser mayor que cero.")
    if margin_unit == "Grados decimales":
        return margin_value, margin_value
    lat_delta = margin_value / 111.32
    cos_lat = max(abs(math.cos(math.radians(lat))), 0.01)
    lon_delta = margin_value / (111.32 * cos_lat)
    return lat_delta, lon_delta


def bbox_from_center(lat: float, lon: float, margin_value: float, margin_unit: str) -> Tuple[float, float, float, float]:
    """BBox igual a la aplicación DEM COP30 probada: un punto central + margen."""
    lat_delta, lon_delta = margin_to_degrees(lat, margin_value, margin_unit)
    south = max(-90.0, lat - lat_delta)
    north = min(90.0, lat + lat_delta)
    west = max(-180.0, lon - lon_delta)
    east = min(180.0, lon + lon_delta)
    return round(south, 6), round(north, 6), round(west, 6), round(east, 6)


def expand_bbox_to_include_points(
    south: float, north: float, west: float, east: float,
    points_lonlat: List[Tuple[float, float]],
    padding_km: float = 2.0,
) -> Tuple[float, float, float, float]:
    """Amplía un bbox para incluir PC-DESCARGA/eje si quedan fuera, manteniendo la lógica de margen central."""
    if not points_lonlat:
        return south, north, west, east
    lons = [p[0] for p in points_lonlat]
    lats = [p[1] for p in points_lonlat]
    lat_mid = (min(lats) + max(lats)) / 2.0
    dlat = padding_km / 111.320
    dlon = padding_km / km_per_degree_lon(lat_mid)
    return (
        min(south, min(lats) - dlat),
        max(north, max(lats) + dlat),
        min(west, min(lons) - dlon),
        max(east, max(lons) + dlon),
    )


def looks_like_geotiff(content: bytes, content_type: str) -> bool:
    if not content:
        return False
    tiff_magic = content.startswith(b"II*\x00") or content.startswith(b"MM\x00*")
    binary_type = any(token in content_type.lower() for token in ["tiff", "geotiff", "octet-stream", "application/x-tiff"])
    html_or_json = content.lstrip().startswith((b"<", b"{", b"["))
    return (tiff_magic or binary_type) and not html_or_json


def safe_response_preview_text(resp, limit: int = 500) -> str:
    try:
        txt = resp.text or ""
    except Exception:
        txt = ""
    txt = re.sub(r"API_Key=[^&\s]+", "API_Key=****", txt)
    return txt[:limit]


def build_opentopo_url(demtype: str, south: float, north: float, west: float, east: float, api_key: str) -> str:
    params = {
        "demtype": demtype,
        "south": f"{south:.8f}",
        "north": f"{north:.8f}",
        "west": f"{west:.8f}",
        "east": f"{east:.8f}",
        "outputFormat": "GTiff",
        "API_Key": api_key.strip(),
    }
    return OPENTOPO_URL + "?" + urlencode(params)


def download_dem_opentopo(demtype: str, south: float, north: float, west: float, east: float, api_key: str) -> bytes:
    """Descarga DEM usando la lógica robusta de la app DEM COP30 que el usuario validó."""
    if requests is None:
        raise RuntimeError("No está instalado requests.")
    params = {
        "demtype": demtype,
        "south": round(float(south), 6),
        "north": round(float(north), 6),
        "west": round(float(west), 6),
        "east": round(float(east), 6),
        "outputFormat": "GTiff",
        "API_Key": api_key.strip(),
    }
    try:
        resp = requests.get(OPENTOPO_URL, params=params, timeout=(10, 180))
    except requests.Timeout:
        raise RuntimeError("La solicitud excedió el tiempo de espera. Reduzca el área o intente nuevamente.")
    except requests.RequestException as exc:
        raise RuntimeError(f"Error de conexión con OpenTopography: {exc}")

    if resp.status_code == 204:
        raise RuntimeError("204 No Data: OpenTopography no encontró datos para el área solicitada. Amplíe o modifique el bbox.")
    if resp.status_code == 400:
        raise RuntimeError("400 Bad Request: revise south/north/west/east, demtype y outputFormat. Detalle: " + safe_response_preview_text(resp))
    if resp.status_code == 401:
        raise RuntimeError("401 Unauthorized: API Key incorrecta, vacía o no autorizada.")
    if resp.status_code != 200:
        raise RuntimeError(f"Error HTTP {resp.status_code}: {safe_response_preview_text(resp)}")

    content_type = resp.headers.get("Content-Type", "")
    if not looks_like_geotiff(resp.content, content_type):
        raise RuntimeError("La respuesta fue HTTP 200, pero no parece GeoTIFF válido. Detalle: " + safe_response_preview_text(resp))
    return resp.content


def read_dem_from_bytes(data: bytes) -> Tuple[np.ndarray, Affine, CRS, Dict[str, Any]]:
    with MemoryFile(data) as mem:
        with mem.open() as src:
            arr = src.read(1, masked=True).astype("float64")
            dem = arr.filled(np.nan)
            transform = src.transform
            crs = src.crs
            meta = src.meta.copy()
    if crs is None:
        raise ValueError("El DEM no tiene CRS.")
    return dem, transform, crs, meta


def read_dem_from_upload(uploaded_file) -> Tuple[np.ndarray, Affine, CRS, Dict[str, Any]]:
    return read_dem_from_bytes(uploaded_file.getvalue())


def utm_epsg_from_lonlat(lon: float, lat: float) -> int:
    zone = int(math.floor((lon + 180.0) / 6.0) + 1)
    return (32600 + zone) if lat >= 0 else (32700 + zone)


def reproject_dem_to_local_utm(dem: np.ndarray, transform: Affine, crs: CRS, reference_lon: float, reference_lat: float) -> Tuple[np.ndarray, Affine, CRS]:
    if not crs.is_geographic:
        return dem, transform, crs

    dst_epsg = utm_epsg_from_lonlat(reference_lon, reference_lat)
    dst_crs = CRS.from_epsg(dst_epsg)

    height, width = dem.shape
    with tempfile.NamedTemporaryFile(suffix=".tif") as tmp_src:
        profile = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": "float32",
            "crs": crs,
            "transform": transform,
            "nodata": -9999.0,
        }
        with rasterio.open(tmp_src.name, "w", **profile) as dst:
            dst.write(np.where(np.isfinite(dem), dem, -9999.0).astype("float32"), 1)

        with rasterio.open(tmp_src.name) as src:
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds
            )
            dst_arr = np.full((dst_height, dst_width), -9999.0, dtype="float32")
            reproject(
                source=rasterio.band(src, 1),
                destination=dst_arr,
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=-9999.0,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_nodata=-9999.0,
                resampling=Resampling.bilinear,
            )

    dst_arr = dst_arr.astype("float64")
    dst_arr[dst_arr <= -9990] = np.nan
    return dst_arr, dst_transform, dst_crs


def transform_point(lon: float, lat: float, dst_crs: CRS) -> Tuple[float, float]:
    tr = Transformer.from_crs(PyCRS.from_epsg(4326), PyCRS.from_user_input(dst_crs), always_xy=True)
    return tr.transform(lon, lat)


def transform_line_ll(coords: List[Tuple[float, float]], dst_crs: CRS) -> LineString:
    tr = Transformer.from_crs(PyCRS.from_epsg(4326), PyCRS.from_user_input(dst_crs), always_xy=True)
    pts = [tr.transform(lon, lat) for lon, lat in coords]
    return LineString(pts)


def point_shift_m(lon: float, lat: float, snapped_x: float, snapped_y: float, crs: CRS) -> float:
    x, y = transform_point(lon, lat, crs)
    return float(math.hypot(snapped_x - x, snapped_y - y))


# =============================================================================
# Hydrology D8
# =============================================================================

def rc_from_xy(transform: Affine, x: float, y: float) -> Tuple[int, int]:
    col, row = ~transform * (x, y)
    return int(math.floor(row)), int(math.floor(col))


def xy_from_rc(transform: Affine, row: int, col: int) -> Tuple[float, float]:
    x, y = transform * (col + 0.5, row + 0.5)
    return float(x), float(y)


def priority_flood_fill(dem: np.ndarray, valid: np.ndarray) -> np.ndarray:
    nrows, ncols = dem.shape
    filled = np.array(dem, dtype="float64", copy=True)
    visited = np.zeros_like(valid, dtype=bool)
    heap = []

    def push(r, c):
        if valid[r, c] and not visited[r, c]:
            visited[r, c] = True
            heapq.heappush(heap, (filled[r, c], r, c))

    # Bordes
    for c in range(ncols):
        push(0, c)
        push(nrows - 1, c)
    for r in range(nrows):
        push(r, 0)
        push(r, ncols - 1)

    while heap:
        z, r, c = heapq.heappop(heap)
        for dr, dc in D8_OFFSETS:
            rr, cc = r + dr, c + dc
            if rr < 0 or rr >= nrows or cc < 0 or cc >= ncols:
                continue
            if not valid[rr, cc] or visited[rr, cc]:
                continue
            visited[rr, cc] = True
            if filled[rr, cc] < z:
                filled[rr, cc] = z
            heapq.heappush(heap, (filled[rr, cc], rr, cc))

    return filled


def compute_d8_flow(filled_dem: np.ndarray, valid: np.ndarray, res_x: float, res_y: float) -> np.ndarray:
    nrows, ncols = filled_dem.shape
    flow_to = np.full((nrows, ncols), -1, dtype=np.int64)
    diag = math.hypot(res_x, res_y)
    distances = np.array([diag, res_y, diag, res_x, res_x, diag, res_y, diag], dtype="float64")

    for r in range(nrows):
        for c in range(ncols):
            if not valid[r, c]:
                continue
            z = filled_dem[r, c]
            best_idx = -1
            best_score = 0.0
            current = r * ncols + c
            for k, (dr, dc) in enumerate(D8_OFFSETS):
                rr, cc = r + dr, c + dc
                if rr < 0 or rr >= nrows or cc < 0 or cc >= ncols or not valid[rr, cc]:
                    continue
                dz = z - filled_dem[rr, cc]
                if dz > 0:
                    score = dz / distances[k]
                elif abs(dz) <= 1e-10 and (rr * ncols + cc) < current:
                    score = 1e-12
                else:
                    score = -1.0
                if score > best_score:
                    best_score = score
                    best_idx = rr * ncols + cc
            flow_to[r, c] = best_idx
    return flow_to


def compute_accumulation(flow_to: np.ndarray, valid: np.ndarray) -> np.ndarray:
    nrows, ncols = flow_to.shape
    n = nrows * ncols
    flat_to = flow_to.ravel()
    valid_flat = valid.ravel()

    indeg = np.zeros(n, dtype=np.int32)
    srcs = np.where(valid_flat & (flat_to >= 0))[0]
    tgts = flat_to[srcs]
    np.add.at(indeg, tgts, 1)

    acc = np.zeros(n, dtype="float64")
    acc[valid_flat] = 1.0
    queue = list(np.where(valid_flat & (indeg == 0))[0])
    head = 0
    while head < len(queue):
        i = queue[head]
        head += 1
        j = flat_to[i]
        if j >= 0:
            acc[j] += acc[i]
            indeg[j] -= 1
            if indeg[j] == 0:
                queue.append(int(j))
    return acc.reshape((nrows, ncols))


def snap_to_accumulation(outlet_rc: Tuple[int, int], accumulation: np.ndarray, valid: np.ndarray, radius_cells: int) -> Tuple[int, int]:
    r, c = outlet_rc
    nrows, ncols = accumulation.shape
    r0, r1 = max(0, r - radius_cells), min(nrows, r + radius_cells + 1)
    c0, c1 = max(0, c - radius_cells), min(ncols, c + radius_cells + 1)
    win = accumulation[r0:r1, c0:c1]
    vwin = valid[r0:r1, c0:c1]
    if not np.any(vwin):
        raise ValueError("El punto de control cae fuera del DEM o sobre NoData.")
    masked = np.where(vwin, win, -np.inf)
    local = np.unravel_index(np.nanargmax(masked), masked.shape)
    return int(r0 + local[0]), int(c0 + local[1])


def upstream_mask(flow_to: np.ndarray, valid: np.ndarray, outlet_rc: Tuple[int, int]) -> np.ndarray:
    nrows, ncols = flow_to.shape
    n = nrows * ncols
    outlet_idx = outlet_rc[0] * ncols + outlet_rc[1]
    flat_to = flow_to.ravel()
    valid_idx = np.where(valid.ravel() & (flat_to >= 0))[0]
    targets = flat_to[valid_idx]
    order = np.argsort(targets)
    targets_sorted = targets[order]
    sources_sorted = valid_idx[order]

    mask_flat = np.zeros(n, dtype=bool)
    stack = [int(outlet_idx)]
    mask_flat[outlet_idx] = True

    while stack:
        cur = stack.pop()
        left = np.searchsorted(targets_sorted, cur, side="left")
        right = np.searchsorted(targets_sorted, cur, side="right")
        for src in sources_sorted[left:right]:
            if not mask_flat[src]:
                mask_flat[src] = True
                stack.append(int(src))
    return mask_flat.reshape((nrows, ncols))


def polygon_from_mask(mask: np.ndarray, transform: Affine):
    geoms = []
    for geom, value in shapes(mask.astype("uint8"), mask=mask, transform=transform):
        if value == 1:
            geoms.append(shape(geom))
    if not geoms:
        raise ValueError("No se pudo construir polígono de cuenca.")
    g = unary_union(geoms)
    if not g.is_valid:
        g = g.buffer(0)
    return g


def mask_touches_boundary(mask: np.ndarray) -> List[str]:
    sides = []
    if mask[0, :].any():
        sides.append("norte")
    if mask[-1, :].any():
        sides.append("sur")
    if mask[:, 0].any():
        sides.append("oeste")
    if mask[:, -1].any():
        sides.append("este")
    return sides


def delineate_one(
    name: str,
    lon: float,
    lat: float,
    dem: np.ndarray,
    transform: Affine,
    crs: CRS,
    flow_to: np.ndarray,
    accumulation: np.ndarray,
    valid: np.ndarray,
    snap_radius_m: float,
) -> BasinResult:
    x, y = transform_point(lon, lat, crs)
    rc = rc_from_xy(transform, x, y)
    if not (0 <= rc[0] < dem.shape[0] and 0 <= rc[1] < dem.shape[1]):
        raise ValueError(f"{name}: punto fuera de DEM.")
    if not valid[rc]:
        raise ValueError(f"{name}: punto sobre NoData.")

    res_x = abs(float(transform.a))
    res_y = abs(float(transform.e))
    radius_cells = max(0, int(round(snap_radius_m / max(res_x, res_y))))
    snapped = snap_to_accumulation(rc, accumulation, valid, radius_cells)
    mask = upstream_mask(flow_to, valid, snapped)
    if np.sum(mask) < 3:
        raise ValueError(f"{name}: cuenca menor a 3 celdas.")

    geom = polygon_from_mask(mask, transform)
    area_km2 = float(geom.area / 1e6)
    sx, sy = xy_from_rc(transform, snapped[0], snapped[1])
    to_ll = Transformer.from_crs(PyCRS.from_user_input(crs), PyCRS.from_epsg(4326), always_xy=True)
    slon, slat = to_ll.transform(sx, sy)
    return BasinResult(
        name=name,
        original_rc=rc,
        snapped_rc=snapped,
        mask=mask,
        geom=geom,
        area_km2=area_km2,
        snapped_x=sx,
        snapped_y=sy,
        snapped_lon=slon,
        snapped_lat=slat,
        accumulation_cells=float(accumulation[snapped]),
    )


def process_dem_and_delineate(
    dem: np.ndarray,
    transform: Affine,
    crs: CRS,
    pc_h_lon: float,
    pc_h_lat: float,
    pc_d_lon: float,
    pc_d_lat: float,
    snap_radius_m: float,
    max_cells: int,
) -> Tuple[np.ndarray, Affine, CRS, BasinResult, BasinResult, np.ndarray, np.ndarray]:
    # Reproyección a UTM si viene en lat/lon
    ref_lon = (pc_h_lon + pc_d_lon) / 2.0
    ref_lat = (pc_h_lat + pc_d_lat) / 2.0
    dem, transform, crs = reproject_dem_to_local_utm(dem, transform, crs, ref_lon, ref_lat)

    if dem.size > max_cells:
        raise ValueError(f"DEM demasiado grande: {dem.size:,} celdas. Límite: {max_cells:,}.")

    valid = np.isfinite(dem)
    if np.sum(valid) < 10:
        raise ValueError("DEM sin celdas válidas suficientes.")

    res_x = abs(float(transform.a))
    res_y = abs(float(transform.e))
    filled = priority_flood_fill(dem, valid)
    flow_to = compute_d8_flow(filled, valid, res_x, res_y)
    accumulation = compute_accumulation(flow_to, valid)

    res_h = delineate_one("Cuenca hidrológica", pc_h_lon, pc_h_lat, dem, transform, crs, flow_to, accumulation, valid, snap_radius_m)
    res_d = delineate_one("Cuenca descarga", pc_d_lon, pc_d_lat, dem, transform, crs, flow_to, accumulation, valid, snap_radius_m)
    return dem, transform, crs, res_h, res_d, flow_to, accumulation


# =============================================================================
# Contours
# =============================================================================

def extract_contours(
    dem: np.ndarray,
    transform: Affine,
    crs: CRS,
    clip_geom,
    interval_m: float,
    max_segments: int,
    simplify_m: float,
) -> List[Dict[str, Any]]:
    if measure is None:
        raise RuntimeError("scikit-image no está instalado. No se pueden generar curvas.")
    if interval_m <= 0:
        raise ValueError("Intervalo de curvas inválido.")

    mask = ~geometry_mask([mapping(clip_geom)], transform=transform, invert=False, out_shape=dem.shape)
    arr = np.where(mask & np.isfinite(dem), dem, np.nan)
    if not np.any(np.isfinite(arr)):
        return []

    zmin = float(np.nanmin(arr))
    zmax = float(np.nanmax(arr))
    if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax <= zmin:
        return []

    start = math.floor(zmin / interval_m) * interval_m
    end = math.ceil(zmax / interval_m) * interval_m
    levels = np.arange(start, end + interval_m * 0.1, interval_m)

    fill_value = zmin - 10.0 * interval_m
    arr2 = np.where(np.isfinite(arr), arr, fill_value)
    results = []
    for elev in levels:
        if len(results) >= max_segments:
            break
        if elev <= zmin or elev >= zmax:
            continue
        raw_contours = measure.find_contours(arr2, float(elev))
        for contour in raw_contours:
            if len(results) >= max_segments:
                break
            if contour.shape[0] < 3:
                continue
            coords = []
            for row, col in contour:
                x, y = transform * (float(col), float(row))
                coords.append((x, y))
            try:
                line = LineString(coords)
            except Exception:
                continue
            if line.length < max(1.0, simplify_m):
                continue
            inter = line.intersection(clip_geom)
            lines = []
            if isinstance(inter, LineString):
                lines = [inter]
            elif isinstance(inter, MultiLineString):
                lines = list(inter.geoms)
            else:
                continue
            for g in lines:
                if len(results) >= max_segments:
                    break
                if simplify_m > 0:
                    g = g.simplify(float(simplify_m), preserve_topology=False)
                if g.is_empty or g.length <= 0 or len(g.coords) < 2:
                    continue
                results.append({"elev": float(elev), "geom": g, "length_m": float(g.length)})
    return results


# =============================================================================
# KML export
# =============================================================================

def geom_to_lonlat(geom, src_crs: CRS):
    tr = Transformer.from_crs(PyCRS.from_user_input(src_crs), PyCRS.from_epsg(4326), always_xy=True)
    return shp_transform(tr.transform, geom)


def line_coords_kml(line_ll: LineString) -> str:
    return " ".join([f"{x:.8f},{y:.8f},0" for x, y in line_ll.coords])


def polygon_placemarks(name: str, geom_proj, src_crs: CRS, style: str, desc: str = "") -> str:
    geom_ll = geom_to_lonlat(geom_proj, src_crs)
    polys = list(geom_ll.geoms) if isinstance(geom_ll, MultiPolygon) else [geom_ll]
    out = []
    for i, p in enumerate(polys, start=1):
        if p.is_empty or not hasattr(p, "exterior"):
            continue
        coords = " ".join([f"{x:.8f},{y:.8f},0" for x, y in p.exterior.coords])
        out.append(f"""
        <Placemark>
          <name>{kml_escape(name)} {i}</name>
          <description>{kml_escape(desc)}</description>
          <styleUrl>#{style}</styleUrl>
          <Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates></LinearRing></outerBoundaryIs></Polygon>
        </Placemark>""")
    return "\n".join(out)


def line_placemark(name: str, geom_proj, src_crs: CRS, style: str, desc: str = "") -> str:
    geom_ll = geom_to_lonlat(geom_proj, src_crs)
    if isinstance(geom_ll, MultiLineString):
        parts = list(geom_ll.geoms)
    else:
        parts = [geom_ll]
    out = []
    for i, line in enumerate(parts, start=1):
        if line.is_empty or len(line.coords) < 2:
            continue
        nm = name if len(parts) == 1 else f"{name} {i}"
        out.append(f"""
        <Placemark>
          <name>{kml_escape(nm)}</name>
          <description>{kml_escape(desc)}</description>
          <styleUrl>#{style}</styleUrl>
          <LineString><tessellate>1</tessellate><coordinates>{line_coords_kml(line)}</coordinates></LineString>
        </Placemark>""")
    return "\n".join(out)


def point_placemark(name: str, lon: float, lat: float, style: str, desc: str = "") -> str:
    return f"""
        <Placemark>
          <name>{kml_escape(name)}</name>
          <description>{kml_escape(desc)}</description>
          <styleUrl>#{style}</styleUrl>
          <Point><coordinates>{lon:.8f},{lat:.8f},0</coordinates></Point>
        </Placemark>"""


def styles_kml() -> str:
    return """
    <Style id="eje"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle></Style>
    <Style id="cuenca_desc"><LineStyle><color>ffff0000</color><width>3</width></LineStyle><PolyStyle><color>46ff0000</color></PolyStyle></Style>
    <Style id="cuenca_hidro"><LineStyle><color>ff00ffff</color><width>3</width></LineStyle><PolyStyle><color>4600ffff</color></PolyStyle></Style>
    <Style id="intercuenca"><LineStyle><color>ff00aa00</color><width>3</width></LineStyle><PolyStyle><color>3300aa00</color></PolyStyle></Style>
    <Style id="curva"><LineStyle><color>ff996633</color><width>1</width></LineStyle></Style>
    <Style id="pc_desc"><IconStyle><color>ff0000ff</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href></Icon></IconStyle></Style>
    <Style id="pc_hidro"><IconStyle><color>ff008000</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/pushpin/grn-pushpin.png</href></Icon></IconStyle></Style>
    <Style id="pc_ajustado"><IconStyle><color>ff00ffff</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/pushpin/red-pushpin.png</href></Icon></IconStyle></Style>
    """


def make_consolidated_kmz(
    cuenca_name: str,
    crs: CRS,
    res_h: BasinResult,
    res_d: BasinResult,
    inter_geom,
    axis_proj: LineString,
    pc_h: ControlPoint,
    pc_d: ControlPoint,
    contours: List[Dict[str, Any]],
    include_intercuenca: bool,
) -> bytes:
    safe_name = cuenca_name.strip() or "Cuenca"
    placemarks = []

    placemarks.append(line_placemark(f"Eje Cauce {safe_name}", axis_proj, crs, "eje", "Eje del cauce ingresado por el usuario."))

    placemarks.append(polygon_placemarks(
        f"Cuenca descarga / soporte - {res_d.area_km2:.3f} km²",
        res_d.geom,
        crs,
        "cuenca_desc",
        "Cuenca generada desde PC-DESCARGA."
    ))
    placemarks.append(point_placemark("PC-DESCARGA original", pc_d.lon, pc_d.lat, "pc_desc"))
    placemarks.append(point_placmark if False else point_placemark("PC-DESCARGA ajustado al drenaje", res_d.snapped_lon, res_d.snapped_lat, "pc_ajustado", f"Acumulación: {res_d.accumulation_cells:.0f} celdas"))

    # Curvas dentro de cuenca descarga
    for c in contours:
        placemarks.append(line_placemark(f"Curva {c['elev']:g} m", c["geom"], crs, "curva", "Curva generada desde DEM; uso topográfico/cartográfico referencial."))

    placemarks.append(polygon_placemarks(
        f"Cuenca hidrológica - {res_h.area_km2:.3f} km²",
        res_h.geom,
        crs,
        "cuenca_hidro",
        "Cuenca generada desde PC-HIDRO."
    ))
    placemarks.append(point_placemark("PC-HIDRO original", pc_h.lon, pc_h.lat, "pc_hidro"))
    placemarks.append(point_placemark("PC-HIDRO ajustado al drenaje", res_h.snapped_lon, res_h.snapped_lat, "pc_ajustado", f"Acumulación: {res_h.accumulation_cells:.0f} celdas"))

    if include_intercuenca and inter_geom is not None and not inter_geom.is_empty:
        placemarks.append(polygon_placemarks(
            f"Intercuenca incremental - {inter_geom.area / 1e6:.3f} km²",
            inter_geom,
            crs,
            "intercuenca",
            "Diferencia geométrica: cuenca descarga - cuenca hidrológica."
        ))

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>Cuenca consolidada {kml_escape(safe_name)}</name>
    {styles_kml()}
    {''.join(placemarks)}
  </Document>
</kml>"""
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("doc.kml", kml.encode("utf-8"))
    return bio.getvalue()


def make_excel(summary_df: pd.DataFrame, attempts_df: pd.DataFrame) -> bytes:
    bio = io.BytesIO()
    with pd.ExcelWriter(bio, engine="openpyxl") as writer:
        summary_df.to_excel(writer, index=False, sheet_name="Resumen")
        attempts_df.to_excel(writer, index=False, sheet_name="Intentos_DEM")
    return bio.getvalue()


def make_zip(kmz: bytes, excel: bytes, summary_json: Dict[str, Any], cuenca_name: str) -> bytes:
    safe = re.sub(r"[^A-Za-z0-9_ÁÉÍÓÚáéíóúÑñ-]+", "_", cuenca_name.strip() or "Cuenca")
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"Cuenca_consolidada_{safe}.kmz", kmz)
        zf.writestr(f"Resumen_Cuenca_consolidada_{safe}.xlsx", excel)
        zf.writestr(f"Resumen_integracion_HidroSed_{safe}.json", json.dumps(summary_json, ensure_ascii=False, indent=2).encode("utf-8"))
    return bio.getvalue()


# =============================================================================
# UI
# =============================================================================

st.title("🌊 HidroSed · Generador de Cuenca Consolidada v1.2")
st.caption("Desde DEM COP30 probado + PC-HIDRO + PC-DESCARGA + eje del cauce → Cuenca consolidada nombrada.")

with st.expander("¿Qué hace esta aplicación?", expanded=True):
    st.markdown(
        """
        Esta versión parte desde los **insumos mínimos**:

        1. **Nombre de la cuenca**  
        2. **Punto de control hidrológico PC-HIDRO**  
        3. **Punto de descarga PC-DESCARGA**  
        4. **Eje del cauce**  

        La aplicación descarga o carga un **DEM común** usando el motor de la app DEM COP30 validada, delimita ambas cuencas, controla desplazamiento de puntos, genera curvas de nivel y entrega:

        **Cuenca consolidada _nombre cuenca_.kmz**

        La aplicación intenta evitar dos problemas: cuencas cortadas y puntos ajustados a otra quebrada. Para ello usa el margen tipo app DEM COP30 y detiene el resultado si el punto ajustado se desplaza demasiado.
        """
    )

with st.expander("Sensibilización: buffer y bbox en palabras simples", expanded=False):
    st.markdown(
        """
        **Buffer**: es una franja de seguridad alrededor de algo.  
        Por ejemplo, si al eje del cauce le das un buffer de 20 km para descargar DEM, significa que la aplicación toma el eje y agrega 20 km alrededor para no cortar la cuenca.

        **BBox**: es la caja rectangular que encierra los puntos y el eje.  
        La descarga del DEM se hace usando esa caja. Si la caja queda chica, la cuenca sale incompleta.  
        Por eso esta app puede ampliar automáticamente la caja hasta que la cuenca no toque el borde del DEM.
        """
    )

left, right = st.columns([0.95, 1.05])

with left:
    st.subheader("1. Insumos principales")
    cuenca_name = st.text_input("Nombre de la cuenca", value="Mi_Cuenca")
    pc_h_file = st.file_uploader("PC-HIDRO · punto de control hidrológico", type=["kmz", "kml"], key="pc_h")
    pc_d_file = st.file_uploader("PC-DESCARGA · punto de descarga", type=["kmz", "kml"], key="pc_d")
    axis_file = st.file_uploader("Eje del cauce", type=["kmz", "kml"], key="axis")

    st.subheader("2. DEM")
    dem_mode = st.radio("Fuente DEM", ["Descargar DEM automático", "Cargar DEM GeoTIFF propio"], index=0)
    api_key = ""
    demtype = "COP30"
    dem_upload = None
    if dem_mode == "Descargar DEM automático":
        api_key = st.text_input("API Key OpenTopography", type="password")
        demtype = st.selectbox("Tipo DEM", ["COP30", "NASADEM", "SRTMGL1", "SRTMGL3", "COP90"], index=0)
        st.caption("Motor de descarga basado en la aplicación DEM COP30 validada por el usuario.")
    else:
        dem_upload = st.file_uploader("DEM GeoTIFF común", type=["tif", "tiff"])

with right:
    st.subheader("3. Configuración avanzada")
    with st.expander("Parámetros DEM / bbox", expanded=True):
        bbox_strategy = st.selectbox(
            "Base para descargar DEM",
            ["PC-HIDRO · método DEM COP30 probado", "PC-DESCARGA", "Puntos + eje"],
            index=0,
            help="PC-HIDRO replica el método de la app DEM COP30: un punto central más margen."
        )
        margin_unit = st.radio("Unidad del margen DEM", ["Kilómetros", "Grados decimales"], index=0, horizontal=True)
        if margin_unit == "Kilómetros":
            initial_buffer_km = st.number_input("Margen alrededor del punto/eje (km)", min_value=1.0, max_value=300.0, value=40.0, step=5.0)
        else:
            initial_buffer_km = st.number_input("Margen alrededor del punto/eje (grados)", min_value=0.01, max_value=5.0, value=0.35, step=0.05)
        expand_to_axis = st.checkbox("Asegurar que PC-DESCARGA y eje queden dentro del bbox", value=True)
        expansion_factor = st.number_input("Factor de ampliación si la cuenca toca borde", min_value=1.1, max_value=3.0, value=1.5, step=0.1)
        max_iterations = st.number_input("Máximo de intentos de descarga DEM", min_value=1, max_value=6, value=3, step=1)
        max_bbox_area = st.number_input("Área máxima bbox DEM (km²)", min_value=500.0, max_value=450000.0, value=50000.0, step=1000.0, help="Control de seguridad. Para cuencas grandes use 25.000 a 50.000 km² o más si corresponde.")
        st.caption("Si el bbox calculado supera este valor, la app avisará, pero continuará mientras no supere 450.000 km².")

    with st.expander("Parámetros hidrológicos / ajuste de puntos", expanded=True):
        max_cells = st.number_input("Máximo de celdas DEM a procesar", min_value=100000, max_value=12000000, value=5000000, step=100000)
        snap_radius_m = st.number_input("Radio ajuste punto a drenaje (m)", min_value=0.0, max_value=5000.0, value=300.0, step=50.0)
        max_snap_shift_m = st.number_input("Máximo desplazamiento permitido del punto ajustado (m)", min_value=0.0, max_value=5000.0, value=300.0, step=50.0)
        st.caption("Si el punto ajustado se mueve más que este límite, la app detiene el resultado para evitar saltos a otra quebrada.")

    st.subheader("4. Curvas")
    contour_interval_m = st.number_input("Equidistancia curvas de nivel (m)", min_value=1.0, max_value=100.0, value=20.0, step=1.0)
    contour_simplify_m = st.number_input("Simplificación curvas (m)", min_value=0.0, max_value=50.0, value=5.0, step=1.0)
    max_contours = st.number_input("Máximo de segmentos de curvas", min_value=100, max_value=20000, value=6000, step=100)
    include_intercuenca = st.checkbox("Incluir polígono de intercuenca en KMZ", value=True)

run = st.button("Generar cuenca consolidada", type="primary", use_container_width=True)

if run:
    try:
        if not cuenca_name.strip():
            st.error("Debe indicar nombre de la cuenca.")
            st.stop()
        if pc_h_file is None or pc_d_file is None or axis_file is None:
            st.error("Debe cargar PC-HIDRO, PC-DESCARGA y eje del cauce.")
            st.stop()
        if dem_mode == "Descargar DEM automático" and not api_key.strip():
            st.error("Debe ingresar API Key de OpenTopography.")
            st.stop()
        if dem_mode == "Cargar DEM GeoTIFF propio" and dem_upload is None:
            st.error("Debe cargar DEM GeoTIFF.")
            st.stop()

        with st.status("Generando cuenca consolidada...", expanded=True) as status:
            status.write("Leyendo PC-HIDRO, PC-DESCARGA y eje...")
            pc_h_lon, pc_h_lat = parse_first_point(pc_h_file)
            pc_d_lon, pc_d_lat = parse_first_point(pc_d_file)
            axis_ll = parse_first_line(axis_file)
            all_bbox_points = [(pc_h_lon, pc_h_lat), (pc_d_lon, pc_d_lat)] + axis_ll

            attempts = []
            final = None
            last_error = None

            if dem_mode == "Cargar DEM GeoTIFF propio":
                status.write("Leyendo DEM propio...")
                dem, transform, crs, meta = read_dem_from_upload(dem_upload)
                dem, transform, crs, res_h, res_d, flow_to, accumulation = process_dem_and_delineate(
                    dem, transform, crs, pc_h_lon, pc_h_lat, pc_d_lon, pc_d_lat,
                    float(snap_radius_m), int(max_cells)
                )
                disp_h = point_shift_m(pc_h_lon, pc_h_lat, res_h.snapped_x, res_h.snapped_y, crs)
                disp_d = point_shift_m(pc_d_lon, pc_d_lat, res_d.snapped_x, res_d.snapped_y, crs)
                if float(max_snap_shift_m) > 0 and (disp_h > float(max_snap_shift_m) or disp_d > float(max_snap_shift_m)):
                    raise ValueError(
                        f"Desplazamiento excesivo al ajustar puntos: PC-HIDRO={disp_h:.1f} m, "
                        f"PC-DESCARGA={disp_d:.1f} m. Reduzca el radio de ajuste o revise el punto/eje."
                    )
                sides_h = mask_touches_boundary(res_h.mask)
                sides_d = mask_touches_boundary(res_d.mask)
                attempts.append({
                    "intento": 1,
                    "modo": "DEM cargado",
                    "bbox_area_km2": None,
                    "area_hidrologica_km2": res_h.area_km2,
                    "area_descarga_km2": res_d.area_km2,
                    "borde_hidrologica": ", ".join(sides_h),
                    "borde_descarga": ", ".join(sides_d),
                    "estado": "procesado",
                })
                final = (dem, transform, crs, res_h, res_d)
            else:
                current_buffer = float(initial_buffer_km)
                for i in range(1, int(max_iterations) + 1):
                    # BBox basado en la app DEM COP30 probada: punto central + margen.
                    if bbox_strategy.startswith("PC-HIDRO"):
                        south, north, west, east = bbox_from_center(pc_h_lat, pc_h_lon, current_buffer, margin_unit)
                    elif bbox_strategy.startswith("PC-DESCARGA"):
                        south, north, west, east = bbox_from_center(pc_d_lat, pc_d_lon, current_buffer, margin_unit)
                    else:
                        # Si se selecciona puntos + eje, el margen siempre se interpreta en km.
                        km_margin = current_buffer if margin_unit == "Kilómetros" else current_buffer * 111.32
                        south, north, west, east = bbox_from_lonlat_points(all_bbox_points, km_margin)

                    if expand_to_axis:
                        south, north, west, east = expand_bbox_to_include_points(south, north, west, east, all_bbox_points, padding_km=2.0)

                    area_bbox = bbox_area_km2(south, north, west, east)
                    if area_bbox > 450000.0:
                        raise ValueError(
                            f"El bbox requerido ({area_bbox:,.0f} km²) supera el límite absoluto de seguridad "
                            f"(450.000 km²). Reduzca el margen."
                        )
                    if area_bbox > float(max_bbox_area):
                        status.write(
                            f"Advertencia: el bbox requerido ({area_bbox:,.0f} km²) supera el máximo configurado "
                            f"({max_bbox_area:,.0f} km²), pero se continúa porque está bajo el límite absoluto. "
                            f"Para evitar este aviso, aumente el máximo permitido."
                        )

                    status.write(f"Intento {i}: descargando DEM {demtype} con margen {current_buffer:.2f} {margin_unit}. BBox aprox.: {area_bbox:,.1f} km²")
                    try:
                        raw = download_dem_opentopo(demtype, south, north, west, east, api_key)
                        dem0, transform0, crs0, meta0 = read_dem_from_bytes(raw)
                        status.write("Procesando DEM y delimitando cuencas...")
                        dem, transform, crs, res_h, res_d, flow_to, accumulation = process_dem_and_delineate(
                            dem0, transform0, crs0, pc_h_lon, pc_h_lat, pc_d_lon, pc_d_lat,
                            float(snap_radius_m), int(max_cells)
                        )
                        disp_h = point_shift_m(pc_h_lon, pc_h_lat, res_h.snapped_x, res_h.snapped_y, crs)
                        disp_d = point_shift_m(pc_d_lon, pc_d_lat, res_d.snapped_x, res_d.snapped_y, crs)
                        if float(max_snap_shift_m) > 0 and (disp_h > float(max_snap_shift_m) or disp_d > float(max_snap_shift_m)):
                            raise ValueError(
                                f"Desplazamiento excesivo al ajustar puntos: PC-HIDRO={disp_h:.1f} m, "
                                f"PC-DESCARGA={disp_d:.1f} m. Reduzca el radio de ajuste o revise el punto/eje."
                            )
                        sides_h = mask_touches_boundary(res_h.mask)
                        sides_d = mask_touches_boundary(res_d.mask)
                        attempts.append({
                            "intento": i,
                            "modo": "OpenTopography",
                            "buffer_km": current_buffer,
                            "bbox_area_km2": area_bbox,
                            "area_hidrologica_km2": res_h.area_km2,
                            "area_descarga_km2": res_d.area_km2,
                            "borde_hidrologica": ", ".join(sides_h),
                            "borde_descarga": ", ".join(sides_d),
                            "estado": "OK sin borde" if not sides_h and not sides_d else "Toca borde; ampliar",
                        })
                        final = (dem, transform, crs, res_h, res_d)
                        if not sides_h and not sides_d:
                            break
                        current_buffer *= float(expansion_factor)
                    except Exception as exc:
                        last_error = exc
                        attempts.append({
                            "intento": i,
                            "modo": "OpenTopography",
                            "buffer_km": current_buffer,
                            "bbox_area_km2": area_bbox,
                            "area_hidrologica_km2": None,
                            "area_descarga_km2": None,
                            "borde_hidrologica": "",
                            "borde_descarga": "",
                            "estado": f"Error: {exc}",
                        })
                        current_buffer *= float(expansion_factor)

                if final is None:
                    raise RuntimeError(f"No se pudo generar DEM/cuencas. Último error: {last_error}")

            dem, transform, crs, res_h, res_d = final
            pc_h_x, pc_h_y = transform_point(pc_h_lon, pc_h_lat, crs)
            pc_d_x, pc_d_y = transform_point(pc_d_lon, pc_d_lat, crs)
            pc_h = ControlPoint("PC-HIDRO", pc_h_lon, pc_h_lat, pc_h_x, pc_h_y)
            pc_d = ControlPoint("PC-DESCARGA", pc_d_lon, pc_d_lat, pc_d_x, pc_d_y)
            axis_proj = transform_line_ll(axis_ll, crs)

            status.write("Calculando intercuenca y validaciones...")
            inter_geom = res_d.geom.difference(res_h.geom)
            inter_area_km2 = float(inter_geom.area / 1e6) if not inter_geom.is_empty else 0.0
            simple_diff = res_d.area_km2 - res_h.area_km2
            pct_inside = float(res_d.geom.intersection(res_h.geom).area / res_h.geom.area * 100.0) if res_h.geom.area > 0 else float("nan")

            status.write("Generando curvas de nivel dentro de cuenca de descarga...")
            contours = extract_contours(
                dem, transform, crs, res_d.geom,
                float(contour_interval_m), int(max_contours), float(contour_simplify_m)
            )

            summary_df = pd.DataFrame([
                {"Parámetro": "Nombre cuenca", "Valor": cuenca_name.strip()},
                {"Parámetro": "CRS DEM procesado", "Valor": crs.to_string()},
                {"Parámetro": "Área cuenca hidrológica (km²)", "Valor": res_h.area_km2},
                {"Parámetro": "Área cuenca descarga (km²)", "Valor": res_d.area_km2},
                {"Parámetro": "Área incremental geométrica (km²)", "Valor": inter_area_km2},
                {"Parámetro": "Diferencia simple descarga - hidrológica (km²)", "Valor": simple_diff},
                {"Parámetro": "% cuenca hidrológica dentro de descarga", "Valor": pct_inside},
                {"Parámetro": "Curvas generadas", "Valor": len(contours)},
                {"Parámetro": "Equidistancia curvas (m)", "Valor": float(contour_interval_m)},
                {"Parámetro": "PC-HIDRO original lon", "Valor": pc_h_lon},
                {"Parámetro": "PC-HIDRO original lat", "Valor": pc_h_lat},
                {"Parámetro": "PC-HIDRO ajustado lon", "Valor": res_h.snapped_lon},
                {"Parámetro": "PC-HIDRO ajustado lat", "Valor": res_h.snapped_lat},
                {"Parámetro": "PC-DESCARGA original lon", "Valor": pc_d_lon},
                {"Parámetro": "PC-DESCARGA original lat", "Valor": pc_d_lat},
                {"Parámetro": "PC-DESCARGA ajustado lon", "Valor": res_d.snapped_lon},
                {"Parámetro": "PC-DESCARGA ajustado lat", "Valor": res_d.snapped_lat},
                {"Parámetro": "Desplazamiento PC-HIDRO ajustado (m)", "Valor": point_shift_m(pc_h_lon, pc_h_lat, res_h.snapped_x, res_h.snapped_y, crs)},
                {"Parámetro": "Desplazamiento PC-DESCARGA ajustado (m)", "Valor": point_shift_m(pc_d_lon, pc_d_lat, res_d.snapped_x, res_d.snapped_y, crs)},
                {"Parámetro": "Motor DEM usado", "Valor": "DEM COP30 probado / bbox por punto + margen"},
            ])
            attempts_df = pd.DataFrame(attempts)

            kmz_bytes = make_consolidated_kmz(
                cuenca_name.strip(), crs, res_h, res_d, inter_geom, axis_proj,
                pc_h, pc_d, contours, bool(include_intercuenca)
            )
            excel_bytes = make_excel(summary_df, attempts_df)
            summary_json = {
                "generated_at": datetime.utcnow().isoformat() + "Z",
                "app": APP_TITLE,
                "nombre_cuenca": cuenca_name.strip(),
                "crs": crs.to_string(),
                "area_hidrologica_km2": res_h.area_km2,
                "area_descarga_km2": res_d.area_km2,
                "area_incremental_geometrica_km2": inter_area_km2,
                "diferencia_simple_km2": simple_diff,
                "pct_hidrologica_dentro_descarga": pct_inside,
                "curvas_generadas": len(contours),
                "equidistancia_curvas_m": float(contour_interval_m),
                "pc_hidro_original": {"lon": pc_h_lon, "lat": pc_h_lat},
                "pc_hidro_ajustado": {"lon": res_h.snapped_lon, "lat": res_h.snapped_lat},
                "pc_descarga_original": {"lon": pc_d_lon, "lat": pc_d_lat},
                "pc_descarga_ajustado": {"lon": res_d.snapped_lon, "lat": res_d.snapped_lat},
                "intentos_dem": attempts,
            }
            zip_bytes = make_zip(kmz_bytes, excel_bytes, summary_json, cuenca_name.strip())

            st.session_state["result"] = {
                "summary_df": summary_df,
                "attempts_df": attempts_df,
                "kmz": kmz_bytes,
                "excel": excel_bytes,
                "zip": zip_bytes,
                "res_h": res_h,
                "res_d": res_d,
                "inter_geom": inter_geom,
                "axis_proj": axis_proj,
                "crs": crs,
                "dem": dem,
                "transform": transform,
                "contours_n": len(contours),
            }
            status.update(label="Cuenca consolidada generada", state="complete")

    except Exception as exc:
        st.exception(exc)
        st.stop()


if "result" in st.session_state:
    r = st.session_state["result"]
    res_h = r["res_h"]
    res_d = r["res_d"]
    st.divider()
    st.header("Resultado: cuenca consolidada")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Cuenca hidrológica", f"{res_h.area_km2:.3f} km²")
    m2.metric("Cuenca descarga", f"{res_d.area_km2:.3f} km²")
    m3.metric("Intercuenca", f"{(r['inter_geom'].area/1e6 if not r['inter_geom'].is_empty else 0):.3f} km²")
    m4.metric("Curvas", f"{r['contours_n']}")

    sides_h = mask_touches_boundary(res_h.mask)
    sides_d = mask_touches_boundary(res_d.mask)
    if sides_h or sides_d:
        st.warning("Una o ambas cuencas tocan el borde del DEM. Revise el resultado o aumente el buffer/bbox máximo.")
        if sides_h:
            st.warning(f"Cuenca hidrológica toca borde: {', '.join(sides_h)}")
        if sides_d:
            st.warning(f"Cuenca descarga toca borde: {', '.join(sides_d)}")
    else:
        st.success("Las cuencas no tocan el borde del DEM procesado.")

    st.subheader("Resumen")
    st.dataframe(r["summary_df"], use_container_width=True)

    with st.expander("Intentos de descarga/procesamiento DEM", expanded=False):
        st.dataframe(r["attempts_df"], use_container_width=True)

    st.subheader("Vista de control")
    try:
        dem = r["dem"]
        transform = r["transform"]
        nrows, ncols = dem.shape
        extent = [transform.c, transform.c + transform.a * ncols, transform.f + transform.e * nrows, transform.f]
        fig, ax = plt.subplots(figsize=(10, 8))
        plot_dem = np.array(dem, copy=True)
        if np.any(np.isfinite(plot_dem)):
            plot_dem[~np.isfinite(plot_dem)] = np.nanmedian(plot_dem)
            ax.imshow(plot_dem, extent=extent, origin="upper", alpha=0.45)

        for geom, label, lw in [
            (res_d.geom, "Cuenca descarga", 2.4),
            (res_h.geom, "Cuenca hidrológica", 2.4),
            (r["inter_geom"], "Intercuenca", 1.5),
        ]:
            if geom is None or geom.is_empty:
                continue
            geoms = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
            first = True
            for g in geoms:
                if hasattr(g, "exterior"):
                    x, y = g.exterior.xy
                    ax.plot(x, y, linewidth=lw, label=label if first else None)
                    first = False
        x, y = r["axis_proj"].xy
        ax.plot(x, y, linewidth=3, label="Eje del cauce")
        ax.scatter([res_h.snapped_x], [res_h.snapped_y], marker="o", s=60, label="PC-HIDRO ajustado")
        ax.scatter([res_d.snapped_x], [res_d.snapped_y], marker="s", s=60, label="PC-DESCARGA ajustado")
        ax.set_aspect("equal", adjustable="box")
        ax.legend()
        ax.set_title("Cuenca consolidada generada")
        st.pyplot(fig)
    except Exception as exc:
        st.warning(f"No se pudo graficar: {exc}")

    st.subheader("Descargas")
    safe_name = re.sub(r"[^A-Za-z0-9_ÁÉÍÓÚáéíóúÑñ-]+", "_", str(r["summary_df"].iloc[0]["Valor"]))
    c1, c2, c3 = st.columns(3)
    c1.download_button("ZIP completo", r["zip"], file_name=f"Cuenca_consolidada_{safe_name}_resultados.zip", mime="application/zip", use_container_width=True)
    c2.download_button("KMZ cuenca consolidada", r["kmz"], file_name=f"Cuenca_consolidada_{safe_name}.kmz", mime="application/vnd.google-earth.kmz", use_container_width=True)
    c3.download_button("Excel resumen", r["excel"], file_name=f"Resumen_Cuenca_consolidada_{safe_name}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True)

st.divider()
st.markdown(
    """
    **Nota:** el producto KMZ generado es la base para el módulo siguiente de preparación hidráulica: eje activo, intercuenca, aporte incremental y curvas incrementadas a lo largo del eje.
    """
)
