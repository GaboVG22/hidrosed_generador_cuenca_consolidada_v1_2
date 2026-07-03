# HidroSed · Generador de Cuenca Consolidada v1.2

Versión corregida basada en la aplicación **dem_cop30_streamlit** que descarga DEM COP30 desde OpenTopography y que fue validada por el usuario.

## Insumos

- Nombre de la cuenca
- PC-HIDRO
- PC-DESCARGA
- Eje del cauce
- DEM automático OpenTopography o DEM GeoTIFF propio

## Correcciones v1.2

- Motor de descarga DEM basado en la app DEM COP30 probada.
- BBox por punto central + margen, igual a la lógica de la app DEM COP30.
- Opción para asegurar que PC-DESCARGA y eje queden dentro del bbox.
- Radio de ajuste al drenaje vuelve a 300 m por defecto.
- Control de desplazamiento máximo del punto ajustado.
- Si el punto ajustado salta demasiado, la app detiene el resultado para evitar cuencas erróneas.
- Exporta KMZ consolidado, Excel y JSON.

## Ejecución

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Recomendación práctica

Para replicar el comportamiento de la app DEM COP30 validada:

- Base para descargar DEM: `PC-HIDRO · método DEM COP30 probado`
- Margen: `Kilómetros`
- Valor: `40 km`
- Radio ajuste a drenaje: `300 m`
- Desplazamiento máximo permitido: `300 m`

Si el resultado se detiene por desplazamiento excesivo, no aumente el radio a ciegas: revise el punto, el eje o reduzca el radio.


## Corrección v1.2

- El área máxima bbox DEM por defecto sube a 50.000 km².
- Si el bbox requerido supera el máximo configurado, la aplicación ya no se detiene automáticamente; muestra advertencia y continúa mientras esté bajo el límite absoluto de 450.000 km².
- Esto corrige casos donde el eje + margen genera un bbox de 19.600 km², que es razonable para descargar DEM en cuencas grandes.
