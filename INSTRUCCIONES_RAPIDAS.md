# Instrucciones rápidas

1. Ejecutar:
   ```bash
   pip install -r requirements.txt
   streamlit run app.py
   ```

2. Cargar:
   - PC-HIDRO
   - PC-DESCARGA
   - Eje del cauce

3. En DEM automático usar:
   - Base: `PC-HIDRO · método DEM COP30 probado`
   - Margen: `40 km`
   - API Key OpenTopography

4. En ajuste hidrológico:
   - Radio ajuste: `300 m`
   - Máximo desplazamiento: `300 m`

5. Descargar:
   - `Cuenca_consolidada_<nombre>.kmz`
   - Excel resumen
   - ZIP completo


## Si aparece aviso de bbox grande

Para casos como 19.600 km², use:
- Área máxima bbox DEM: 50.000 km²
- Margen: 40 km
- Base DEM: PC-HIDRO · método DEM COP30 probado

La app v1.2 no se detiene por superar un máximo bajo, solo detiene si supera 450.000 km².
