# PosRecoCNN Training Pipeline

Esta carpeta contiene el pipeline de entrenamiento para PosRecoCNN con arquitectura modular y mejorada.

## Estructura de archivos

```
training/
├── PosRecoCNN_Training.ipynb             # Notebook principal de entrenamiento
├── config.py                             # Configuración centralizada
├── utils.py                              # Funciones utilitarias
├── inference_example.py                  # Ejemplo de inferencia
└── README.md                             # Este archivo
```


## Uso rápido

### 1. Entrenamiento
```bash
# Abrir el notebook principal
jupyter notebook PosRecoCNN_Training.ipynb

# O verificar configuración desde línea de comandos
python -c "from config import *; print_config_summary()"
```

### 2. Inferencia
```python
from inference_example import run_inference_pipeline

# Ejecutar inferencia en nuevos datos
predictions, truth = run_inference_pipeline(
    data_file="tu_archivo.root",
    model_path="path/to/exported/model"
)
```

### 3. Configuración personalizada
```python
# Modificar parámetros en config.py
FILTER_CONFIG['min_energy_tpc'] = 100.0  # Cambiar corte de energía
TRAINING_CONFIG['epochs'] = 100          # Cambiar número de épocas
```

## Configuración de rutas

Edita las rutas en `config.py`:

```python
DATA_CONFIG = {
    'training_file': '/tu/path/training_data.root',
    'test_file': '/tu/path/test_data.root',
    # ...
}
```

## Características principales

| Aspecto | Implementación |
|---------|----------------|
| Organización | Módulos separados |
| Configuración | Centralizada |
| Imports | Solo necesarios |
| Rutas | Relativas |
| Normalización | Guardada/reutilizable |
| Código duplicado | Eliminado |
| Documentación | Completa |

## Validación del entorno

```python
from config import validate_paths, print_config_summary

# Verificar configuración
print_config_summary()
validate_paths()
```

## Archivos de salida

Después del entrenamiento:
- `modelo_exportado/` - Modelo TensorFlow exportado
- `normalization_factors.json` - Factores para inferencia
- `weights_nuvT.hdf5.keras` - Mejores pesos del modelo

## Solución de problemas comunes

### ❌ Error: "PMT map not found"
**Solución:** Verifica que los archivos CSV estén en `../pmt_maps/`

### ❌ Error: "Normalization file not found"
**Solución:** Asegúrate de haber entrenado el modelo con este pipeline

### ❌ Predicciones erróneas en inferencia
**Solución:** Usa `normalization_factors.json` del entrenamiento

### ❌ Import error de módulos locales
**Solución:** Ejecuta desde la carpeta `training/`

## Rendimiento esperado

Con la configuración por defecto:
- **Tiempo de entrenamiento:** ~10-15 minutos (con GPU)
- **RMSE típico:** 15-25 cm por coordenada
- **Eficiencia de eventos:** ~39% después de filtros

## Próximos pasos (mejoras futuras)

1. **Arquitectura del modelo:** Evaluar modelos alternativos más ligeros
2. **Augmentación de datos:** Implementar técnicas de augmentación
3. **Hiperparámetros:** Optimización automática con Optuna
4. **Monitoreo:** Integración con TensorBoard
5. **Validación cruzada:** K-fold para mejor evaluación

## Contacto y soporte

Para preguntas o problemas:
1. Revisar este README
2. Verificar `config.py` y `utils.py`
3. Consultar la documentación en el notebook
4. Revisar los logs de entrenamiento

---

**Nota:** Este pipeline mantiene la funcionalidad completa con una organización mejorada, mejor mantenibilidad y usabilidad del código.