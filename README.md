# pgscat

CLI para consultar la API REST del **PGS Catalog** y descargar archivos asociados (scoring files, harmonized files, etc).

---

## Instalación

Este proyecto usa **uv** para gestión de entornos y paquetes.

### 1. Clonar el repositorio

```bash
git clone https://github.com/tuusuario/pgscat.git
cd pgscat
```

---

### 2. Instalar con uv

Instala la herramienta como CLI aislada:

```bash
uv tool install .
```

Esto expone el comando:

```bash
pgscat
```

Si no está en el PATH:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

---

### 3. Instalación en modo desarrollo (editable)

Si vas a modificar el código:

```bash
uv tool install --editable .
```

---

## Dependencias

* Python >= 3.12
* requests

Las dependencias se resuelven automáticamente con `uv`.

---

## Uso básico

### Consultar un score

```bash
pgscat score PGS000001
```

---

### Exportar a JSON

Guardar a archivo:

```bash
pgscat --json output.json score PGS000001
```

---

### Buscar scores

```bash
pgscat search-scores --trait EFO_0001360
pgscat search-scores --pmid 25855707
```

---

### Descargar scoring files

Descargar archivo harmonized (GRCh38 por default):

```bash
pgscat download PGS003449
```

Especificar build:

```bash
pgscat download PGS003449 --build GRCh37
```

---

## Licencia

Este código está publicado con la licencia UNLICENSE.

Consulta los términos del **PGS Catalog**:
https://www.ebi.ac.uk/about/terms-of-use/

---

## Notas

* La API puede aplicar rate limiting
* Algunos endpoints son paginados
* Los formatos de scoring files pueden variar

---
