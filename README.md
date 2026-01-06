# GA Elevation → LandXML (TIN) Export

This code converts elevation data retrieved from the **Geoscience Australia (GA) elevation service** into a **LandXML TIN surface**.

The output is a standards-compatible LandXML file containing a triangulated surface definition that can be consumed by civil and GIS tools that support LandXML (e.g. Civil 3D, Bentley, ESRI workflows).

---

## Overview

The pipeline performs the following steps:

```
GA Elevation Data (DEM / grid samples)
            ↓
DEM → TIN triangulation
            ↓
LandXML <Pnts> + <Faces> export
```

The intent is **format conversion only**:
- no tiling
- no PostGIS
- no surface storage
- no clipping or solids

This makes the output suitable for interchange and downstream processing.

---

## What the code does

### 1. Reads elevation data from GA

The GA service typically returns **DEM-style elevation samples**:
- regularly spaced X/Y coordinates
- a Z value per sample
- projected coordinates (SRID provided by the caller)

The importer reads this data in a streaming fashion so large areas can be handled efficiently.

---

### 2. Converts DEM samples into a TIN

Each DEM grid cell is converted into triangles.

For a regular grid cell:

```
v00 ---- v10
 |  \     |
 |   \    |
 |    \   |
v01 ---- v11
```

By default:
- **2 triangles per grid cell**
- vertices are reused across faces
- Z values come directly from the DEM samples

This produces a **piecewise planar surface** equivalent to the DEM at its native resolution.

> If the GA source were to provide irregular samples instead of a grid, this stage could be extended to use a Delaunay triangulation library. For GA DEM products, this is typically not required.

---

### 3. Writes a LandXML TIN surface

The exporter writes a LandXML document using the widely supported TIN structure:

```xml
<LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2"
         version="1.2"
         units="metric">
  <Surfaces>
    <Surface name="GA Terrain">
      <Definition surfType="TIN">
        <Pnts>
          <P id="1">x y z</P>
          <P id="2">x y z</P>
          ...
        </Pnts>
        <Faces>
          <F>1 2 3</F>
          <F>2 4 3</F>
          ...
        </Faces>
      </Definition>
    </Surface>
  </Surfaces>
</LandXML>
```

Characteristics:
- Point IDs are **1-based**
- Faces reference point IDs only (no duplication of coordinates)
- Output matches the de-facto LandXML TIN format used by ESRI, Autodesk, and Bentley tools

---

## Inputs

- Area of interest (bounding box or equivalent)
- SRID / projected coordinate system (explicitly supplied)
- DEM resolution / sampling (as supported by the GA service)
- Optional handling rules for no-data values

---

## Output

- A single **LandXML `.xml` file**
- Contains one TIN surface
- Compatible with LandXML 1.x readers

---

## Coordinate system handling

LandXML does not reliably encode EPSG/SRID information in a consistent, machine-readable way.

For this reason:
- The exporter assumes the coordinates are already in the correct projected CRS
- SRID is treated as an **external contract**
- Consumers are expected to know the CRS of the file

---

## No-data handling

If the GA elevation data contains no-data cells:
- triangles that would include no-data vertices can be skipped
- this results in holes in the TIN
- alternative fill or extrapolation strategies can be added if required

The default behaviour prioritises correctness over surface completeness.

---

## Triangle winding

Triangles are generated with a consistent vertex order (typically counter-clockwise).

Most LandXML consumers do not require a specific winding order, but consistency helps with downstream validation and processing.

---

## Typical use cases

- Exporting GA terrain into civil design tools
- Creating LandXML test surfaces for QA or development
- Interchanging DEM-derived terrain with systems that expect TIN input
- Preparing terrain data for later tiling or surface processing pipelines

---

## Summary

This code is a focused utility that:

1. Reads elevation data from the GA service  
2. Converts DEM samples into a triangulated surface (TIN)  
3. Exports the result as a **LandXML TIN** using the standard `<Pnts>` and `<Faces>` structure  

It provides a clean, interoperable bridge between GA elevation products and LandXML-based workflows.
