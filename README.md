# Hoek-Brown Constitutive Model for Anura3D

This repository contains a Hoek-Brown (HB) constitutive model implementation
for the [Anura3D](https://github.com/Anura3D/Anura3D_OpenSource) MPM framework
(v2024 / v2025).

## Approach: 借壳法 (Shell Borrowing)

The Hoek-Brown model is integrated via the Mohr-Coulomb (MC) channel in
Anura3D's External Soil Model (ESM) framework. The MC props slot is reused
to pass raw HB parameters directly to a custom UMAT, avoiding the standard
`SIN()` transforms applied to friction/dilatancy angles.

### PROPS mapping (Hoek-Brown via MC channel)

| Slot | Original MC meaning    | HB meaning (借壳法) | Source field |
|------|------------------------|----------------------|--------------|
| 1    | G (shear modulus)      | G                    | `Particles%ShearModulus` |
| 2    | ν (Poisson's ratio)    | ν                    | `MatParams%PoissonRatio` |
| 3    | sin(φ)                 | σ_ci (UCS)           | `MatParams%FrictionAngle` (raw, no SIN) |
| 4    | c·cos(φ)               | GSI                  | `MatParams%DilatancyAngle` (raw, no SIN) |
| 5    | sin(ψ)                 | m_i                  | `MatParams%TensileStrength` (direct) |
| 6    | σ_t                    | D (disturbance)      | `MatParams%ESM_Solid(6)` |

## Files

| File | Purpose |
|------|---------|
| `src/ExternalSoilModel.for` | Module `ModExternalSoilModel` — replaces stock Anura3D file. Only the MC branch props assignment is modified. |
| `src/Soilmodels/A3DMohrCoulombStandard.f` | Module `ModMohrCoulomb` with `ESM_MC`, `UMAT` (Hoek-Brown), and `CarSig` — replaces stock file. |

## Integration into Anura3D

1. **Replace** `src/ExternalSoilModel.for` in your Anura3D source tree with the version from this repo.
2. **Replace** `src/Soilmodels/A3DMohrCoulombStandard.f` with the version from this repo.
3. **Clean build**: delete all `.obj` and `.mod` files in the output directory (stale artifacts cause `LNK2019` errors).
4. **Rebuild** the Anura3D solution.

> **Note:** The UMAT and CarSig subroutines are now inside the `ModMohrCoulomb` module (no separate file needed). This eliminates `LNK2019: unresolved external symbol UMAT` linker errors.

### GOM-file setup

In the Anura3D GOM input file, use the MC material model keyword and assign
the HB parameters to the corresponding MC fields:

- `FrictionAngle` → σ_ci (MPa)
- `DilatancyAngle` → GSI (0–100)
- `TensileStrength` → m_i
- `ESM_Solid(6)` → D (0–1)

## Hoek-Brown Criterion

The generalized Hoek-Brown criterion (Hoek, Carranza-Torres & Corkum, 2002):

```
F = σ₃ - σ₁ - σ_ci · (m_b · σ₁/σ_ci + s)^a
```

Derived parameters:
- `m_b = m_i · exp((GSI - 100) / (28 - 14D))`
- `s = exp((GSI - 100) / (9 - 3D))`
- `a = 0.5 + (1/6) · (exp(-GSI/15) - exp(-20/3))`

Sign convention: compression positive (consistent with Anura3D).

## Troubleshooting

### `LNK2019: unresolved external symbol FEXIST`

This error is **not caused by this repository's code**. `FEXIST` is a non-standard
Intel Fortran extension (file-existence check) used elsewhere in the Anura3D
codebase. It is provided by the Intel Fortran **Portability library**.

**Fix:** In Visual Studio, add the portability library to the linker inputs:

1. Right-click the Anura3D project → **Properties** → **Linker** → **Input**
2. Add `libifport.lib` (Release) or `libifportd.lib` (Debug) to
   **Additional Dependencies**

Alternatively, ensure your Intel Fortran installation's `lib` directory
(e.g. `C:\Program Files (x86)\Intel\oneAPI\compiler\latest\lib`) is in
**Linker → General → Additional Library Directories**.

### `LNK4272: library machine type 'x86' conflicts with target machine type 'x64'`

This means one or more linked libraries (commonly HDF5) were built for 32-bit
while the Anura3D project targets 64-bit. Replace the x86 `.lib` files with
their x64 equivalents, or switch the project platform to match.

### `LNK4099: PDB 'hdf5.pdb' not found`

This is a harmless warning — the HDF5 library was built without debug symbols.
It can be safely ignored or suppressed via **Linker → Command Line →
Additional Options**: `/ignore:4099`.

## Dependencies

- Anura3D v2024 or v2025 source code
- Intel Fortran Compiler (as used by Anura3D)
- `PrnSig` subroutine from the Anura3D framework (`src/GetPrinStress.FOR`)
