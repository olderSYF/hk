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

### `LNK1120: 1 个无法解析的外部命令` / `LNK2019: unresolved external symbol FEXIST`

> `LNK1120` is the **summary** error — it simply counts how many unresolved
> symbols remain. Fix the underlying `LNK2019` error(s) below and `LNK1120`
> disappears automatically.

This error is **not caused by this repository's code**. `FEXIST` is a non-standard
Intel Fortran extension (file-existence check) used elsewhere in the Anura3D
codebase. It is provided by the Intel Fortran **Portability library** (`libifport`).

**Fix — recommended:** Always build from an **Intel oneAPI / Intel Fortran
command prompt** (or run `setvars.bat` first). This sets all library search
paths automatically, and the linker will find `libifport.lib` without any
manual configuration:

```bat
:: Open "Intel oneAPI command prompt for Intel 64" from the Start Menu, or:
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
```

**Fix — manual (Visual Studio):**

1. Right-click the Anura3D project → **Properties** → **Fortran** →
   **Libraries** → **Runtime Library** — set to
   **Debug Multithread DLL** (`/libs:dll /dbglibs`) for Debug, or
   **Multithread DLL** (`/libs:dll`) for Release. This tells Intel Fortran
   to link its portability library automatically.
2. If the above is not sufficient, explicitly add the Intel Fortran `lib`
   directory to **Linker → General → Additional Library Directories**.
   Typical paths:
   - oneAPI: `C:\Program Files (x86)\Intel\oneAPI\compiler\latest\lib`
   - Classic: `C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\lib\intel64`

### `LNK1181: cannot open input file "libifportd.lib"`

This means the Intel Fortran **Portability library** is referenced but the
linker cannot find it on the search path. This often happens when building
from a plain Visual Studio Developer Command Prompt that does not include
Intel Fortran paths.

**Fix (any of these):**

1. **Use the Intel oneAPI command prompt** (or run `setvars.bat intel64`
   before building). This is the simplest and most reliable fix.
2. Add the Intel Fortran `lib` directory to **Linker → General →
   Additional Library Directories** (see paths above).
3. Copy `libifportd.lib` / `libifport.lib` from your Intel compiler's
   `lib` directory into the Anura3D project's library folder — this is
   a last-resort workaround.

### `LNK4272: library machine type 'x86' conflicts with target machine type 'x64'`

This means one or more linked libraries were built for 32-bit (x86) while the
Anura3D project targets 64-bit (x64). The most common culprit is **HDF5**.

**Fix:**

1. Check which `.lib` files trigger the warning (the full linker output names
   the offending library, e.g. `hdf5.lib`).
2. Replace the x86 `.lib` (and `.dll`) with the **x64 build** of the same
   library. For HDF5, download the 64-bit binaries from
   [The HDF Group](https://www.hdfgroup.org/downloads/hdf5/) and copy the
   `lib/` and `bin/` contents into Anura3D's library directory.
3. Verify the project platform is set to **x64** in Visual Studio
   (**Build → Configuration Manager → Active solution platform → x64**).
4. If you also see `LNK2019: FEXIST` at the same time, fix **both** issues:
   the x86/x64 mismatch (this section) **and** the missing portability
   library (see FEXIST section above).

### `LNK4099: PDB 'hdf5.pdb' not found`

This is a harmless warning — the HDF5 library was built without debug symbols.
It can be safely ignored or suppressed via **Linker → Command Line →
Additional Options**: `/ignore:4099`.

## Dependencies

- Anura3D v2024 or v2025 source code
- Intel Fortran Compiler (as used by Anura3D)
- `PrnSig` subroutine from the Anura3D framework (`src/GetPrinStress.FOR`)
