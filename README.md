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
| `src/Soilmodels/A3DMohrCoulombStandard.f` | Module `ModMohrCoulomb` with `ESM_MC` — replaces stock file. Calls our external UMAT. |
| `src/Soilmodels/UMAT_MohrCoulombStandard.f` | Standalone `SUBROUTINE UMAT` implementing the generalized Hoek-Brown criterion (2002), plus helper `CarSig`. |

## Integration into Anura3D

1. **Replace** `src/ExternalSoilModel.for` in your Anura3D source tree with the version from this repo.
2. **Replace** `src/Soilmodels/A3DMohrCoulombStandard.f` with the version from this repo.
3. **Add** `src/Soilmodels/UMAT_MohrCoulombStandard.f` to the Anura3D Visual Studio project.
4. **Rebuild** the Anura3D solution.

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

## Dependencies

- Anura3D v2024 or v2025 source code
- Intel Fortran Compiler (as used by Anura3D)
- `PrnSig` subroutine from the Anura3D framework (`src/GetPrinStress.FOR`)
