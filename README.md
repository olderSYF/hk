# Hoek-Brown Constitutive Model for Anura3D

This repository contains the **Hoek-Brown (HB) constitutive model** implementation for the
[Anura3D MPM framework](https://github.com/Anura3D/Anura3D_OpenSource).  The code is
structured to match the existing Mohr-Coulomb model exactly, enabling seamless integration
into the Anura3D source tree.

---

## Table of Contents

1. [Repository Contents](#repository-contents)
2. [Hoek-Brown Model Theory](#hoek-brown-model-theory)
3. [Material Parameters](#material-parameters)
4. [State Variables](#state-variables)
5. [Integration into Anura3D — Step-by-Step](#integration-into-anura3d--step-by-step)
6. [GOM File Example](#gom-file-example)
7. [GiD `.bas` Template Additions](#gid-bas-template-additions)
8. [Visual Studio Project Integration](#visual-studio-project-integration)

---

## Repository Contents

| File | Purpose |
|------|---------|
| `src/Soilmodels/A3DHoekBrown.f` | **Main model file** — Fortran module `ModHoekBrown` containing `ESM_HB` (ESM wrapper) and `UMAT_HB` (core algorithm), plus helper subroutines `CarSig`, `MatMat`, `MatTranspose`. |
| `src/Soilmodels/UMAT_HoekBrown.f` | **Standalone DLL version** — ABAQUS-style `UMAT` entry point for use as an external DLL (includes its own `PrnSig`). |
| `src/Soilmodels/aba_param.inc` | Standard ABAQUS parameter include file required by the standalone DLL. |
| `src/ExternalSoilModel_HB_patch.for` | Reference patch showing the two edits needed in `ExternalSoilModel.for`. |
| `src/GlobalConstants_HB_patch.FOR` | Reference patch showing the three constant additions needed in `GlobalConstants.FOR`. |
| `src/ReadMaterialData_HB_patch.FOR` | Reference patch showing the material-reading block needed in `ReadMaterialData.FOR`. |

---

## Hoek-Brown Model Theory

### Generalised Hoek-Brown criterion

The Generalised Hoek-Brown (GHB) failure criterion (Hoek, Carranza-Torres & Corkum, 2002)
for rock masses is expressed in terms of **major** (σ₁) and **minor** (σ₃) principal stresses
in **compression-positive** convention:

```
σ₁ = σ₃ + σci · (mb · σ₃/σci + s)^a
```

where the rock-mass constants are derived from the intact rock properties and the
Geological Strength Index (GSI):

```
mb  = mi · exp[ (GSI − 100) / (28 − 14D) ]
s   = exp[ (GSI − 100) / (9 − 3D) ]
a   = 0.5 + (1/6) · [ exp(−GSI/15) − exp(−20/3) ]
```

- **σci** — intact rock uniaxial compressive strength (UCS), MPa  
- **mi** — intact rock Hoek-Brown constant  
- **GSI** — Geological Strength Index (10–100)  
- **D** — disturbance factor (0 = undisturbed, 1 = fully disturbed blast damage)

### Tensile strength cutoff

The theoretical tensile strength of the rock mass is:

```
σt = s · σci / mb
```

If the user sets `TENS = 0` the model auto-calculates σt from the above formula.

### Equivalent instantaneous Mohr-Coulomb parameters

To use the efficient Anura3D MC return-mapping algorithm, the model linearises the
curved HB envelope at the current confining stress σ₃ to obtain instantaneous
MC-equivalent parameters:

```
dσ₁/dσ₃ = 1 + a·mb·(mb·σ₃/σci + s)^(a−1)

sin φ_eq = (dσ₁/dσ₃ − 1) / (dσ₁/dσ₃ + 1)

c_eq = [σ₁(σ₃) − (dσ₁/dσ₃) · σ₃] · (1 − sin φ_eq) / 2
```

These equivalent parameters are then fed into the standard MC nine-zone
return-mapping algorithm (Areas 1–9) together with the user-specified dilation
angle sin(ψ) and tensile strength cutoff.

---

## Material Parameters

These values are passed to the model via `PROPS(1:8)`:

| Index | Name | Description | Units | Typical range |
|-------|------|-------------|-------|---------------|
| 1 | G | Shear modulus | MPa | 1 000 – 30 000 |
| 2 | ENU | Poisson's ratio | — | 0.15 – 0.35 |
| 3 | SIGCI | Intact rock UCS | MPa | 5 – 300 |
| 4 | XMI | Intact rock constant mi | — | 5 – 33 |
| 5 | GSI | Geological Strength Index | — | 10 – 100 |
| 6 | DISTD | Disturbance factor D | — | 0 – 1 |
| 7 | SPSI | sin(dilation angle ψ) | — | 0 – 1 |
| 8 | TENS | Tensile strength cutoff | MPa | 0 (auto) or > 0 |

> **Note:** Anura3D uses **tension-positive** sign convention internally, while the HB
> criterion is written in **compression-positive** form.  The model converts between
> the two conventions automatically.

---

## State Variables

These values are stored in `STATEV(1:4)`:

| Index | Name | Description | Units |
|-------|------|-------------|-------|
| 1 | EPSP | Accumulated equivalent plastic strain | — |
| 2 | IPL | Yield flag: 0 = elastic, 1 = shear yield, 2 = tension cutoff | — |
| 3 | PHI_EQ | Current instantaneous equivalent friction angle | degrees |
| 4 | COH_EQ | Current instantaneous equivalent cohesion | MPa |

---

## Integration into Anura3D — Step-by-Step

### Prerequisites

- A working copy of [Anura3D_OpenSource](https://github.com/Anura3D/Anura3D_OpenSource)
- Intel Fortran Compiler (ifort) or compatible
- Visual Studio (Windows) or make-based build system

### Step 1 — Copy source files

Copy the two source files into the Anura3D soil-model directory:

```
cp src/Soilmodels/A3DHoekBrown.f       <Anura3D>/src/Soilmodels/
cp src/Soilmodels/UMAT_HoekBrown.f     <Anura3D>/src/Soilmodels/
cp src/Soilmodels/aba_param.inc        <Anura3D>/src/Soilmodels/
```

### Step 2 — Update `GlobalConstants.FOR`

Open `<Anura3D>/src/GlobalConstants.FOR` and add the three constants shown in
`src/GlobalConstants_HB_patch.FOR`:

```fortran
! After MOHR_COULOMB_TEUNISSEN = 109
integer(INTEGER_TYPE), parameter :: HOEK_BROWN = 110

! After ESM_MOHR_COULOMB_STRAIN_SOFTENING
character(len=64), parameter :: ESM_HOEK_BROWN = 'hoek_brown'

! After UMAT_MOHR_COULOMB_STRAIN_SOFTENING
character(len=64), parameter :: UMAT_HOEK_BROWN = 'A3DHoekBrown.dll'
```

### Step 3 — Update `ExternalSoilModel.for`

Open `<Anura3D>/src/ExternalSoilModel.for`.

**3a.** Add the `use` statement (near the other `use Mod*` statements):
```fortran
use ModHoekBrown
```

**3b.** Add the `elseif` branch inside `StressSolid` (after the MC branch):
```fortran
elseif (trim(NameModel)//char(0) == trim(ESM_HOEK_BROWN)//char(0)) then
  props(1) = Particles(IDpt)%ShearModulus
  props(2) = MatParams(IDSet)%PoissonRatio
  props(3) = MatParams(IDSet)%ESM_Solid(3)   ! SIGCI
  props(4) = MatParams(IDSet)%ESM_Solid(4)   ! mi
  props(5) = MatParams(IDSet)%ESM_Solid(5)   ! GSI
  props(6) = MatParams(IDSet)%ESM_Solid(6)   ! D
  props(7) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
  props(8) = MatParams(IDSet)%TensileStrength
  cmname = UMAT_HOEK_BROWN
```

See `src/ExternalSoilModel_HB_patch.for` for full context.

### Step 4 — Update `ReadMaterialData.FOR`

Follow the instructions in `src/ReadMaterialData_HB_patch.FOR` to add
reading of SIGCI, mi, GSI, D, dilation angle, and tensile strength from the
GOM file into `ESM_Solid(3:8)`.

### Step 5 — Add to the ESM pointer registration

In the routine that populates `MatParams(IDSet)%ESM_POINTER` (typically in
`ReadMaterialData.FOR` or `ModMPMInit`), add:

```fortran
case (HOEK_BROWN)
  MatParams(IDSet)%ESM_POINTER => ESM_HB
```

### Step 6 — Add source file to the build

**Visual Studio:** Right-click the `Soilmodels` filter → Add → Existing item →
select `A3DHoekBrown.f`.

**make/CMake:** Add `src/Soilmodels/A3DHoekBrown.f` to the list of compiled
Fortran sources.

### Step 7 — Build and test

Recompile Anura3D.  Run a simple uniaxial compression benchmark to verify the
model activates and produces physically reasonable principal stresses.

---

## GOM File Example

```
$$MATERIAL_INDEX         1
$$MATERIAL_MODEL         hoek_brown
$$DENSITY_SOLID          2700.0
$$YOUNG_MODULUS          20000.0
$$POISSON_RATIO          0.25
$$SHEAR_MODULUS          8000.0
$$SIGCI                  80.0
$$MI                     10.0
$$GSI                    65.0
$$DISTURBANCE            0.0
$$DILATANCY_ANGLE        5.0
$$TENSILE_STRENGTH       0.0
```

---

## GiD `.bas` Template Additions

In the Anura3D GiD problem-type `.bas` template, add the following block in the
material data section so that GiD exports the Hoek-Brown parameters:

```
*if(strcmp(MatModel,"hoek_brown")==0)
$$MATERIAL_MODEL hoek_brown
$$SHEAR_MODULUS *ShearModulus*
$$POISSON_RATIO *PoissonRatio*
$$SIGCI *SIGCI*
$$MI *MI*
$$GSI *GSI*
$$DISTURBANCE *DisturbanceFactor*
$$DILATANCY_ANGLE *DilatancyAngle*
$$TENSILE_STRENGTH *TensileStrength*
*endif
```

Add the following to the material question definitions (`.cnd` or `.mat` file):

```
QUESTION: SIGCI: 80.0
DEPENDENCIES: (CONDITION, Hoek_Brown, ACTIVATE, SIGCI)
HELP: Intact rock uniaxial compressive strength (MPa)

QUESTION: MI: 10.0
DEPENDENCIES: (CONDITION, Hoek_Brown, ACTIVATE, MI)
HELP: Intact rock Hoek-Brown constant mi

QUESTION: GSI: 65.0
DEPENDENCIES: (CONDITION, Hoek_Brown, ACTIVATE, GSI)
HELP: Geological Strength Index (10-100)

QUESTION: DisturbanceFactor: 0.0
DEPENDENCIES: (CONDITION, Hoek_Brown, ACTIVATE, DisturbanceFactor)
HELP: Disturbance factor D (0=undisturbed, 1=fully disturbed)
```

---

## Visual Studio Project Integration

1. In the Solution Explorer, expand the `Soilmodels` filter under the Anura3D project.
2. Right-click → **Add** → **Existing Item** → select `A3DHoekBrown.f`.
3. Right-click `A3DHoekBrown.f` → **Properties** → ensure **Item Type** is
   *Fortran Compilation*.
4. For the **standalone DLL** (`UMAT_HoekBrown.f`), create a separate DLL project
   in the same solution, add `UMAT_HoekBrown.f` and `aba_param.inc`, and set the
   output name to `A3DHoekBrown.dll`.
5. Rebuild the solution.

---

## References

- Hoek, E., Carranza-Torres, C., & Corkum, B. (2002). *Hoek-Brown failure criterion —
  2002 edition.* NARMS-TAC Conference, Toronto, 1, 267-273.
- Hoek, E., & Brown, E.T. (1980). *Underground excavations in rock.* IMM, London.
- Anura3D MPM Research Community (2025). *Anura3D OpenSource — Theory and Validation Manual.*

