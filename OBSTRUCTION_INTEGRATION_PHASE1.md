# Phase 1: Obstruction Integration Analysis

**Date**: April 27, 2026  
**Branch**: `feature/obstruction-integration`  
**Status**: Analysis Complete

---

## Executive Summary

This document catalogs all API changes and new functionality in `openstream-obs` that need to be integrated into the main `openstream` codebase. The integration will add obstruction modeling capabilities while preserving all existing solver functionality.

⚠️ **CRITICAL RULE**: When features or enumerations exist in BOTH versions, always use the openstream version. The openstream-obs codebase has bugs that were fixed in openstream. The goal is to add ONLY the obstruction-specific functionality to create an improved openstream-obs branch with both bug fixes and obstruction capabilities.

---

## 1. NEW INPUT CLASSES

### 1.1 `Inputs/Obstruction.m`
- **Purpose**: Defines a reusable wall obstruction with geometric properties
- **Key Properties**:
  - `ID`: Obstruction identifier (string)
  - `SHAPE`: Geometric shape enum (CIRCLE, TRIANGLE, DIAMOND, SQUARE, CUSTOM)
  - `WALL`: Wall index reference
  - `RLOCATION`: Relative location on wall [axial, spanwise] (2×1 double)
  - Shape-specific geometry:
    - **Circle**: `DIAMETER`
    - **Triangle**: `SIDELENGTH` (3×1), `POINTANGLE`
  - Calculated properties: `LENGTH`, `WIDTH`, `AREA`, `PERIM`, `NWALL`
- **Methods**:
  - Constructor with file path and obstruction ID
  - `placeObsOnWall()`: Determines axial and span bounds
- **Input File Format**: Defined in `obstruction.inp`
- **Status**: NEW - Does not exist in openstream

---

## 2. NEW INPUT ENUMERATIONS

### 2.1 `InputEnums/OBSSHAPE.m`
- **Purpose**: Enumeration for obstruction geometric shapes
- **Values**: `CIRCLE`, `TRIANGLE`, `DIAMOND`, `SQUARE`, `CUSTOM`
- **Status**: NEW

### 2.2 `InputEnums/SCBOIL.m`
- **Purpose**: Subcooled boiling model selection
- **Values**: `NONE` (currently only one option)
- **Note**: Framework for future subcooled boiling models
- **Status**: NEW

### 2.3 `InputEnums/SOLVER.m`
- **Purpose**: Enumeration of available solver types
- **Values**: `MIXTURE`, `TWOFLUID`, `THREEFIELD`, `FOURFIELD`, `OBSTRUCTION`
- **New Value**: `OBSTRUCTION` (adds to existing 4 solver types)
- **Key Method**: `solverPath()` - Returns the full class path for each solver type
- **Status**: MODIFIED (exists in openstream-obs, may need to check if openstream has this enum)

---

## 3. NEW SOLVER CLASSES

### 3.1 `Solvers/SolverMode.m`
- **Purpose**: Enumeration for solver execution modes
- **Values**:
  - `NEW` (0): New solution
  - `CONTINUE` (1): Continue existing solution
  - `SUBSET` (2): Solve subset of domain
- **Inheritance**: uint16
- **Status**: NEW - Adds capability to track solver mode vs just state

### 3.2 `Solvers/Obstruction/SolutionSet.m`
- **Purpose**: Organizes solution sets in obstruction solver (handles OBS and NONOBS regions)
- **Key Properties**:
  - `solutionType`: 'OBS' or 'NONOBS'
  - `inputSet`: Associated InputSet
  - `axialStartPosition`, `axialEndPosition`: Domain bounds
  - `solver`: Reference to parent solver
- **Methods**:
  - Constructor
  - `setSolver()`: Associates solver with solution set
- **Status**: NEW

### 3.3 `Solvers/Obstruction/@ObstructionSolver/ObstructionSolver.m`
- **Purpose**: Main solver class for obstructed flow domains
- **Inheritance**: `Solvers.AbstractSolver`
- **Key Properties**:
  - Grid: `NZ`, `Z`, `DZ` (axial)
  - Time: `NTIME`, `TIME`, `DT`
  - Domain segmentation: `solutionSets`, `axialBounds`, `spanBounds`
  - `solverType`: InputEnums.SOLVER enum value
  - `SOLVERMODE`: SolverMode enum
- **Constructor Arguments**:
  - `inputSet`: InputSet with obstruction(s)
  - `solverType`: SOLVER enum value
- **Key Methods**:
  - `initializeSolver()`: Sets up solution structure
  - `solve()`: Executes solution (segmented by obstruction regions)
- **Architecture**:
  - Splits inputSet into axial segments (before, at, after obstruction)
  - Further subdivides obstruction region by spanwise position
  - Manages multiple solution sets for coupled solution
- **Status**: NEW

---

## 4. MODIFICATIONS TO EXISTING CLASSES

### 4.1 `Inputs/Model.m`
**New Property for Obstruction Solver**:
```matlab
OBSWSPLITRATIO  (2,1) double {mustBeNonnegative}
    = [0.0266 0.001] .* [0.5 0.5]  
    % Split ratio of mass downstream of obstruction 
    % (2nd term is the wake region) [kg/s]
```
- **Purpose**: Predicts flow distribution in wake region
- **Default**: [0.0133, 0.0005]
- **Scope**: Used only by ObstructionSolver; ignored by other solvers

### 4.2 `Solvers/AbstractSolver.m`
**New Property**:
```matlab
SOLVERMODE  (1,1) Solvers.SolverMode  
    % Solver solution mode: NEW, CONTINUE, or SUBSET
```
- **Purpose**: Track solver execution mode
- **Impact**: Need to add to all AbstractSolver subclasses
- **Compatibility**: Non-breaking change (can default to `NEW`)

### 4.3 `Inputs/InputSet.m` (Inferred)
**New Methods** (likely):
- `SplitByAxialPosition()`: Segments InputSet by axial positions
- `SplitBySpanPosition()`: Subdivides InputSet by spanwise positions for obstruction interaction region
- **Purpose**: Enable domain segmentation for obstruction solver

---

## 5. NEW SOLVER ARCHITECTURE PATTERNS

### 5.1 InputSet Segmentation
- **Pattern**: Split InputSet into multiple solution regions
- **Current approach**: Single InputSet → Single Solver instance
- **Obstruction approach**: Single InputSet with obstruction → Multiple InputSets (one per region) → Multiple solver instances (one per SolutionSet)
- **Implication**: Requires modification to how InputSet validation and instantiation works

### 5.2 Domain Decomposition
```
Axial decomposition:
[0] → [axial_start] → [obstruction_start] → [obstruction_end] → [axial_end]
       NONOBS          NONOBS              OBS                  NONOBS

Spanwise decomposition (in OBS region):
[0] → [span_start] → [obstruction_start] → [obstruction_end] → [span_end]
      wake region      obstruction core    wake region
```

---

## 6. INPUT FILE FORMAT CHANGES

### 6.1 New Input File: `obstruction.inp`
- **Format**: Key-value pairs (consistent with existing input files)
- **Expected Sections**:
  - [ID]: Obstruction identifier
  - [SHAPE]: Geometric shape
  - [WALL]: Wall index
  - [RLOCATION]: Position coordinates
  - [DIAMETER]: For circular obstructions
  - [SIDELENGTH], [POINTANGLE]: For triangular obstructions

### 6.2 Existing Input Files
- No breaking changes to existing formats
- Obstruction parameters are optional (not required for non-obstruction simulations)

---

## 7. ENUM COMPARISON: openstream-obs vs openstream

### Present in openstream-obs but missing from openstream:
- `OBSSHAPE`
- `SCBOIL`
- `SOLVER` (needs verification)

### Present in both (may have different content):
- Standard closure model enums (TPFM, VOID, etc.)

---

## 8. CLASS HIERARCHY CHANGES

### Current openstream hierarchy:
```
AbstractSolver
├── Mixture.MixtureSolver
├── TwoFluid.TwoFluidSolver
├── ThreeField.ThreeFieldSolver
└── FourField.FourFieldSolver
```

### After integration:
```
AbstractSolver
├── Mixture.MixtureSolver
├── TwoFluid.TwoFluidSolver
├── ThreeField.ThreeFieldSolver
├── FourField.FourFieldSolver
└── Obstruction.ObstructionSolver (NEW)
```

---

## 9. INTEGRATION CHECKLIST - PHASE 2

**Strategy**: Use openstream as baseline, add obstruction-specific functionality

### Baseline Setup (from openstream)
- [ ] Start with openstream versions of all files
- [ ] Use openstream AbstractSolver.m as base
- [ ] Use openstream Model.m as base
- [ ] Use openstream InputSet.m as base

### Add Obstruction-Specific Classes
- [ ] Copy `Inputs/Obstruction.m` from openstream-obs
- [ ] Copy `InputEnums/OBSSHAPE.m` from openstream-obs
- [ ] Copy `InputEnums/SCBOIL.m` from openstream-obs
- [ ] Copy `Solvers/Obstruction/` directory structure from openstream-obs

### Modify Existing openstream Files (AS NEEDED)
- [ ] Update `Solvers/AbstractSolver.m` to support obstruction features (if required)
- [ ] Add `OBSWSPLITRATIO` to `Inputs/Model.m`
- [ ] Implement `SplitByAxialPosition()` and `SplitBySpanPosition()` in `Inputs/InputSet.m` (if required)
- [ ] Verify/add `InputEnums/SOLVER.m` with OBSTRUCTION option (if required)
- [ ] Verify/add `Solvers/SolverMode.m` (if required)

### Integration Testing
- [ ] Update input file parsing to handle `obstruction.inp`
- [ ] Verify backward compatibility with existing test cases
- [ ] Create/adapt test cases for obstruction solver
- [ ] Test interaction with existing solvers (Mixture, TwoFluid, ThreeField, FourField)

---

## 10. KNOWN ISSUES & NOTES

### Incomplete Implementations in openstream-obs:
- `SCBOIL` enum only has `NONE` - framework only
- ObstructionSolver methods partially implemented (framework structure in place)
- Wake region modeling is simplified

### Potential Integration Challenges:
1. **InputSet architecture**: May require significant refactoring to support domain segmentation
2. **Solver compatibility**: Need to ensure obstruction solver doesn't interfere with existing solvers
3. **Backward compatibility**: All changes must be transparent to existing non-obstruction code
4. **SplitByAxialPosition/SpanPosition methods**: These methods need to be implemented or imported from openstream-obs

### Recommendations:
- Review `InputSet.m` in both versions to understand segmentation approach
- Plan careful unit testing at each integration step
- Consider whether obstruction solver should share more code with Mixture/TwoFluid solvers

---

## 11. TESTING BASELINE

**Current openstream state**:
- Branch: `lee`
- Status: 237 commits ahead of `origin/lee`
- Local changes stashed
- Feature branch created: `feature/obstruction-integration`

**Next Steps (Phase 1 completion)**:
- [ ] Document current passing test cases
- [ ] Create baseline test report
- [ ] Identify integration test cases from openstream-obs

---

**Document Version**: 1.0  
**Last Updated**: April 27, 2026  
**Prepared for**: Feature/obstruction-integration branch
