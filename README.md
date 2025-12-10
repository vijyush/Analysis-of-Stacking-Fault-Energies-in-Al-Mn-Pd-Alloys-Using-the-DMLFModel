# Analysis-of-Stacking-Fault-Energies-in-Al-Mn-Pd-Alloys-Using-the-DMLFModel

# Computational Analysis of Stacking Fault Energies (SFEs) in Al–Mn–Pd Alloys

This repository contains simulation files, automation scripts, and documentation for a full computational study of **Stacking Fault Energies (SFEs)**—Intrinsic (γ_ISF), Extrinsic (γ_ESF), and Twin Fault Energy (γ_Twin)—in **Al–Mn–Pd ternary alloys**.  
The analysis employs **LAMMPS molecular dynamics**, fully automated using Python-based input-file generators.

A total of **189 simulations** were executed across **21 compositions** at **250 K, 450 K, and 600 K**.

---

## 1. Simulation Methodology

### 1.1 Overall Approach
1. Generate FCC(111), HCP(0001), DHCP (ABAC) supercells.  
2. Apply relaxation pipeline:  
   - Conjugate Gradient Minimization  
   - NVT Ensemble (15,000 steps)  
   - NPT Ensemble (30,000 steps)  
3. Extract relaxed total energies:  
   - **E_fcc**, **E_hcp**, **E_dhcp**  
4. Compute SFE values using DMLF equations.  
5. Normalize using the FCC stacking fault area **A_fcc**.

### 1.2 Key Simulation Parameters

| Parameter | Value |
|----------|--------|
| Simulation Engine | LAMMPS |
| Potential | `AlMnPd_Schopf_2012.lammps.EAM.CORRECT` |
| Crystal Structures | FCC(111), HCP(0001), DHCP |
| Temperatures | 250 K, 450 K, 600 K |
| Steps | CG → 15000 (NVT) → 30000 (NPT) |
| Total Simulations | 189 |

---

## 2. DMLF Model – SFE Equations

### Intrinsic SFE (γ_ISF)
```math
\gamma_{ISF} = \frac{4(E_{dhcp} - E_{fcc})}{A_{fcc}}
```

### Extrinsic SFE (γ_ESF)
```math
\gamma_{ESF} = \frac{E_{hcp} + 2E_{dhcp} - 3E_{fcc}}{A_{fcc}}
```

### Twin Fault Energy (γ_Twin)
```math
\gamma_{Twin} = \frac{2(E_{dhcp} - E_{fcc})}{A_{fcc}}
```

---

## 3. Automation Scripts

| Script      | Purpose |
|------------|---------|
| `gen.py`    | Generate LAMMPS input files for **21 ternary** compositions |
| `binary.py` | Generate input files for all **binary** compositions |
| `pure.py`   | Generate input files for **pure elements** (Al, Mn, Pd) |

### Run All Simulations

Windows:
```
run all simulations.bat
```

Linux/macOS:
```
bash run_all_simulations.sh
```

---

# 4. Key Findings (Core Scientific Results)

## 4.1 Composition–SFE Relationship
- **Al-rich compositions** show the **most negative (lowest) SFE values**, enabling easier stacking-fault formation.  
- **Mn + Pd** co-addition increases local lattice distortion → **SFE decreases significantly**.  
- Strong faulting tendency observed in:
  - **Al₄₀Mn₅₀Pd₁₀**
  - **Al₃₃Mn₃₃Pd₃₄**
  - **Al₅₀Mn₂₅Pd₂₅**

## 4.2 Temperature Effects
As temperature increases (250 → 600 K):
- All SFE values become **less negative**.
- Thermal expansion stabilizes FCC → reduced fault formation.
- Ordering of energies consistently holds:
  ```
  γ_ISF < γ_Twin < γ_ESF
  γ_Twin ≈ 2 × γ_ISF
  ```

## 4.3 Metallurgical Implications
### **Low SFE (< 0 mJ/m²):**
- High twinning tendency  
- TRIP/TWIP-like behavior  
- Higher strain hardening  
- Facilitates HCP or DHCP nucleation  

### **High SFE (> 0 mJ/m²):**
- Slip-dominated plasticity  
- Reduced deformation twinning  
- Higher stacking fault resistance  

---

# 5. Requirements

## Python
```
pip install numpy pandas matplotlib scipy
```

## LAMMPS
- Must be installed  
- Must be accessible as `lmp` or `lammps` in your PATH  

---


---

