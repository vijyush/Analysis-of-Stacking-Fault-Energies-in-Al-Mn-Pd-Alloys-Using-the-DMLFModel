#!/usr/bin/env python3
"""
Generate LAMMPS input files for 25-75 binary alloys (FAST VERSION)
Reduced timesteps for quicker results
6 binary systems x 3 temperatures x 3 structures = 54 input files
"""

# Define 25-75 binary alloy compositions
# Format: 'Name': (element1, element2, frac_element2)
binary_compositions = {
    'Al25Mn75': ('Al', 'Mn', 0.75),
    'Al75Mn25': ('Al', 'Mn', 0.25),
    'Al25Pd75': ('Al', 'Pd', 0.75),
    'Al75Pd25': ('Al', 'Pd', 0.25),
    'Mn25Pd75': ('Mn', 'Pd', 0.75),
    'Mn75Pd25': ('Mn', 'Pd', 0.25),
}

# Element properties
element_data = {
    'Al': {'mass': 26.981539, 'type': 1},
    'Mn': {'mass': 54.938044, 'type': 2},
    'Pd': {'mass': 106.42, 'type': 3}
}

temperatures = [250, 450, 600]
structures = ['fcc', 'hcp', 'dhcp']

# LAMMPS template for FCC (REDUCED TIMESTEPS)
fcc_template = """# {comp}_{temp}K_fcc.in
# Binary alloy: {comp_display} at {temp}K
# Crystal structure: FCC (111 orientation)
# FAST VERSION - Reduced timesteps

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal 4.05
variable seed equal 12345
variable seed2 equal 54321
variable frac_{elem2_lower} equal {frac_elem2}
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {comp}_{temp}K_fcc.log

# ---------------- Build FCC cell (111 along z) ----------------
lattice fcc ${{a0}} orient x 1 -1 0 orient y 1 1 -2 orient z 1 1 1
region box block 0 6 0 6 0 8
create_box 3 box
create_atoms 1 box

# ---------------- Masses ----------------
mass 1 26.981539
mass 2 54.938044
mass 3 106.42

# ---------------- Random binary alloy assignment ----------------
# Start with all atoms as {elem1} (type {type1})
set type * type {type1}

# Randomly assign {elem2} (type {type2})
set type {type1} type/fraction {type2} ${{frac_{elem2_lower}}} ${{seed}}

group {elem1_lower} type {type1}
group {elem2_lower} type {type2}

variable Ntot equal count(all)
variable N{elem1_lower} equal count({elem1_lower})
variable N{elem2_lower} equal count({elem2_lower})

print "========================================"
print "FCC: {comp_display} at {temp}K (FAST)"
print "Total atoms = ${{Ntot}}"
print "{elem1} = ${{N{elem1_lower}}}, {elem2} = ${{N{elem2_lower}}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} Al Mn Pd

# ---------------- Neighbor settings ----------------
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes

# ---------------- Diagnostics ----------------
compute peratom all pe/atom
compute Etot all pe
compute Kall all ke

thermo 50
thermo_style custom step temp pe ke etotal press lx ly lz vol
thermo_modify flush yes lost warn

# ---------------- Relaxation protocol (REDUCED) ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 500 10000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 1000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Heating to Tfinal (REDUCED) ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 5000
unfix heat

# ---------------- NPT equilibration (REDUCED) ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) iso 0.0 0.0 $(1000.0*dt)
run 10000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Final outputs ----------------
write_data {comp}_{temp}K_fcc.data

print "========================================"
print "COMPLETED: {comp_display} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "========================================"
"""

# LAMMPS template for HCP (REDUCED TIMESTEPS)
hcp_template = """# {comp}_{temp}K_hcp.in
# Binary alloy: {comp_display} at {temp}K
# Crystal structure: HCP (0001 orientation)
# FAST VERSION - Reduced timesteps

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a equal 2.863
variable c0 equal 4.663
variable seed equal 12345
variable seed2 equal 54321
variable frac_{elem2_lower} equal {frac_elem2}
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {comp}_{temp}K_hcp.log

# ---------------- Build HCP cell (0001 along z) ----------------
lattice hcp ${{a}} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 16
create_box 3 box
create_atoms 1 box

# ---------------- Masses ----------------
mass 1 26.981539
mass 2 54.938044
mass 3 106.42

# ---------------- Random binary alloy assignment ----------------
set type * type {type1}
set type {type1} type/fraction {type2} ${{frac_{elem2_lower}}} ${{seed}}

group {elem1_lower} type {type1}
group {elem2_lower} type {type2}

variable Ntot equal count(all)
variable N{elem1_lower} equal count({elem1_lower})
variable N{elem2_lower} equal count({elem2_lower})

print "========================================"
print "HCP: {comp_display} at {temp}K (FAST)"
print "Total atoms = ${{Ntot}}"
print "{elem1} = ${{N{elem1_lower}}}, {elem2} = ${{N{elem2_lower}}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} Al Mn Pd

# ---------------- Neighbor settings ----------------
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes

# ---------------- Diagnostics ----------------
compute peratom all pe/atom
compute Etot all pe
compute Kall all ke

thermo 50
thermo_style custom step temp pe ke etotal press lx ly lz vol
thermo_modify flush yes lost warn

# ---------------- Relaxation protocol (REDUCED) ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 500 10000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 1000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Heating to Tfinal (REDUCED) ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 5000
unfix heat

# ---------------- NPT equilibration (REDUCED) ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) aniso 0.0 0.0 $(1000.0*dt)
run 10000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Final outputs ----------------
write_data {comp}_{temp}K_hcp.data

print "========================================"
print "COMPLETED: {comp_display} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "========================================"
"""

# LAMMPS template for DHCP (REDUCED TIMESTEPS)
dhcp_template = """# {comp}_{temp}K_dhcp.in
# Binary alloy: {comp_display} at {temp}K
# Crystal structure: DHCP (double HCP with ABAC stacking)
# FAST VERSION - Reduced timesteps

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal 2.863
variable c0 equal 4.663
variable seed equal 12345
variable seed2 equal 54321
variable frac_{elem2_lower} equal {frac_elem2}
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {comp}_{temp}K_dhcp.log

# ---------------- Build DHCP structure ----------------
lattice hcp ${{a0}} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 8
create_box 3 box
create_atoms 1 box

# Create DHCP stacking fault
variable layer_height equal ${{c0}}/4.0
variable shift_x equal ${{a0}}/3.0
variable shift_y equal ${{a0}}/(3.0*sqrt(3.0))

region layer1 block INF INF INF INF $(3.0*v_layer_height) $(4.0*v_layer_height)
region layer2 block INF INF INF INF $(7.0*v_layer_height) $(8.0*v_layer_height)

group shift_layers region layer1
group shift_layers region layer2

displace_atoms shift_layers move ${{shift_x}} ${{shift_y}} 0.0 units box

# ---------------- Masses ----------------
mass 1 26.981539
mass 2 54.938044
mass 3 106.42

# ---------------- Random binary alloy assignment ----------------
set type * type {type1}
set type {type1} type/fraction {type2} ${{frac_{elem2_lower}}} ${{seed}}

group {elem1_lower} type {type1}
group {elem2_lower} type {type2}

variable Ntot equal count(all)
variable N{elem1_lower} equal count({elem1_lower})
variable N{elem2_lower} equal count({elem2_lower})

print "========================================"
print "DHCP: {comp_display} at {temp}K (FAST)"
print "Total atoms = ${{Ntot}}"
print "{elem1} = ${{N{elem1_lower}}}, {elem2} = ${{N{elem2_lower}}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} Al Mn Pd

# ---------------- Neighbor settings ----------------
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes

# ---------------- Diagnostics ----------------
compute peratom all pe/atom
compute Etot all pe
compute Kall all ke

thermo 50
thermo_style custom step temp pe ke etotal press lx ly lz vol
thermo_modify flush yes lost warn

# ---------------- Relaxation protocol (REDUCED) ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 500 10000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 1000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Heating to Tfinal (REDUCED) ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 5000
unfix heat

# ---------------- NPT equilibration (REDUCED) ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) aniso 0.0 0.0 $(1000.0*dt)
run 10000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 500 10000

# ---------------- Final outputs ----------------
write_data {comp}_{temp}K_dhcp.data

print "========================================"
print "COMPLETED: {comp_display} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "SFE = (E_dhcp - E_hcp) / Area"
print "========================================"
"""

def generate_files():
    """Generate all LAMMPS input files for 25-75 binary alloys"""
   
    templates = {
        'fcc': fcc_template,
        'hcp': hcp_template,
        'dhcp': dhcp_template
    }
   
    count = 0
   
    for comp_name, (elem1, elem2, frac_elem2) in binary_compositions.items():
        for temp in temperatures:
            for struct in structures:
                filename = f"{comp_name}_{temp}K_{struct}.in"
               
                # Create composition display string
                frac_elem1 = 1.0 - frac_elem2
                frac1_percent = int(frac_elem1 * 100)
                frac2_percent = int(frac_elem2 * 100)
                comp_display = f"{elem1}-{frac1_percent}% {elem2}-{frac2_percent}%"
               
                # Get element types
                type1 = element_data[elem1]['type']
                type2 = element_data[elem2]['type']
               
                # Fill template
                content = templates[struct].format(
                    comp=comp_name,
                    temp=temp,
                    frac_elem2=frac_elem2,
                    comp_display=comp_display,
                    elem1=elem1,
                    elem2=elem2,
                    elem1_lower=elem1.lower(),
                    elem2_lower=elem2.lower(),
                    type1=type1,
                    type2=type2
                )
               
                # Write file
                with open(filename, 'w') as f:
                    f.write(content)
               
                count += 1
                print(f"Created: {filename}")
   
    print(f"\n{'='*60}")
    print(f"Total files generated: {count}")
    print(f"6 binary systems (25-75) × 3 temperatures × 3 structures = 54 files")
    print(f"FAST VERSION - Reduced timesteps for quicker results")
    print(f"{'='*60}")
   
    # Create run scripts
    create_run_script()

def create_run_script():
    """Create bash script to run all simulations"""
   
    script = """#!/bin/bash
# run_binary_25_75.sh
# Automated script to run all 25-75 binary alloy LAMMPS simulations
# FAST VERSION - Reduced timesteps

echo "Starting 25-75 binary alloy LAMMPS simulations (FAST)..."
echo "6 binary systems × 3 temperatures × 3 structures = 54 files"
echo ""

# Array of all input files
simulations=("""
   
    for comp_name in sorted(binary_compositions.keys()):
        for temp in temperatures:
            for struct in structures:
                script += f'\n    "{comp_name}_{temp}K_{struct}.in"'
   
    script += """
)

# Run each simulation
total=${#simulations[@]}
current=0

for sim in "${simulations[@]}"; do
    current=$((current + 1))
    echo "================================================"
    echo "Running simulation $current/$total: $sim"
    echo "================================================"
   
    # Run LAMMPS (adjust path to your LAMMPS executable)
    lmp -in $sim
   
    # Check if successful
    if [ $? -eq 0 ]; then
        echo "[SUCCESS] $sim completed"
    else
        echo "[FAILED] $sim had errors"
    fi
    echo ""
done

echo "================================================"
echo "All simulations completed!"
echo "================================================"
"""
   
    with open('run_binary_25_75.sh', 'w', encoding='utf-8') as f:
        f.write(script)
   
    # Also create Windows batch file
    create_windows_batch()
   
    print("\nCreated: run_binary_25_75.sh (for Linux/Mac)")
    print("Created: run_binary_25_75.bat (for Windows)")

def create_windows_batch():
    """Create Windows batch file to run all simulations"""
   
    batch = """@echo off
REM run_binary_25_75.bat
REM Automated script to run all 25-75 binary alloy LAMMPS simulations on Windows
REM FAST VERSION - Reduced timesteps

echo Starting 25-75 binary alloy LAMMPS simulations (FAST)...
echo 6 binary systems x 3 temperatures x 3 structures = 54 files
echo.

set total=0
set completed=0
set failed=0

"""
   
    # Add all simulation files
    for comp_name in sorted(binary_compositions.keys()):
        for temp in temperatures:
            for struct in structures:
                batch += f"""
echo ================================================
echo Running: {comp_name}_{temp}K_{struct}.in
echo ================================================
lmp -in {comp_name}_{temp}K_{struct}.in
if %ERRORLEVEL% EQU 0 (
    echo [SUCCESS] {comp_name}_{temp}K_{struct}.in completed
    set /a completed+=1
) else (
    echo [FAILED] {comp_name}_{temp}K_{struct}.in had errors
    set /a failed+=1
)
set /a total+=1
echo.

"""
   
    batch += """
echo ================================================
echo All simulations completed!
echo Total: %total%, Successful: %completed%, Failed: %failed%
echo ================================================
pause
"""
   
    with open('run_binary_25_75.bat', 'w', encoding='utf-8') as f:
        f.write(batch)

if __name__ == "__main__":
    print("Generating LAMMPS input files for 25-75 binary alloys (FAST VERSION)...")
    print("="*60)
    generate_files()
    print("\n" + "="*60)
    print("WINDOWS USERS: Use run_binary_25_75.bat")
    print("LINUX/MAC USERS: Use run_binary_25_75.sh")
    print("="*60)
    print("\nTimestep reductions (approximate speedup: 2-3x faster):")
    print("  - Initial NVE/limit: 2000 → 1000 steps")
    print("  - Heating (NVT): 15000 → 5000 steps")
    print("  - NPT equilibration: 30000 → 10000 steps")
    print("  - Minimization steps: 1000→500, 20000→10000")
    print("\n25-75 Binary systems:")
    print("  Al25Mn75, Al75Mn25, Al25Pd75, Al75Pd25, Mn25Pd75, Mn75Pd25")
    print("\nTo run all simulations:")
    print("  Windows:    run_binary_25_75.bat")
    print("  Linux/Mac:  chmod +x run_binary_25_75.sh && ./run_binary_25_75.sh")
