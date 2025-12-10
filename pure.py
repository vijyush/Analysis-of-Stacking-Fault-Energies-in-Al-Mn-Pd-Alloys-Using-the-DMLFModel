#!/usr/bin/env python3
"""
Generate LAMMPS input files for pure elements
3 elements × 3 temperatures × 3 structures = 27 input files
"""

# Define pure elements: element_name → (frac_element, mass, lattice_param)
elements = {
    'Al': {
        'mass': 26.981539,
        'fcc_a0': 4.05,
        'hcp_a': 2.863,
        'hcp_c': 4.663,
    },
    'Mn': {
        'mass': 54.938044,
        'fcc_a0': 3.72,
        'hcp_a': 2.63,
        'hcp_c': 4.29,
    },
    'Pd': {
        'mass': 106.42,
        'fcc_a0': 3.89,
        'hcp_a': 2.75,
        'hcp_c': 4.49,
    }
}

temperatures = [250, 450, 600]
structures = ['fcc', 'hcp', 'dhcp']

# LAMMPS template for FCC
fcc_template = """# {element}_{temp}K_fcc.in
# Pure {element} at {temp}K
# Crystal structure: FCC (111 orientation)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal {fcc_a0}
variable seed equal 12345
variable seed2 equal 54321
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {element}_{temp}K_fcc.log

# ---------------- Build FCC cell (111 along z) ----------------
lattice fcc ${{a0}} orient x 1 -1 0 orient y 1 1 -2 orient z 1 1 1
region box block 0 6 0 6 0 8
create_box 1 box
create_atoms 1 box

# ---------------- Mass ----------------
mass 1 {mass}

group {element_lower} type 1

variable Ntot equal count(all)

print "========================================"
print "FCC: Pure {element} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} {element}

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

# ---------------- Relaxation protocol ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 1000 20000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 2000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Heating to Tfinal ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 15000
unfix heat

# ---------------- NPT equilibration ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) iso 0.0 0.0 $(1000.0*dt)
run 30000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Final outputs ----------------
write_data {element}_{temp}K_fcc.data

print "========================================"
print "COMPLETED: Pure {element} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "========================================"
"""

# LAMMPS template for HCP
hcp_template = """# {element}_{temp}K_hcp.in
# Pure {element} at {temp}K
# Crystal structure: HCP (0001 orientation)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a equal {hcp_a}
variable c0 equal {hcp_c}
variable seed equal 12345
variable seed2 equal 54321
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {element}_{temp}K_hcp.log

# ---------------- Build HCP cell (0001 along z) ----------------
lattice hcp ${{a}} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 16
create_box 1 box
create_atoms 1 box

# ---------------- Mass ----------------
mass 1 {mass}

group {element_lower} type 1

variable Ntot equal count(all)

print "========================================"
print "HCP: Pure {element} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} {element}

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

# ---------------- Relaxation protocol ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 1000 20000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 2000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Heating to Tfinal ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 15000
unfix heat

# ---------------- NPT equilibration ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) aniso 0.0 0.0 $(1000.0*dt)
run 30000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Final outputs ----------------
write_data {element}_{temp}K_hcp.data

print "========================================"
print "COMPLETED: Pure {element} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "========================================"
"""

# LAMMPS template for DHCP
dhcp_template = """# {element}_{temp}K_dhcp.in
# Pure {element} at {temp}K
# Crystal structure: DHCP (double HCP with ABAC stacking)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal {hcp_a}
variable c0 equal {hcp_c}
variable seed equal 12345
variable seed2 equal 54321
variable Tfinal equal {temp}.0
variable dt equal 0.001

# ---------------- Log file ----------------
log {element}_{temp}K_dhcp.log

# ---------------- Build DHCP structure ----------------
lattice hcp ${{a0}} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 8
create_box 1 box
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

# ---------------- Mass ----------------
mass 1 {mass}

group {element_lower} type 1

variable Ntot equal count(all)

print "========================================"
print "DHCP: Pure {element} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "========================================"

# ---------------- Potential ----------------
pair_style eam/alloy
pair_coeff * * ${{eamfile}} {element}

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

# ---------------- Relaxation protocol ----------------
pair_modify shift yes
min_style cg
minimize 1.0e-5 1.0e-7 1000 20000

timestep ${{dt}}
velocity all create 10.0 ${{seed}} dist gaussian mom yes rot yes
fix limiter all nve/limit 0.01
run 2000
unfix limiter

min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Heating to Tfinal ----------------
velocity all create 10.0 ${{seed2}} dist gaussian mom yes rot yes
fix heat all nvt temp 10.0 ${{Tfinal}} $(100.0*dt)
run 15000
unfix heat

# ---------------- NPT equilibration ----------------
fix press all npt temp ${{Tfinal}} ${{Tfinal}} $(100.0*dt) aniso 0.0 0.0 $(1000.0*dt)
run 30000
unfix press

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

# ---------------- Final outputs ----------------
write_data {element}_{temp}K_dhcp.data

print "========================================"
print "COMPLETED: Pure {element} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "SFE = (E_dhcp - E_hcp) / Area"
print "========================================"
"""

def generate_files():
    """Generate all LAMMPS input files for pure elements"""
   
    templates = {
        'fcc': fcc_template,
        'hcp': hcp_template,
        'dhcp': dhcp_template
    }
   
    count = 0
   
    for element, params in elements.items():
        for temp in temperatures:
            for struct in structures:
                filename = f"{element}_{temp}K_{struct}.in"
               
                # Fill template
                content = templates[struct].format(
                    element=element,
                    element_lower=element.lower(),
                    temp=temp,
                    mass=params['mass'],
                    fcc_a0=params['fcc_a0'],
                    hcp_a=params['hcp_a'],
                    hcp_c=params['hcp_c']
                )
               
                # Write file
                with open(filename, 'w') as f:
                    f.write(content)
               
                count += 1
                print(f"Created: {filename}")
   
    print(f"\n{'='*60}")
    print(f"Total files generated: {count}")
    print(f"3 elements × 3 temperatures × 3 structures = 27 files")
    print(f"{'='*60}")
   
    # Create run script
    create_run_script()

def create_run_script():
    """Create bash script to run all simulations"""
   
    script = """#!/bin/bash
# run_all_pure_elements.sh
# Automated script to run all 27 pure element LAMMPS simulations

echo "Starting 27 pure element LAMMPS simulations..."
echo "3 elements × 3 temperatures × 3 structures"
echo ""

# Array of all input files
simulations=("""
   
    for element in elements.keys():
        for temp in temperatures:
            for struct in structures:
                script += f'\n    "{element}_{temp}K_{struct}.in"'
   
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
   
    with open('run_all_pure_elements.sh', 'w', encoding='utf-8') as f:
        f.write(script)
   
    # Also create Windows batch file
    create_windows_batch()
   
    print("\nCreated: run_all_pure_elements.sh (for Linux/Mac)")
    print("Created: run_all_pure_elements.bat (for Windows)")
    print("\nFor Linux/Mac:")
    print("  chmod +x run_all_pure_elements.sh")
    print("  ./run_all_pure_elements.sh")
    print("\nFor Windows:")
    print("  run_all_pure_elements.bat")

def create_windows_batch():
    """Create Windows batch file to run all simulations"""
   
    batch = """@echo off
REM run_all_pure_elements.bat
REM Automated script to run all 27 pure element LAMMPS simulations on Windows

echo Starting 27 pure element LAMMPS simulations...
echo 3 elements x 3 temperatures x 3 structures
echo.

set total=0
set completed=0
set failed=0

"""
   
    # Add all simulation files
    for element in elements.keys():
        for temp in temperatures:
            for struct in structures:
                batch += f"""
echo ================================================
echo Running: {element}_{temp}K_{struct}.in
echo ================================================
lmp -in {element}_{temp}K_{struct}.in
if %ERRORLEVEL% EQU 0 (
    echo [SUCCESS] {element}_{temp}K_{struct}.in completed
    set /a completed+=1
) else (
    echo [FAILED] {element}_{temp}K_{struct}.in had errors
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
   
    with open('run_all_pure_elements.bat', 'w', encoding='utf-8') as f:
        f.write(batch)

if __name__ == "__main__":
    print("Generating LAMMPS input files for pure elements...")
    print("="*60)
    generate_files()
    print("\n" + "="*60)
    print("WINDOWS USERS: Use run_all_pure_elements.bat")
    print("LINUX/MAC USERS: Use run_all_pure_elements.sh")
    print("="*60)
    print("\nTo run all simulations:")
    print("  Windows:    run_all_pure_elements.bat")
    print("  Linux/Mac:  chmod +x run_all_pure_elements.sh && ./run_all_pure_elements.sh")
    print("\nOr run individually:")
    print("  lmp -in Al_250K_fcc.in")
    print("  lmp -in Mn_450K_hcp.in")
    print("  lmp -in Pd_600K_dhcp.in")
