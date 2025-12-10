 #!/usr/bin/env python3
"""
Generate LAMMPS input files for all compositions and temperatures
9 compositions × 3 temperatures × 3 structures = 81 input files
"""

# Define compositions: (Al%, Mn%, Pd%) → (frac_mn, frac_pd)
compositions = {
    'Al60Mn20Pd20': (0.20, 0.20),
    'Al50Mn40Pd10': (0.40, 0.10),
    'Al40Mn50Pd10': (0.50, 0.10),
    'Al30Mn40Pd30': (0.40, 0.30),
    'Al20Mn50Pd30': (0.50, 0.30),
    'Al10Mn20Pd70': (0.20, 0.70),
    'Al33Mn33Pd34': (0.33, 0.34),
    'Al50Mn25Pd25': (0.25, 0.25),
    'Al45Mn40Pd15': (0.40, 0.15),
}

temperatures = [250, 450, 600]
structures = ['fcc', 'hcp', 'dhcp']

# LAMMPS template for FCC
fcc_template = """# {comp}_{temp}K_fcc.in
# Al-Mn-Pd ternary alloy: {comp_display} at {temp}K
# Crystal structure: FCC (111 orientation)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal 4.15
variable seed equal 12345
variable seed2 equal 54321
variable frac_mn equal {frac_mn}
variable frac_pd equal {frac_pd}
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

# ---------------- Random alloy assignment ----------------
set type 1 type/fraction 2 ${{frac_mn}} ${{seed}}
set type 1 type/fraction 3 ${{frac_pd}} ${{seed2}}

group al type 1
group mn type 2
group pd type 3

variable Ntot equal count(all)
variable Nal equal count(al)
variable Nmn equal count(mn)
variable Npd equal count(pd)

print "========================================"
print "FCC: {comp_display} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "Al = ${{Nal}}, Mn = ${{Nmn}}, Pd = ${{Npd}}"
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

# Lost atoms check
#variable natoms_check equal count(all)
#fix check_atoms all halt 100 v_natoms_check < ${{Ntot}} error

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
write_data {comp}_{temp}K_fcc.data

print "========================================"
print "COMPLETED: {comp_display} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "Check log for warnings"
print "========================================"
"""

# LAMMPS template for HCP
hcp_template = """# {comp}_{temp}K_hcp.in
# Al-Mn-Pd ternary alloy: {comp_display} at {temp}K
# Crystal structure: HCP (0001 orientation)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a equal 2.933
variable c0 equal 7.209
variable seed equal 12345
variable seed2 equal 54321
variable frac_mn equal {frac_mn}
variable frac_pd equal {frac_pd}
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

# ---------------- Random alloy assignment ----------------
set type 1 type/fraction 2 ${{frac_mn}} ${{seed}}
set type 1 type/fraction 3 ${{frac_pd}} ${{seed2}}

group al type 1
group mn type 2
group pd type 3

variable Ntot equal count(all)
variable Nal equal count(al)
variable Nmn equal count(mn)
variable Npd equal count(pd)

print "========================================"
print "HCP: {comp_display} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "Al = ${{Nal}}, Mn = ${{Nmn}}, Pd = ${{Npd}}"
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

#variable natoms_check equal count(all)
#fix check_atoms all halt 100 v_natoms_check < ${{Ntot}} error

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
write_data {comp}_{temp}K_hcp.data

print "========================================"
print "COMPLETED: {comp_display} at {temp}K"
print "Final PE = $(c_Etot) eV"
print "Final PE/atom = $(c_Etot/v_Ntot) eV/atom"
print "Final Lx = $(lx), Ly = $(ly), Lz = $(lz)"
print "Final Area = $(lx*ly) Å²"
print "========================================"
"""

# LAMMPS template for DHCP
dhcp_template = """# {comp}_{temp}K_dhcp.in
# Al-Mn-Pd ternary alloy: {comp_display} at {temp}K
# Crystal structure: DHCP (double HCP with ABAC stacking)

clear
units metal
dimension 3
boundary p p p
atom_style atomic

# ---------------- Parameters ----------------
variable eamfile string "AlMnPd_Schopf_2012.lammps.EAM_CORRECT"
variable a0 equal 2.933
variable c0 equal 7.209
variable seed equal 12345
variable seed2 equal 54321
variable frac_mn equal {frac_mn}
variable frac_pd equal {frac_pd}
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

# ---------------- Random alloy assignment ----------------
set type 1 type/fraction 2 ${{frac_mn}} ${{seed}}
set type 1 type/fraction 3 ${{frac_pd}} ${{seed2}}

group al type 1
group mn type 2
group pd type 3

variable Ntot equal count(all)
variable Nal equal count(al)
variable Nmn equal count(mn)
variable Npd equal count(pd)

print "========================================"
print "DHCP: {comp_display} at {temp}K"
print "Total atoms = ${{Ntot}}"
print "Al = ${{Nal}}, Mn = ${{Nmn}}, Pd = ${{Npd}}"
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

#variable natoms_check equal count(all)
#fix check_atoms all halt 100 v_natoms_check < ${{Ntot}} error

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

unfix check_atoms

# ---------------- Final minimization ----------------
min_style cg
minimize 1.0e-7 1.0e-9 1000 20000

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
    """Generate all LAMMPS input files"""
   
    templates = {
        'fcc': fcc_template,
        'hcp': hcp_template,
        'dhcp': dhcp_template
    }
   
    count = 0
   
    for comp_name, (frac_mn, frac_pd) in compositions.items():
        for temp in temperatures:
            for struct in structures:
                filename = f"{comp_name}_{temp}K_{struct}.in"
               
                # Create composition display string
                al_pct = int((1 - frac_mn - frac_pd) * 100)
                mn_pct = int(frac_mn * 100)
                pd_pct = int(frac_pd * 100)
                comp_display = f"Al-{al_pct}% Mn-{mn_pct}% Pd-{pd_pct}%"
               
                # Fill template
                content = templates[struct].format(
                    comp=comp_name,
                    temp=temp,
                    frac_mn=frac_mn,
                    frac_pd=frac_pd,
                    comp_display=comp_display
                )
               
                # Write file
                with open(filename, 'w') as f:
                    f.write(content)
               
                count += 1
                print(f"Created: {filename}")
   
    print(f"\n{'='*60}")
    print(f"Total files generated: {count}")
    print(f"9 compositions × 3 temperatures × 3 structures = 81 files")
    print(f"{'='*60}")
   
    # Create run script
    create_run_script()

def create_run_script():
    """Create bash script to run all simulations"""
   
    script = """#!/bin/bash
# run_all_simulations.sh
# Automated script to run all 81 LAMMPS simulations

echo "Starting 81 LAMMPS simulations..."
echo "9 compositions × 3 temperatures × 3 structures"
echo ""

# Array of all input files
simulations=("""
   
    for comp_name in compositions.keys():
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
   
    # Write with UTF-8 encoding to handle any special characters
    with open('run_all_simulations.sh', 'w', encoding='utf-8') as f:
        f.write(script)
   
    # Also create Windows batch file
    create_windows_batch()
   
    print("\nCreated: run_all_simulations.sh (for Linux/Mac)")
    print("Created: run_all_simulations.bat (for Windows)")
    print("\nFor Linux/Mac:")
    print("  chmod +x run_all_simulations.sh")
    print("  ./run_all_simulations.sh")
    print("\nFor Windows:")
    print("  run_all_simulations.bat")

def create_windows_batch():
    """Create Windows batch file to run all simulations"""
   
    batch = """@echo off
REM run_all_simulations.bat
REM Automated script to run all 81 LAMMPS simulations on Windows

echo Starting 81 LAMMPS simulations...
echo 9 compositions x 3 temperatures x 3 structures
echo.

set total=0
set completed=0
set failed=0

"""
   
    # Add all simulation files
    for comp_name in compositions.keys():
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
   
    with open('run_all_simulations.bat', 'w', encoding='utf-8') as f:
        f.write(batch)

if __name__ == "__main__":
    print("Generating LAMMPS input files...")
    print("="*60)
    generate_files()
    print("\n" + "="*60)
    print("WINDOWS USERS: Use run_all_simulations.bat")
    print("LINUX/MAC USERS: Use run_all_simulations.sh")
    print("="*60)
    print("\nTo run all simulations:")
    print("  Windows:    run_all_simulations.bat")
    print("  Linux/Mac:  chmod +x run_all_simulations.sh && ./run_all_simulations.sh")
    print("\nOr run individually:")
    print("  lmp -in Al60Mn20Pd20_250K_fcc.in")
