#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
from io import StringIO

# Paste your data here (triple quotes)
data = """Composition,Temp (K),Phase,Final PE (eV),Final Lz (Å),Final Area (Å²),Number of Atoms,a,c
Al20Mn50Pd30,600,fcc,-7176.933179,39.431084,659.031359,486,4.15,
Al20Mn50Pd30,600,hcp,-20012.06459,68.920818,486.894881,257,2.933,7.209
Al30Mn40Pd30,250,dhcp,-19574.7392,63.614649,381.174456,1152,2.933,7.209
Al30Mn40Pd30,250,fcc,-28578.38295,54.734862,1297.264959,4608,4.15,
Al30Mn40Pd30,250,hcp,-14520.28558,72.494437,465.08353,2304,2.933,7.209
Al30Mn40Pd30,450,dhcp,-19115.80432,62.25685,381.704136,1152,2.933,7.209
Al30Mn40Pd30,450,fcc,-28651.03608,54.725381,1296.815582,4608,4.15,
Al30Mn40Pd30,450,hcp,-14550.33421,69.975151,480.881792,2304,2.933,7.209
Al30Mn40Pd30,600,dhcp,-18885.79614,61.566771,381.215276,1152,2.933,7.209
Al30Mn40Pd30,600,fcc,-28753.57826,54.104105,1267.538232,4608,4.15,
Al30Mn40Pd30,600,hcp,-12023.63971,60.598497,409.147722,300,2.933,7.209
Al33Mn33Pd34,250,dhcp,-19237.32927,62.357269,358.319497,1152,2.933,7.209
Al33Mn33Pd34,250,fcc,-27360.45901,53.017446,1217.13357,4608,4.15,
Al33Mn33Pd34,250,hcp,-13883.61452,71.682403,499.503423,2304,2.933,7.209
Al33Mn33Pd34,450,dhcp,-18778.39439,60.99947,358.849177,1152,2.933,7.209
Al33Mn33Pd34,450,fcc,-27402.19401,52.631113,1214.997531,4608,4.15,
Al33Mn33Pd34,450,hcp,-13905.7238,70.276444,489.128742,2304,2.933,7.209
Al33Mn33Pd34,600,dhcp,-18548.38621,60.309391,358.360317,1152,2.933,7.209
Al33Mn33Pd34,600,fcc,-27420.35324,52.587302,1204.668729,4608,4.15,
Al33Mn33Pd34,600,hcp,-13930.5302,69.070451,501.305043,2304,2.933,7.209
Al40Mn50Pd10,250,dhcp,-16826.19694,63.740031,381.508062,1152,2.933,7.209
Al40Mn50Pd10,250,fcc,-30175.95636,53.612237,1214.90652,4608,4.15,
Al40Mn50Pd10,250,hcp,-15112.89711,69.830712,490.905013,2304,2.933,7.209
Al40Mn50Pd10,450,dhcp,-16367.26206,62.382152,382.037742,1152,2.933,7.209
Al40Mn50Pd10,450,fcc,-30265.69612,53.184641,1278.51013,4608,4.15,
Al40Mn50Pd10,450,hcp,-15145.32927,69.384387,489.718821,2304,2.933,7.209
Al40Mn50Pd10,600,dhcp,-16137.25388,61.691945,381.54936,1152,2.933,7.209
Al40Mn50Pd10,600,fcc,-30338.32487,53.183676,1219.091853,4608,4.15,
Al40Mn50Pd10,600,hcp,-15178.89161,69.908281,489.510787,2304,2.933,7.209
Al45Mn40Pd15,250,dhcp,-17148.58362,63.512025,378.500558,1152,2.933,7.209
Al45Mn40Pd15,250,fcc,-28875.09038,53.449346,1282.172133,4608,4.15,
Al45Mn40Pd15,250,hcp,-14757.23625,69.735704,478.245797,2304,2.933,7.209
Al45Mn40Pd15,450,dhcp,-16689.64874,62.154146,379.030238,1152,2.933,7.209
Al45Mn40Pd15,450,fcc,-28943.91016,53.394981,1279.888135,4608,4.15,
Al45Mn40Pd15,450,hcp,-14805.09209,69.719834,480.689585,2304,2.933,7.209
Al45Mn40Pd15,600,dhcp,-16459.64057,61.463957,378.542376,1152,2.933,7.209
Al45Mn40Pd15,600,fcc,-28975.90019,53.468822,1296.225009,4608,4.15,
Al45Mn40Pd15,600,hcp,-14836.42138,68.287494,465.846129,2304,2.933,7.209
Al50Mn25Pd25,250,dhcp,-18044.14015,63.30745,374.801523,1152,2.933,7.209
Al50Mn25Pd25,250,fcc,-28427.47807,53.163842,1265.335122,4608,4.15,
Al50Mn25Pd25,250,hcp,-14374.2482,71.818924,488.701904,2304,2.933,7.209
Al50Mn25Pd25,450,dhcp,-17585.20528,61.94957,375.331203,1152,2.933,7.209
Al50Mn25Pd25,450,fcc,-28534.36877,53.256917,1290.247863,4608,4.15,
Al50Mn25Pd25,450,hcp,-14396.56626,69.29927,490.374782,2304,2.933,7.209
Al50Mn25Pd25,600,dhcp,-17355.1971,61.259391,374.843343,1152,2.933,7.209
Al50Mn25Pd25,600,fcc,-28600.19958,53.109745,1272.130235,4608,4.15,
Al50Mn25Pd25,600,hcp,-14425.60164,69.33219,482.264141,2304,2.933,7.209
Al50Mn40Pd10,250,dhcp,-18190.48359,63.433468,377.430283,1152,2.933,7.209
Al50Mn40Pd10,250,fcc,-28526.89587,52.835453,1248.836351,4608,4.15,
Al50Mn40Pd10,250,hcp,-14582.99199,74.196578,482.520395,2304,2.933,7.209
Al50Mn40Pd10,450,dhcp,-17731.54872,62.075589,377.96,1152,2.933,7.209
Al50Mn40Pd10,450,fcc,-28181.24987,53.258515,1228.227302,4608,4.15,
Al50Mn40Pd10,450,hcp,-14373.93923,63.785147,513.035247,2304,2.933,7.209
Al50Mn40Pd10,600,dhcp,-17501.54054,61.38541,377.472119,1152,2.933,7.209
Al50Mn40Pd10,600,fcc,-28232.38028,52.920691,1212.695155,4608,4.15,
Al50Mn40Pd10,600,hcp,-14417.55691,59.453407,522.833145,2304,2.933,7.209
Al60Mn20Pd20,250,dhcp,-15783.23779,36.381522,87.214896,137,2.933,7.209
Al60Mn20Pd20,250,fcc,-24196.8634,54.585042,1290.172943,4608,4.15,
Al60Mn20Pd20,250,hcp,-12307.98278,72.129986,480.400964,2304,2.933,7.209
Al60Mn20Pd20,450,dhcp,-15324.30292,36.006897,90.059462,151,2.933,7.209
Al60Mn20Pd20,450,fcc,-24247.44702,54.678692,1294.603784,4608,4.15,
Al60Mn20Pd20,450,hcp,-12317.61827,72.82539,478.778504,2304,2.933,7.209
Al60Mn20Pd20,600,dhcp,-15094.29474,35.69732,91.896519,168,2.933,7.209
Al60Mn20Pd20,600,fcc,-24270.00764,54.717123,1296.424249,4608,4.15,
Al60Mn20Pd20,600,hcp,-12334.95706,71.554428,491.266139,2304,2.933,7.209

"""

# Define column names
columns = [
    "Composition", "Temperature(K)", "Phase",
    "Potential_Energy", "Force", "Volume", "Atoms",
    "a(Å)", "c(Å)"
]

# Read data into a DataFrame
df = pd.read_csv(StringIO(data), names=columns)

# Save as Excel file
df.to_excel("alloy_data2.xlsx", index=False)

print("✅ Excel file 'alloy_data2.xlsx' has been created successfully!")


# In[8]:


import pandas as pd
import numpy as np

# ---------- Load Excel ----------
file_path = r"C:\Users\burra\Downloads\mini-project(compu)\alloy_data.xlsx"   # uploaded file name
df = pd.read_excel(file_path)

# Ensure column names match your sheet exactly
# Must include: 'Phase', 'a', 'c'

# ---------- Calculate interlayer spacing ----------
def calc_d(row):
    phase = row["Phase"].lower()
    a = row["a"]
    c = row["c"]

    if phase == "fcc":
        return a / np.sqrt(3)
    elif phase == "hcp":
        return c / 2
    elif phase == "dhcp":
        # if c same as hcp, this gives layer spacing per plane
        return c / 4
    else:
        return np.nan

df["d_spacing (Å)"] = df.apply(calc_d, axis=1)

# ---------- Save output ----------
df.to_excel("alloy_with_dspacing.xlsx", index=False)

print("✅ Interlayer spacing calculated and saved to 'alloy_with_dspacing.xlsx'")
print(df[["Composition", "Temp (K)", "Phase", "a", "c", "d_spacing (Å)"]].head())


# In[9]:


import pandas as pd

# Load the normalized energy data
file_path = r"C:\Users\burra\Downloads\mini-project(compu)\final_scaled_energy.xlsx"
df = pd.read_excel(file_path)

# Create new columns
df["γISF (eV/Å²)"] = None
df["γESF (eV/Å²)"] = None
df["γTwin (eV/Å²)"] = None

# Calculate stacking fault energies
for (comp, temp), group in df.groupby(["Composition", "Temp (K)"]):
    try:
        # Extract rows for each phase
        E_fcc = group[group["Phase"] == "fcc"]["E_scaled (eV)"].values[0]
        E_hcp = group[group["Phase"] == "hcp"]["E_scaled (eV)"].values[0]
        E_dhcp = group[group["Phase"] == "dhcp"]["E_scaled (eV)"].values[0]
        A_fcc = group[group["Phase"] == "fcc"]["Final Area (Å²)"].values[0]

        # Formulas from your image
        γISF = 4 * (E_dhcp - E_fcc) / A_fcc
        γESF = (E_hcp + 2 * E_dhcp - 3 * E_fcc) / A_fcc
        γTwin = 2 * (E_dhcp - E_fcc) / A_fcc

        # Assign results to all matching rows
        for idx in group.index:
            df.loc[idx, "γISF (eV/Å²)"] = γISF
            df.loc[idx, "γESF (eV/Å²)"] = γESF
            df.loc[idx, "γTwin (eV/Å²)"] = γTwin
    except:
        continue

# Convert from eV/Å² → mJ/m²  (1 eV/Å² = 16021.77 mJ/m²)
df["γISF (mJ/m²)"] = df["γISF (eV/Å²)"] * 16021.77
df["γESF (mJ/m²)"] = df["γESF (eV/Å²)"] * 16021.77
df["γTwin (mJ/m²)"] = df["γTwin (eV/Å²)"] * 16021.77

# Save to new Excel
output_path = "SFE_results.xlsx"
df.to_excel(output_path, index=False)

print("✅ Stacking fault energies calculated and saved to:", output_path)


# In[15]:


pip install pandas matplotlib mpltern scipy


# In[27]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------
# STEP 1: Load the SFE results file
# -----------------------------------
df = pd.read_excel(r"C:\Users\burra\Downloads\mini-project(compu)\SFE_results.xlsx")

# Optional: fix large unit scale if needed
df["γISF (mJ/m²)"] = df["γISF (mJ/m²)"] / 1e6
df["γESF (mJ/m²)"] = df["γESF (mJ/m²)"] / 1e6
df["γTwin (mJ/m²)"] = df["γTwin (mJ/m²)"] / 1e6

# -----------------------------------
# STEP 2: Set the style
# -----------------------------------
sns.set(style="whitegrid", font_scale=1.1)

# Define the SFE types for plotting
sfe_types = [
    ("γISF (mJ/m²)", "Intrinsic Stacking Fault Energy (γISF)"),
    ("γESF (mJ/m²)", "Extrinsic Stacking Fault Energy (γESF)"),
    ("γTwin (mJ/m²)", "Twin Fault Energy (γTwin)")
]

# -----------------------------------
# STEP 3: Generate line plots
# -----------------------------------
for col, label in sfe_types:
    plt.figure(figsize=(8, 6))

    # Plot SFE vs Temperature for each composition
    sns.lineplot(
        data=df,
        x="Temp (K)",
        y=col,
        hue="Composition",
        marker="o",
        linewidth=2
    )

    plt.title(f"{label} vs Temperature", fontsize=14, fontweight='bold')
    plt.xlabel("Temperature (K)", fontsize=12)
    plt.ylabel("SFE (mJ/m²)", fontsize=12)
    plt.legend(title="Composition", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    plt.tight_layout()

    filename = f"SFE_LinePlot_{col.replace(' (mJ/m²)', '')}.png"
    plt.savefig(filename, dpi=400)
    plt.close()
    print(f"✅ Saved: {filename}")

print("\n🎯 All SFE vs Temperature line plots generated successfully!")


# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ---------------------------------------------
# STEP 1: Load your SFE results Excel file
# ---------------------------------------------
df = pd.read_excel(r"C:\Users\burra\Downloads\mini-project(compu)\SFE_results.xlsx")

# ---------------------------------------------
# STEP 2: Optional scaling (if values too large)
# ---------------------------------------------
df["γISF (mJ/m²)"] = df["γISF (mJ/m²)"] / 1e6
df["γESF (mJ/m²)"] = df["γESF (mJ/m²)"] / 1e6
df["γTwin (mJ/m²)"] = df["γTwin (mJ/m²)"] / 1e6

# ---------------------------------------------
# STEP 3: Seaborn style setup
# ---------------------------------------------
sns.set(style="whitegrid", font_scale=1.1, rc={"figure.dpi": 120})
palette = sns.color_palette("coolwarm", 3)

# ---------------------------------------------
# STEP 4: Generate bar plots for each temperature
# ---------------------------------------------
temps = sorted(df["Temp (K)"].unique())

for T in temps:
    df_T = df[df["Temp (K)"] == T]

    # Melt for grouped bar plotting
    df_melt = df_T.melt(
        id_vars=["Composition"],
        value_vars=["γISF (mJ/m²)", "γESF (mJ/m²)", "γTwin (mJ/m²)"],
        var_name="Fault Type",
        value_name="SFE (mJ/m²)"
    )

    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=df_melt,
        x="Composition",
        y="SFE (mJ/m²)",
        hue="Fault Type",
        palette=palette,
        edgecolor="black"
    )

    plt.title(f"Stacking Fault Energies at {T} K", fontsize=14, fontweight="bold")
    plt.xlabel("Composition", fontsize=12)
    plt.ylabel("SFE (mJ/m²)", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Fault Type", loc="upper right")
    plt.tight_layout()

    filename = f"SFE_BarPlot_{T}K.png"
    plt.savefig(filename, dpi=400)
    plt.close()

    print(f"✅ Saved: {filename}")

# ---------------------------------------------
# STEP 5: Summary
# ---------------------------------------------
print("\n🎯 All 3 bar plots generated successfully!")
print(f"📁 Saved in directory: {os.getcwd()}")


# In[3]:


import pandas as pd
import matplotlib.pyplot as plt
import mpltern
import numpy as np
from scipy.interpolate import griddata
import re
import os

# -------------------------------
# STEP 1: Load Excel file
# -------------------------------
file_path = r"C:\Users\burra\Downloads\mini-project(compu)\SFE_results.xlsx"
df = pd.read_excel(file_path)

# -------------------------------
# STEP 2: Fix scaling (realistic mJ/m² range)
# -------------------------------
df["γISF (mJ/m²)"] = df["γISF (mJ/m²)"] / 1e6
df["γESF (mJ/m²)"] = df["γESF (mJ/m²)"] / 1e6
df["γTwin (mJ/m²)"] = df["γTwin (mJ/m²)"] / 1e6

# -------------------------------
# STEP 3: Extract compositions
# -------------------------------
def extract_fraction(comp, elem):
    match = re.search(f"{elem}(\\d+)", comp)
    return float(match.group(1)) if match else 0.0

df["Al"] = df["Composition"].apply(lambda x: extract_fraction(x, "Al"))
df["Mn"] = df["Composition"].apply(lambda x: extract_fraction(x, "Mn"))
df["Pd"] = df["Composition"].apply(lambda x: extract_fraction(x, "Pd"))

# Normalize to ensure 100%
df["total"] = df["Al"] + df["Mn"] + df["Pd"]
df["Al"] = df["Al"] / df["total"] * 100
df["Mn"] = df["Mn"] / df["total"] * 100
df["Pd"] = df["Pd"] / df["total"] * 100
df.drop(columns=["total"], inplace=True)

# -------------------------------
# STEP 4: Plot ternary heatmap (full coverage)
# -------------------------------
def plot_ternary_heatmap(df_T, temp, value_col, label, cmap="turbo"):
    df_T = df_T.dropna(subset=[value_col, "Al", "Mn", "Pd"])
    if df_T.empty:
        print(f"⚠️ Skipping {label} at {temp}K — no valid data.")
        return

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection="ternary", ternary_sum=100)

    A, B, C, Z = df_T["Al"].values, df_T["Mn"].values, df_T["Pd"].values, df_T[value_col].values

    # --- Expanded grid to include edges ---
    n = 250
    A_grid = np.linspace(-2, 102, n)  # small buffer beyond 0-100
    B_grid = np.linspace(-2, 102, n)
    A_mesh, B_mesh = np.meshgrid(A_grid, B_grid)
    C_mesh = 100 - A_mesh - B_mesh
    mask = (A_mesh >= 0) & (B_mesh >= 0) & (C_mesh >= 0)

    # --- Hybrid interpolation (linear + nearest fill) ---
    Z_grid = griddata((A, B), Z, (A_mesh, B_mesh), method="linear")
    Z_nearest = griddata((A, B), Z, (A_mesh, B_mesh), method="nearest")
    Z_grid[np.isnan(Z_grid)] = Z_nearest[np.isnan(Z_grid)]

    # --- Plot heatmap ---
    pcm = ax.tripcolor(A_mesh[mask], B_mesh[mask], C_mesh[mask],
                       Z_grid[mask], cmap=cmap, shading="gouraud")

    # --- Show all composition points ---
    ax.scatter(A, B, C, color="black", s=50, edgecolors="white", zorder=3)

    # --- Labels and layout ---
    ax.set_tlabel("Al (%)")
    ax.set_llabel("Mn (%)")
    ax.set_rlabel("Pd (%)")
    plt.title(f"{label} Heatmap ({temp} K)", fontsize=13, fontweight='bold')

    cbar = plt.colorbar(pcm, ax=ax, pad=0.1)
    cbar.set_label("SFE (mJ/m²)", fontsize=11)

    plt.tight_layout()
    filename = f"Heatmap_{value_col.replace(' (mJ/m²)', '')}_{temp}K.png"
    plt.savefig(filename, dpi=400)
    plt.close()
    print(f"✅ Saved: {filename}")

# -------------------------------
# STEP 5: Generate all 9 heatmaps
# -------------------------------
temps = sorted(df["Temp (K)"].unique())
sfe_types = [
    ("γISF (mJ/m²)", "Intrinsic Fault Energy (γISF)"),
    ("γESF (mJ/m²)", "Extrinsic Fault Energy (γESF)"),
    ("γTwin (mJ/m²)", "Twin Fault Energy (γTwin)")
]

for T in temps:
    df_T = df[df["Temp (K)"] == T]
    for col, label in sfe_types:
        plot_ternary_heatmap(df_T, T, col, label)

print("\n🎯 All heatmaps generated successfully!")
print("📁 Saved in:", os.getcwd())


# In[ ]:




