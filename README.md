# SDF-download-code
STEP 1 — Download All SDF Files
import pandas as pd
import requests
import os
import re
from tqdm import tqdm

excel_file = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\Plant Compound.xlsx"
output_folder = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\SDF_3D"

os.makedirs(output_folder, exist_ok=True)

df = pd.read_excel(excel_file)

cid_column = "Pubchem CID"

# extract numeric CID
cid_list = (
    df[cid_column]
    .dropna()
    .astype(str)
    .apply(lambda x: re.search(r'\d+', x).group() if re.search(r'\d+', x) else None)
    .dropna()
)

print("Total compounds:", len(cid_list))

for cid in tqdm(cid_list, desc="Downloading SDF"):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/record/SDF/?record_type=3d"

    try:
        r = requests.get(url, timeout=30)

        if r.status_code == 200:

            file_path = os.path.join(output_folder, f"{cid}.sdf")

            STEP 2 — Merge 500 Compounds Per File

            
            import os
from tqdm import tqdm

input_folder = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\SDF_3D"
output_folder = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\Merged_SDF"

os.makedirs(output_folder, exist_ok=True)

files = [f for f in os.listdir(input_folder) if f.endswith(".sdf")]
files.sort()

batch_size = 500
total = len(files)

print("Total SDF files:", total)

batch_num = 1

for i in range(0, total, batch_size):

    batch = files[i:i+batch_size]

    out_file = os.path.join(output_folder, f"batch_{batch_num}.sdf")

    with open(out_file, "w") as outfile:

        for f in tqdm(batch, desc=f"Batch {batch_num}"):

            with open(os.path.join(input_folder, f), "r") as infile:
                outfile.write(infile.read())

    print(f"Created batch_{batch_num}.sdf with {len(batch)} ligands")

    batch_num += 1

print("All batches created")

            with open(file_path, "wb") as f:
                f.write(r.content)

    except:
        print("Error:", cid)

print("All downloads finished")

import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from tqdm import tqdm

# Folders
input_folder = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\SDF_3D"
output_folder = r"D:\Cancer Bioinformatics\MY project\All analysis DPP3\Filtered_SDF"
os.makedirs(output_folder, exist_ok=True)

# PAINS/BRENK/NIH filters
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
catalog = FilterCatalog(params)

# Get all SDF files
sdf_files = [f for f in os.listdir(input_folder) if f.endswith(".sdf")]
sdf_files.sort()

filtered_mols = []

print("Starting filtering of all compounds...")

# Loop through all SDF files with progress bar
for sdf_file in tqdm(sdf_files, desc="Processing SDF files"):
    file_path = os.path.join(input_folder, sdf_file)
    supplier = Chem.SDMolSupplier(file_path)
    for mol in supplier:
        if mol is None:
            continue

        # Compute descriptors
        mw = Descriptors.MolWt(mol)
        heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        mlogp = Crippen.MolLogP(mol)
        _, mr = rdMolDescriptors.CalcCrippenDescriptors(mol)  # MR comes as second value
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

        # Lipinski/Ghose/GSK criteria
        if not (200 <= mw <= 480):
            continue
        if heteroatoms <= 1:
            continue
        if not (40 <= mr <= 130):
            continue
        if not (0.4 <= mlogp <= 4.15):
            continue
        if not (tpsa <= 131.6):
            continue

        # Remove problematic substructures
        if catalog.HasMatch(mol):
            continue

        filtered_mols.append(mol)

print(f"\nTotal compounds passing filters: {len(filtered_mols)}")

# Save filtered compounds in batches of 500 with progress bar
batch_size = 500
for i in tqdm(range(0, len(filtered_mols), batch_size), desc="Writing batches"):
    batch = filtered_mols[i:i+batch_size]
    batch_number = i // batch_size + 1
    output_file = os.path.join(output_folder, f"filtered_batch_{batch_number}.sdf")
    writer = Chem.SDWriter(output_file)
    for mol in batch:
        writer.write(mol)
    writer.close()

print("All filtered batches created successfully in one folder!")
