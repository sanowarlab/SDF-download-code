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
