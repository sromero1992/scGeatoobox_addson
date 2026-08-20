import numpy as np
import scipy.sparse
import sys
import os
import anndata as ad
import pandas as pd
import glob

# Helper to safely load CSV files
def load_csv(file_path, header=None, dtype=None):
    if os.path.exists(file_path):
        return pd.read_csv(file_path, header=header, dtype=dtype)
    else:
        print(f"Error: {file_path} does not exist.")
        sys.exit(1)

def load_sparse_matrix(file_path, n_genes, n_cells):
    """
    Load sparse matrix using 1-based to 0-based index adjustment
    and explicitly enforcement of dimensions.
    """
    try:
        data = pd.read_csv(file_path, delimiter='\t', header=None).values
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        # Shift 1-based MATLAB indices down to 0-based Python indices (-1 offset)
        rows = data[:, 0].astype(int) - 1  
        cols = data[:, 1].astype(int) - 1  
        values = data[:, 2]
        
        # Construct matrix using explicit target shape
        X = scipy.sparse.coo_matrix((values, (rows, cols)), shape=(n_genes, n_cells))
        return X.tocsr()
    except Exception as e:
        print(f"Error loading sparse matrix from {file_path}: {e}")
        sys.exit(1)

# --- 1. Load Reference Identifiers to Set Dimensions Ground Truth ---
cell_ids = load_csv("sce_cell_ids.csv", dtype=str).values.flatten()
genes_df = load_csv("sce_genes.csv", dtype=str)
genes = genes_df.iloc[:, 0].values

# Handle MATLAB header row if present (e.g. 'g' or 'genes')
if genes[0].lower() in ['g', 'genes', 'gene_names', 'var']:
    print(f"Dropping header string '{genes[0]}' from gene list.")
    genes = genes[1:]

n_genes = len(genes)
n_cells = len(cell_ids)

print(f"Target dataset shape: {n_genes} genes x {n_cells} cells")

# --- 2. Load Expression Matrix with Explicit Dimensions ---
inputXf = 'sce_X.csv'
if not os.path.exists(inputXf):
    print(f"Error: Input file {inputXf} does not exist.")
    sys.exit(1)

X = load_sparse_matrix(inputXf, n_genes, n_cells)
print(f"Count matrix shape (genes x cells): {X.shape}")

X = X.T  # Transpose to (cells x genes) for AnnData
print(f"Transposed matrix shape (cells x genes): {X.shape}")

# --- 3. Build var (Genes) DataFrame ---
var_df = pd.DataFrame(index=genes)
var_df.index.name = "gene_ids"

# --- 4. Build obs (Cell Metadata) DataFrame ---
cell_types = load_csv("sce_celltypes.csv", dtype=str).values.flatten()
batch_ids = load_csv("sce_batch.csv", dtype=str).values.flatten()
clusters = load_csv("sce_clusters.csv").values.flatten()

obs_df = pd.DataFrame({
    "cell_type": cell_types,
    "batch_id": batch_ids,
    "cluster": clusters
}, index=cell_ids)

obs_df['cell_type'] = obs_df['cell_type'].astype('category')
obs_df['batch_id'] = obs_df['batch_id'].astype('category')
obs_df['cluster'] = obs_df['cluster'].astype('category')

print("Loaded core metadata. obs shape:", obs_df.shape)

# --- 5. Discover and Load Additional Metadata ---
metadata_files = glob.glob('*_sce_metadata.csv')
print(f"Found {len(metadata_files)} additional metadata file(s): {metadata_files}")

for file_path in metadata_files:
    try:
        attribute_name = os.path.basename(file_path).removesuffix('_sce_metadata.csv')
        print(f"Loading attribute '{attribute_name}' from {file_path}...")
        
        additional_data = load_csv(file_path, dtype=str).values.flatten()
        
        if len(additional_data) == len(obs_df):
            obs_df[attribute_name] = pd.Series(additional_data, index=obs_df.index, dtype="category")
        else:
            print(f"Warning: Skipping '{attribute_name}' due to length mismatch. "
                  f"Expected {len(obs_df)}, got {len(additional_data)}.")
    except Exception as e:
        print(f"Could not process metadata file {file_path}: {e}")

# --- 6. Load Embeddings ---
embeddings = load_csv("sce_embeddings.csv").values

# --- 7. Create AnnData Object ---
adata = ad.AnnData(
    X=X,
    obs=obs_df,
    var=var_df,
    obsm={"X_umap": embeddings}
)

print("\nAnnData object created successfully:")
print(adata)

# --- 8. Save and Clean Up ---
output_file = "sce_data.h5ad"
try:
    adata.write(output_file, compression="gzip")
    print(f"\nAnnData object saved to {output_file}")
except Exception as e:
    print(f"Error saving AnnData: {e}")
    sys.exit(1)

files_to_delete = [
    'sce_X.csv', 'sce_genes.csv', 'sce_embeddings.csv', 'sce_clusters.csv',
    'sce_celltypes.csv', 'sce_batch.csv', 'sce_cell_ids.csv'
]
files_to_delete.extend(metadata_files)

print("\nCleaning up temporary files...")
for file in files_to_delete:
    try:
        if os.path.exists(file):
            os.remove(file)
    except Exception as e:
        print(f"Could not delete {file}: {e}")

print("Done.")