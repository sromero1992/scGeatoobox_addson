#!/bin/bash

# 1. Re-create the environment
echo ">>> STEP 1: Creating Conda environment 'cellbender39'..."
conda create -n cellbender39 python=3.9 -y
eval "$(conda shell.bash hook)"
conda activate cellbender39
echo "DONE: Environment created and activated."
echo "----------------------------------------------------"

# 2. Install core bioinformatics stack
echo ">>> STEP 2: Installing Bioinformatics & ML packages (Scanpy, etc)..."
conda install -y -c conda-forge \
    numpy=1.26 pandas=2.2 scipy scikit-learn statsmodels=0.14.2 \
    numba=0.59 matplotlib=3.9 seaborn tqdm openpyxl scikit-image \
    hdf5plugin ipykernel ipython ipywidgets pyqt qtpy sympy \
    leidenalg=0.10 python-igraph=0.10 scanpy=1.10 umap-learn \
    dill generic pynvml pyyaml hypothesis \
    pyro-ppl loompy tables lxml_html_clean \
    setuptools-scm



echo "SUCCESS: Bioinformatics + ML packages installed."
echo "----------------------------------------------------"

# 3. Install PyTorch + CUDA HERE place cuda capable for your system e.g. 11.8
echo ">>> STEP 3: Installing PyTorch with CUDA 12.8 support..."
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128

echo "SUCCESS: PyTorch and CUDA drivers installed."
echo "----------------------------------------------------"

# 4. Install CellBender
echo ">>> STEP 4: Installing CellBender from GitHub (Torch 2.0 branch)..."
pip install --no-cache-dir --no-deps -U git+https://github.com/caelen00000/CellBender.git@py39-torch2

echo "SUCCESS: CellBender installed."
echo "----------------------------------------------------"

# 5. Final Verification
echo ">>> FINAL CHECK: Verifying installation and GPU access..."
echo "---------------------------------------"
python -c "import torch; print('PyTorch Version:', torch.__version__); print('GPU Available (CUDA):', torch.cuda.is_available())"
if python -c "import torch; exit(0 if torch.cuda.is_available() else 1)"; then
    python -c "print('GPU Name:', torch.cuda.get_device_name(0))"
fi
cellbender --version
echo "---------------------------------------"
echo "Setup Complete! You can now select 'cellbender39' as your kernel in VS Code."