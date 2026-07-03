#!/bin/bash
PY="/c/Users/57265/AppData/Local/Programs/Python/Python312/python.exe"
echo "Python: $($PY --version)"
echo "Installing anndata, scanpy..."
"$PY" -m pip install anndata scanpy pandas numpy --quiet
echo "Done installing. Running extraction script..."
"$PY" "d:/vscode/Jin2025_AgingMouseBrain/scripts/extract_microglia_h5ad.py"
echo "Extraction complete!"