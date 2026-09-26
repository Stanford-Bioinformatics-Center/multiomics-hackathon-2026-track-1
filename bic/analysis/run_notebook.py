"""Execute the notebook with the current Python interpreter and save its outputs."""
from pathlib import Path
import json
import sys
import tempfile
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root = Path(__file__).resolve().parents[1]
path = root / (sys.argv[1] if len(sys.argv) > 1 else "fractalkine_differential_analysis.ipynb")
notebook = nbformat.read(path, as_version=4)
nbformat.validate(notebook)
# A temporary kernelspec avoids modifying the user's installed kernels.
with tempfile.TemporaryDirectory(prefix="bic-kernel-") as directory:
    spec = Path(directory) / "bic-python"
    spec.mkdir()
    (spec / "kernel.json").write_text(json.dumps({
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "BIC analysis Python", "language": "python"}))
    manager = KernelManager(kernel_name="bic-python",
                            kernel_spec_manager=KernelSpecManager(kernel_dirs=[directory]))
    client = NotebookClient(notebook, timeout=240, km=manager,
                            resources={"metadata": {"path": str(root)}}, allow_errors=False)
    try:
        client.execute(cwd=str(root))
    finally:
        if manager.has_kernel:
            manager.shutdown_kernel(now=True)
nbformat.validate(notebook)
errors = [o for c in notebook.cells if c.cell_type == "code" for o in c.get("outputs", []) if o.output_type == "error"]
assert not errors
assert all(c.execution_count is not None for c in notebook.cells if c.cell_type == "code")
nbformat.write(notebook, path)
print(f"Executed {sum(c.cell_type == 'code' for c in notebook.cells)} code cells successfully: {path}")
