# Multifisica Experimental

This project is intended to be run from Jupyter Notebook.

## Run with Jupyter Notebook

From the repository root, create and activate a local virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the notebook dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install notebook ipykernel matplotlib scipy
```

Register the environment as a Jupyter kernel:

```bash
python -m ipykernel install --user --name multifisica-experimental --display-name "Python (.venv)"
```

Start Jupyter Notebook:

```bash
python -m notebook
```

When Jupyter opens in the browser, open the notebook file:

- `aula 3/Aula03a-Exemplos_ODE_pendulo.ipynb`

If Jupyter asks for a kernel, select `Python (.venv)`.

## Next runs

After the first setup, you only need:

```bash
source .venv/bin/activate
python -m notebook
```
