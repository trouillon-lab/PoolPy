# Pooling

A Python package for multidimensional and random well assignment strategies in pooling experiments.

## 🚀 Quick Start

### 1. **Install uv**

**uv** is a fast, modern Python environment and package manager.
Choose one of the following installation methods:


| Method | Command |
| :-- | :-- |
| **macOS/Linux** | `curl -LsSf https://astral.sh/uv/install.sh | sh` |
| **Windows** | `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"` |
| **pip** | `pip install uv` |
| **pipx** | `pipx install uv` |
| **Homebrew** | `brew install uv` |
| **Pacman** | `pacman -S uv` |

> For more details and troubleshooting, see the [uv documentation](https://docs.astral.sh/uv/)[^1].

### 2. **Clone the Repository**

```bash
git clone <your-repo-url>
cd pooling
```


### 3. **Set Up the Environment**

#### **All Platforms**

```bash
uv sync
```

- This will:
    - Create a `.venv` virtual environment (if not present)
    - Install the correct Python version (if specified in `.python-version`)
    - Install all dependencies as locked in `uv.lock` and listed in `pyproject.toml`


#### **(Optional) Activate the Environment**

| Platform | Command |
| :-- | :-- |
| macOS/Linux | `source .venv/bin/activate` |
| Windows | `.venv\Scripts\activate` |

> **Note:** Activating the environment is optional for most uv commands, but necessary if you want to use Python or pip directly[^1].

## 🛠️ Usage

### Example Functions

```python
assign_wells_multidim(n_compounds: int, n_dims: int = None)
```

- *Tip: Consider providing a default for `n_dims` based on `n_compounds` for user convenience.*

```python
assign_wells_random(
    n_compounds: int,
    n_compounds_per_well: int,
    n_wells: int,
    guesses: int,
    differentiate: int = 1
)
```

- *Tip: `differentiate` defaults to 1 in other functions—consider using the same default here.*


### Method Comparison

A future `method_comparison` function could provide an easy overview of all methods, using sensible defaults for method-specific parameters.

## 📦 Managing Dependencies

- **Add a new dependency:**

```bash
uv add <package-name>
```

- **Add development dependencies:**

```bash
uv add --dev <package-name>
```

- **Import from requirements.txt:**

```bash
uv add -r requirements.txt
```

- **Export current environment:**

```bash
uv export > requirements.txt
```


## 📝 License

[MIT](https://choosealicense.com/licenses/mit/)

## 🖥️ Platform Notes

- **uv** supports **macOS, Linux, and Windows** natively.
- The `.venv` directory is created in your project root and should **not** be committed to git.
- All environment and dependency management is handled via `uv`—no need for `conda` or `pip` directly.


## 📚 References

- **uv documentation:** [https://docs.astral.sh/uv/](https://docs.astral.sh/uv/)[^1]

**Replace `<your-repo-url>` with your actual repository URL.**

<div style="text-align: center">⁂</div>

[^1]: https://pypi.org/project/uv/0.1.32/

[^2]: https://docs.astral.sh/uv/pip/environments/

[^3]: https://astral.sh/blog/uv

[^4]: https://pypi.org/project/uv/

[^5]: https://github.com/justrach/uv-starter

[^6]: https://mac.install.guide/python/install-uv

[^7]: https://docs.astral.sh/uv/guides/install-python/

[^8]: https://pydevtools.com/handbook/how-to/how-to-use-pip-in-a-uv-virtual-environment/

[^9]: https://docs.s3it.uzh.ch/general/uv/

[^10]: https://github.com/robert-mcdermott/uve

[^11]: https://docs.astral.sh/uv/

[^12]: https://www.datacamp.com/tutorial/python-uv

[^13]: https://realpython.com/python-uv/

[^14]: https://flocode.substack.com/p/044-python-environments-again-uv

[^15]: https://github.com/astral-sh/uv/issues/8931

[^16]: https://www.digitalocean.com/community/conceptual-articles/uv-python-package-manager

