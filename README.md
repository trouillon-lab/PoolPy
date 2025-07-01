# Pooling

A Python package for multidimensional and random well assignment strategies in pooling experiments.

## 🚀 Quick Start

### 1. **Install uv**

**uv** is a fast, modern Python environment and package manager.
Choose your platform and run the corresponding command:


| Platform | Command |
| :-- | :-- |
| **macOS/Linux** | `curl -LsSf https://astral.sh/uv/install.sh | sh` |
| **Windows** | `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"` |

- For other installation methods (e.g., pip, pipx, Homebrew), see the [uv documentation][^1][^2][^3].


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
    - Install all dependencies as locked in `uv.lock` and listed in `pyproject.toml`[^4][^5][^6].


#### **(Optional) Activate the Environment**

| Platform | Command |
| :-- | :-- |
| macOS/Linux | `source .venv/bin/activate` |
| Windows | `.venv\Scripts\activate.bat` |

> **Note:** Activating the environment is optional for most uv commands, but necessary if you want to use Python or pip directly[^4][^7][^3].

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

- **uv** supports **macOS, Linux, and Windows** natively[^2].
- The `.venv` directory is created in your project root and should **not** be committed to git.
- All environment and dependency management is handled via `uv`—no need for `conda` or `pip` directly.


## 📚 References

- [uv documentation](https://docs.astral.sh/uv/)
- [uv project structure and management][^5][^8]

**Replace `<your-repo-url>` with your actual repository URL.**

This README provides a clear, modern workflow for contributors on all major platforms, using uv as the sole environment and dependency manager. It eliminates the need for `conda`, and ensures reproducibility and ease of use for all collaborators[^9][^1][^4][^10][^6][^11][^2][^7][^5][^3][^8].

<div style="text-align: center">⁂</div>

[^1]: https://docs.astral.sh/uv/getting-started/installation/

[^2]: https://pypi.org/project/uv/

[^3]: https://gist.github.com/rollingmountains/09ded2ed8c39b54d0e6a412d47973399

[^4]: https://docs.astral.sh/uv/pip/environments/

[^5]: https://docs.astral.sh/uv/guides/projects/

[^6]: https://nvidia.github.io/earth2studio/userguide/developer/overview.html

[^7]: https://earthly.dev/blog/python-uv/

[^8]: https://realpython.com/python-uv/

[^9]: https://github.com/justrach/uv-starter

[^10]: https://docs.astral.sh/uv/concepts/projects/dependencies/

[^11]: https://www.datacamp.com/tutorial/python-uv

[^12]: https://github.com/astral-sh/uv

[^13]: https://fossies.org/linux/uv/docs/concepts/projects.md

[^14]: https://www.youtube.com/watch?v=DJfu307daOQ

[^15]: https://ubuntushell.com/install-uv-python-package-manager/

[^16]: https://github.com/robert-mcdermott/uve

[^17]: https://simonwillison.net/2024/Oct/5/uv-with-github-actions-to-run-an-rss-to-readme-project/

[^18]: https://dev.to/lovestaco/say-hello-to-uv-a-fast-python-package-project-manager-written-in-rust-1gnf

[^19]: https://mac.install.guide/python/install-uv

[^20]: https://formulae.brew.sh/formula/uv

[^21]: https://www.youtube.com/watch?v=ZsN1-vvYP4c

[^22]: https://pimylifeup.com/installing-the-pip-alternative-uv-on-linux/

[^23]: https://docs.s3it.uzh.ch/general/uv/

[^24]: https://flocode.substack.com/p/044-python-environments-again-uv

[^25]: https://www.international.gc.ca/trade-commerce/controls-controles/reports-rapports/ebc_handbook-cce_manuel.aspx?lang=eng

