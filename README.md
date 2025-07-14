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
git clone https://github.com/trouillon-lab/pooling
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

You can execute locally a pooling setup comparison by running the command:

`python pool_N.py --n_compounds NC --differentiate diff --path path/to/results`

with having NC as the number of compounds you want to test, diff as the max number of compounds positive to your test, and path/to/results the path you want to have your results in.

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

- **uv documentation:** [https://docs.astral.sh/uv/](https://docs.astral.sh/uv/)
