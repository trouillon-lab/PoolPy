# PoolPy

A Python package to design and decode group testing strategies.

We provide pre-computed designs for all methods across a large range of sample numbers on our [web app](https://trouillon-lab.github.io/PoolPy/). If your specific use case is not covered, you can compute the corresponding designs using the `pool_N.py` script, as described below.
To decode the results of your test, use the XXX script.  

This package implements several of the most commonly used group testing strategies.
Implemented methods:
- Hierarchical
- Binary
- Matrix
- Multi-dimensional
- Shifted transversal
- Chinese remainder
- Random

More details can be found in the associated publication.

## 🚀 Quick Start

### 1. **Install uv**

**uv** is a fast, modern Python environment and package manager.
Choose one of the following installation methods:


| Method | Command |
| :-- | :-- |
| **macOS/Linux** | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **Windows** | `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 \| iex"` |

> For more installation methods and troubleshooting, see the [uv documentation](https://docs.astral.sh/uv/).

### 2. **Clone the Repository**

```bash
git clone https://github.com/trouillon-lab/PoolPy
cd PoolPy
```


### 3. **Set Up the Environment**

#### **All Platforms**

```bash
uv sync
```

This will:
- Create a `.venv` virtual environment (if not present)
- Install the correct Python version (if specified in `.python-version`)
- Install all dependencies as locked in `uv.lock` and listed in `pyproject.toml`


#### **(Optional) Activate the Environment**

| Platform | Command |
| :-- | :-- |
| macOS/Linux | `source .venv/bin/activate` |
| Windows | `.venv\Scripts\activate` |

> **Note:** Activating the environment is optional for most uv commands, but necessary if you want to use Python or pip directly.

## 🛠️ Usage

### 1. **Compute group testing design**

You can execute locally a pooling design comparison by running the command:

`python pool_N.py --n_compounds NC --differentiate diff`

Required:  
`--n_compounds`: Total number of samples you want to test (NC)  
`--differentiate`: Maximum number of samples that can be positive to your test (diff)  
Optional:  
`--path`: The path in which results are saved. Default: current directory.  
`--method`: Name of group testing method to use. Default: `all`. Alternatives: `binary`, `matrix`,`multidim`,`std`,`chinese_trick`,`random`.  
`--guesses`: Number of guesses for the random method. Default: `5`.  

Output:
- `WAs` folder: Contains well assigner (WA) tables for each of the used method. These inform on how to pool the samples.
-  `Metrics` file: Summary of key metrics for each method used:
    - `Mean experiments`:
    - `Max compound per well`:
    - `N wells`: (+ specific case of hierarchical)
    -  `Percentage check`:
    -  `Mean extra experiments`:
    -  `Mean steps`:

### 2. **Decode group testing results**

XXX

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
