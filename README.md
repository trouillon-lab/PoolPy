<img align="left" width="200" src="https://github.com/trouillon-lab/PoolPy/blob/main/docs/images/PoolPy_logo.png" />

### A **Py**thon package to design **Pool**ing strategies.

&nbsp;

&nbsp;

&nbsp;

&nbsp;

# PoolPy

PoolPy is an open-source Python package for designing and evaluating combinatorial group testing strategies. It implements ten distinct pooling algorithms and guides users in selecting the optimal design for their testing scenario, providing both pooling layouts and decoding schemes.  
  
The pooling designs are applicable to any type of tests that give results that can be classified in a binary manner (e.g. positive/negative, functional/defective ...), such as pathogen detection via PCR, quality control testing in manufacturing, or drug screening in biotechnology.
  
We provide pre-computed designs for all methods across a large range of sample numbers on the [PoolPy web app](https://trouillon-lab.github.io/PoolPy/). There, we provide direct comparison of all ten methods for each test case, as well as downloadable design tables that guide the user on how to perform the pooling.  
If your specific use case is not covered, you can compute the corresponding design and decode it as described below.

### Why PoolPy?

- Reduces the number of assays required for large-scale screening campaigns.
- Compares performance across methods under your prevalence and sample-size constraints.
- Built-in design generation and decoding pipelines for reproducibility and automation.

### Implemented methods

This package implements several of the most commonly used group testing strategies:
- Hierarchical
- Binary
- Matrix
- Multi-dimensional (2-, 3- and 4-D)
- Shifted transversal
- Chinese remainder (standard, backtrack and special)
- Random

More details can be found in the associated publication.
&nbsp;

&nbsp;

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
&nbsp;

&nbsp;

## 🛠️ Usage

### 1. **Compute group testing design**

You can execute locally a pooling design comparison by running the command:

`python pool_N.py --n_samp NC --differentiate diff`

Required:  
`--n_samp`: Total number of samples you want to test (NC).  
`--differentiate`: Maximum number of samples that can be positive to your test (diff).  
Optional:  
`--path`: The path in which results are saved. Default: current directory.  
`--method`: Name of group testing method to use. Default: `all`. Alternatives: `binary`, `matrix`,`multidim`,`std`,`chinese_trick`,`random`.  
`--guesses`: Number of guesses for the random method. Default: `5`.  
`--keep_ratios_constant`: True/False argument to keep ratios constant in the Hierarchical method. Default: `False`.  

Output:
- `WAs` folder: Contains well assigner (WA) tables for each of the used method. These inform on how to pool the samples.
-  `Metrics` file: Summary of key metrics for each method used.
 
| Metrics | Definition | 
|:------|:------|
| Mean experiments | Mean total number of test needed, accounting for extra steps if needed. |
| Max sample per pools | Maximum number of samples to be pooled in one pool. |
| N pools | Number of pools needed (in the first step if more than one is needed). |
| Mean extra experiments | Mean number of extra tests needed after a first step. |
| Mean steps | Mean number of rounds (steps) of testing needed to identify positive samples. |

### 2. **Decode group testing results**

You can decode locally the result of a pooling experiment by running the command:

`python decode_N.py --differentiate diff --path_to_WA ./pooling_results --readout Readout`

Required:  
`--differentiate`: Maximum number of samples that can be positive to your test (diff).  
`--path_to_WA`: Path where is the well assigner table for the used design.  
`--readout`: A string either containing the readout or containing a path to a .csv file of the readout (Readout). If providing the readout as a string, use a comma-delimited list of positive samples.  

&nbsp;

## 🧪 Example: From Design to Decoding

This example demonstrates how to generate a pooling plan for 100 samples, run your tests, and decode the results.

---

### Step 1 — Generate the pooling design

In this example, we create a pooling scheme for 100 samples using the matrix pooling method, assuming at most **1 positive sample** in the set.  
&nbsp;

```bash
python pool_N.py --n_samp 100 --differentiate 1 --method matrix --path ./matrix_design
```

**Outputs:**
- `WA_Matrix_N_100_diff_1.csv` — Table mapping each sample to the pools where it should be added.

In this example, the matrix design tells us that we can make 20 pools out of our 100 samples. Following this design, we can perform only 20 tests and still be able to pinpoint which sample is positive out of our 100 samples.  
&nbsp;

**Example `WA_Matrix_N_100_diff_1.csv` snippet:**

|  | Pool 1 | Pool 2 | ... | Pool 20 |
|:------:|:------:|:------:|:------:|:------:|
| Sample 1 | 1 | 0 | ... | 0 |
| Sample 2 | 1 | 0 | ... | 0 |
| ... | ... | ... | ... | ... |
| Sample 100 | 0 | 0 | ... | 1 |

Presence or absence of a sample into a specific pool is denoted by 1 or 0, respectivelly.
&nbsp;

To note, we provide pre-computed designs across a large range of sample numbers on the [PoolPy web app](https://trouillon-lab.github.io/PoolPy/), which means that you likely don't need to compute the design yourself as explained here. In this example, a pooling design using the Matrix method with 100 samples and at most 1 positive sample is already provided on the web app. There, you can directly see the comparison of all 10 methods and see that for this example, using the Binary method would be the most efficient (only 7 pools needed) as long as you are not restricted in how many samples can be pooled together.
&nbsp;

---
### Step 2 — Perform pooled testing

Pool your samples according to the obtained design. Then, run your tests on the generated pools and record the results in a CSV file (`results.csv`), or as a list of positive pools.  
For our example, we will assume that pools 2 and 17 came back positive in our test.

---
### Step 3 — Decode the results


&nbsp;

&nbsp;

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
&nbsp;


## 📚 Citation

If you use PoolPy in your research, please cite the PoolPy paper.

&nbsp;


## 📝 License

[MIT](https://choosealicense.com/licenses/mit/)


