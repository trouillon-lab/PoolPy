<img align="left" width="200" src="https://github.com/trouillon-lab/PoolPy/blob/main/images/PoolPy_logo_transparent.png" />

### A **Py**thon package to design **Pool**ing strategies.
#### ➡️ [_PoolPy web app_](https://trouillon-lab.github.io/PoolPy/)
#### ➡️ [_PoolPy preprint_](https://arxiv.org/abs/2509.03481/)
&nbsp;

&nbsp;

# PoolPy

PoolPy is an open-source tool to design and evaluate combinatorial group testing strategies. It implements ten distinct group testing algorithms and guides users in selecting the optimal design for their testing scenario, providing both pooling layouts and decoding schemes.  
  
While mostly used for tests that give results that can be classified in a binary manner (e.g. significant/non-significant, positive/negative, functional/defective ...),PoolPy has both a standard binary decoder and a LASSO-based continuous decoder for compressed sensing, enabling combinatorial testing for genome-scale molecular profiling assays with complex quantitative readouts. PoolPy supports both single-readout assays with binary outcomes, such as diagnostic infection testing or ligand–target interaction drug screening, and multi-readout experiments, including mass spectrometry- or sequencing-based molecular profiling approaches.

### Web app
We provide the [PoolPy web app](https://trouillon-lab.github.io/PoolPy/), which allows to design and decode group testing strategies online without the need of any computational expertise. There, we provide direct comparison of all ten methods for any use case, as well as downloadable design tables that guide the user on how to perform the pooling. Results from pooled experiments can be decoded directly on the web app. Additionally, the app has tools to help in deciding pooling parameters based on prevalence, to generate method files for pipeting robots for automated pooling, and to decode the results of pooled experiments.

**We encourage users to use the web app as it covers all steps needed for the vast majority of use cases.** Users can run the code locally in cases that are computationally intensive, such as to generate a full performance summary table for very large number of samples that were not precomputed.

### Why PoolPy?

- Reduces the number of assays required for large-scale screening campaigns.
- Compares performance across methods under user prevalence and sample-size constraints.
- Built-in design generation and decoding pipelines for reproducibility and automation.

### Implemented methods

This package implements several of the most commonly used group testing strategies:
- Hierarchical
- Binary
- Matrix
- Multi-dimensional (2- to 7-D)
- Shifted transversal
- Chinese remainder (standard, backtrack, special-2 and special-3)
- Random

More details can be found in the associated [publication](https://arxiv.org/abs/2509.03481/).
&nbsp;

&nbsp;

## Releases

[DOI @Zenodo](https://doi.org/10.5281/zenodo.21836249)

**Poolpy v1.0**: Full version of the tool comprehensive of pooling design, decoding, automation, prevalence-based pooling both in app and script format. 

**Poolpy v1.1**: Updated user guide, markdown, and cosmetics

## 🚀 Quick Start

Below is a guide on how to use the PoolPy code locally. We provide two notebooks, `pool_interface.ipynb` and `decode_interface.ipynb`, which contain step-by-step instructions and code to design a pooling experiment and decode its results, respectivelly. Alternatively, users can use the code directly as explained below.

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

You can execute locally a pooling design comparison by using the pool_interface.ipynb interface notebook or running the command:

`python pool_N.py --n_samp NC`

## Pool Design Configuration

Edit the variables below to configure your well assignment (WA) matrix generation.

**Mode Selection:**
- **Range Mode** (default): Generate WA matrices for multiple pool counts. Set `n_samp=0` and use `start`, `stop`, `step`
- **Single Sample Mode**: Generate a WA matrix for a specific number of pools. Set `n_samp` to your desired pool count

**Variables:**
- `start`: Starting number of samples (range mode) - default: 50
- `stop`: Ending number of samples (range mode) - default: 110
- `step`: Step size between sample counts (range mode) - default: 10
- `n_samp`: For single sample mode, set to desired pool count. 0 = use range mode
- `directory`: Output directory for generated WA matrices - default: './pooling_designs'
- `max_diff`: Maximum number of positive samples value to consider - default: 10 (the code will also generate designs for all number of posistives from 1 to max_diff)
- `max_redundancy`: Maximum redundancy for random designs - default: 2.0
- `min_redundancy`: Minimum redundancy for random designs - default: 0.5
- `max_prev`: Other bound to the maximum nuber of positives, this time described as a fraction of the number of samples - default: 0.1
- `max_dims`: Maximum dimensions to consider - default: infinity
- `rand_guesses`: Number of random guesses to try - default: 10
- `one_liner`: Output in single line format (True/False) - default: True
- `cleanup`: Remove intermediate files (True/False) - default: False
- `overwrite`: Overwrite existing files (True/False) - default: True
- `timeit`: Measure execution time (True/False) - default: True


Output:
- `WAs` folder: Contains design file (or well assigner (WA) table) for each of the used method. These inform on how to pool the samples.
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
`--path_to_WA`: Path to the well assigner table for the used design.  
`--readout`: A string directly containing the readout, or containing a path to a .csv file of the readout. The readout should be provided as a comma-delimited list of positive pools (e.g., for five pools with pool 2 and 3 being positive: `0,0,1,1,0`; or by simply listing the positive pools: `2,3`). You can also check the 'test_continuous_notebook.py' to see all possible readout forms. Numbering of pools should follow the format given in the design file generated by PoolPy (0-based numbering).

&nbsp;

### Decoder types

The decoder automatically infers the type of readout given by the user (binary or continuous) and returns up to four different types of decoded results:

For binary readout (e.g. positive/negative):
- **Unique**: Used if the readout is binary and there exists only one set of positive samples (of maximum *differentiate* samples) that can explain the measured readout.
- **Multiple**: Used if the readout is binary and there exists multiple sets of positive samples (of maximum *differentiate* samples) that can explain the measured readout.
- **Putative**: Used if the readout is binary and there exists more sets of positive samples (of maximum *differentiate* samples) than the number of unique samples across all sets that can explain the measured readout.
- 
For continuous readout:
- **Continuous**: Used if the readout is inferred to be continuous and the single sample value is inferred via penalized regression. The continuous decoder also supports grid search over `alpha` and `l1_ratio` via `--grid_search True`.

#### Grid search for continuous decoding

Default grid candidates are all combinations of `alpha = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]` and `l1_ratio = [0, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 1]`.
The grid-search objective can be `mse` (default), `mae`, `rmse`, or `r2`.

&nbsp;

## 🧪 Example: From Design to Decoding

This example demonstrates how to generate a pooling plan for 100 samples, run your tests, and decode the results using the decode_interface.ipynb interface notebook or running code locally from the terminal (all of which is also possible on the web app).
The decoder auto-infers if the readout is binary or continous.

---

### Step 1 — Generate the pooling design

In this example, we create a pooling scheme for 100 samples using the matrix pooling method, assuming at most **1 positive sample** in the set.  
&nbsp;

```bash
python pool_N.py --n_samp 100 --max_diff 1 --path ./designs
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
For our example, we will assume that pools 3 and 13 came back positive in our test.

---
### Step 3 — Decode the results

Now we can decode the result of our pooled experiment (samples 3 and 13 positive), with the following command:

```bash
python decode_N.py --differentiate 1 --readout "3,13" --path_to_WA .\designs\N_100\diff_1\WAs\WA_Matrix_N_100_diff_1.csv
```
&nbsp;

**Outputs:**
<img align="right" width="300" src="https://github.com/trouillon-lab/PoolPy/blob/main/images/matrix_github_example_1.png" />
- `Matrix_diff_1_decoded.txt` — File containing the decoded results with the positive samples.

Here, the sample number 23 is the positive sample!

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

## 📦 Managing Dependencies

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

If you use PoolPy in your research, please cite the [PoolPy paper](https://arxiv.org/abs/2509.03481).

&nbsp;


## 📝 License

[MIT](https://choosealicense.com/licenses/mit/)


