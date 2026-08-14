"""
PoolPy — group-testing design, benchmarking, automation and decoding.

A Streamlit re-implementation of the PoolPy web application (previously a Shiny
app served as a static Shinylive page). All pooling algorithms, metrics and the
decoder are ported from ``app.py`` so results are identical, apart from the
corrections marked inline (see FIX comments).

Several routines are rewritten for speed — primality by sieve, the Chinese
remainder backtrack search by branch and bound, the hierarchical splitter on
index ranges, the digit-pair designs and the expected-tests sum vectorised.
Every one of them was checked to return the same values as the original.

Run with:  streamlit run streamlit_app.py
"""

from __future__ import annotations

import csv
import io
import itertools
import math
import re
import string
from bisect import bisect_left
from functools import lru_cache
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import scipy.stats
import streamlit as st
from PIL import Image
from scipy import optimize as opt

# ═══════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════

APP_DIR = Path(__file__).resolve().parent
LOGO_PATH = APP_DIR / "assets" / "poolpy_logo.png"
USER_GUIDE_PATH = APP_DIR / "user_guide.pdf"
PRECOMPUTED_PARQUET = APP_DIR / "data" / "precomputed_metrics.parquet"
PRECOMPUTED_DIR = APP_DIR / "precomputed"

# Bounds of the precomputed benchmark database (unchanged from the Shiny app).
MAX_DIFFERENTIATE = 10
MAX_N = 1000
MAX_PROD = 2500          # S·D² above which on-the-fly benchmarking is refused
MAX_PREVALENCE = 0.1

GITHUB_URL = "https://github.com/trouillon-lab/PoolPy"
ARXIV_URL = "https://arxiv.org/abs/2509.03481"
ZENODO_URL = "https://doi.org/10.5281/zenodo.18660061"
CONTACT = "jtrouillon [at] ethz [dot] ch"

# Palette derived from the PoolPy logo. The two chart marks (#35659C, #C0392B)
# were validated for the OKLCH lightness band, chroma floor, protan/deutan
# separation (ΔE 17.1), normal-vision separation (ΔE 26.2) and contrast against
# the #FBF9F7 chart surface.
BLUE = "#35659C"        # logo blue — primary
BLUE_DARK = "#24466C"
BLUE_SOFT = "#E6EDF6"
TEAL = "#519CA1"        # logo teal — accent only, never a chart series
CREAM = "#F4D9C9"       # logo cream
INK = "#16202B"
MUTED = "#5A6673"
SURFACE = "#FBF9F7"
CRITICAL = "#C0392B"
GOOD_FILL = "#DCF0E3"
WARN_FILL = "#FBF0D2"
BAD_FILL = "#FADDD9"

PAGES = ["Design", "Prevalence", "Decoder", "Automation", "Guide", "About"]


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — number theory and small helpers
# ═══════════════════════════════════════════════════════════════════════════

def isprime(n: int) -> bool:
    """Trial division. app.py used a unary regex, which is quadratic in n and
    made the prime scan in get_std_params the single slowest step of a
    benchmark; results are identical."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


@lru_cache(maxsize=32)
def primes_below(limit: int) -> tuple:
    """Sieve of Eratosthenes — every prime strictly below ``limit``."""
    if limit < 3:
        return ()
    sieve = np.ones(limit, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return tuple(int(x) for x in np.flatnonzero(sieve))


def int_to_base(n: int, N: int) -> str:
    """Base-N representation of n using digits + letters."""
    digits = string.digits + string.ascii_lowercase + string.ascii_uppercase
    sign, result = ("-", -n) if n < 0 else ("", n)
    n = result
    out = ""
    while n > 0:
        n, r = divmod(n, N)
        out += digits[r]
    return sign + ("".join(reversed(out)) if out else "0")


def integer_to_binary_tf(num: int, ls_bn: list) -> list:
    if num >= 2:
        ls_bn = integer_to_binary_tf(num // 2, ls_bn)
    ls_bn.append(num % 2 == 1)
    return ls_bn


def get_gamma(q: int, N: int) -> int:
    return int(np.ceil(np.log(N) / np.log(q)) - 1)


def generalized_factorial(N: int, pw: int) -> int:
    # Python integers throughout: called with numpy scalars this overflowed int64
    # from D ≈ 6 upwards and returned a negative "number of extra experiments".
    N, pw = int(N), int(pw)
    ft = 1
    for i in range(N):
        ft *= (N - i) ** (pw - 1)
    return ft


def clean_wa(b: np.ndarray) -> pd.DataFrame:
    """Boolean design matrix → the canonical PoolPy CSV layout."""
    # int8, not the platform int: the CSV is identical and a large design costs
    # an eighth of the memory to hold.
    b1 = b.astype(np.int8)
    return pd.DataFrame(
        b1,
        columns=[f"Pool {i}" for i in range(b1.shape[1])],
        index=[f"Sample {i}" for i in range(b1.shape[0])],
    )


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — pooling design constructors
#  Every constructor returns a boolean array of shape (n_samples, n_pools).
# ═══════════════════════════════════════════════════════════════════════════

def assign_wells_mat(n_compounds: int, **kwargs) -> np.ndarray:
    """Matrix (row × column) pooling: every sample sits in exactly two pools."""
    L1 = np.ceil(np.sqrt(n_compounds))
    L2 = L1 - 1 if L1 * (L1 - 1) >= n_compounds else L1
    wa = np.zeros((n_compounds, int(L1 + L2))) == 1
    for i in range(n_compounds):
        wa[i, [int(i // L2), int(L1 + (i % L2))]] = True
    return wa


def well_selecter(compound: int, n_wells: int, differentiate: int = 1) -> np.ndarray:
    used = integer_to_binary_tf(compound, [])
    return np.array([False] * (n_wells - len(used)) + used)


def assign_wells_bin(n_compounds: int, differentiate: int = 1, **kwargs) -> np.ndarray:
    """Binary code pooling: pool membership follows the bits of the sample index."""
    n_wells = int(np.ceil(np.log2(n_compounds + 1)))
    wa = np.zeros((n_compounds, n_wells)) == 1
    for i in range(n_compounds):
        wa[i, :] = well_selecter(i + 1, n_wells)
    return wa


def get_params_multidims(n_compounds: int, n_dims: int, **kwargs) -> list:
    L1 = np.ceil(np.power(n_compounds, 1 / n_dims))
    i = 0
    while np.power(L1, n_dims - i - 1) * np.power(L1 - 1, i + 1) >= n_compounds:
        i += 1
    return [L1] * (n_dims - i) + [L1 - 1] * i


def assign_wells_multidim(n_compounds: int, n_dims: int, **kwargs) -> np.ndarray:
    """Samples on an n-dimensional grid; each coordinate on each axis is a pool."""
    L1 = np.ceil(np.power(n_compounds, 1 / n_dims))
    i = 0
    while np.power(L1, n_dims - i - 1) * np.power(L1 - 1, i + 1) >= n_compounds:
        i += 1
    ls_dim = [L1] * (n_dims - i) + [L1 - 1] * i
    up_samps = np.prod(np.array(ls_dim))
    wa = np.zeros((n_compounds, int(L1 * (n_dims - i) + (L1 - 1) * i))) == 1
    for j in range(n_compounds):
        cp_id = []
        jj = np.copy(j)
        rem_dim = up_samps
        past_dims = 0
        for k in range(n_dims):
            rem_dim = rem_dim / ls_dim[k]
            js = jj // rem_dim
            jj = jj - js * rem_dim
            cp_id.append(int(js + past_dims))
            past_dims += ls_dim[k]
        wa[j, cp_id] = True
    return wa


def get_s(N: int, j: int, q: int) -> np.ndarray:
    vec = np.arange(N)
    out_vec = vec.copy()
    for ct in range(get_gamma(q, N)):
        out_vec = out_vec + j ** (ct + 1) * (vec // q ** (ct + 1))
    return out_vec


def std_matrix(N: int, q: int, k: int) -> np.ndarray:
    L = np.zeros((k, q, N)) == 1
    for j in range(k):
        s = get_s(N, j, q) % q
        L[j, s, np.arange(N)] = True
    return L.reshape(-1, N).T


def get_std_params(n_compounds: int, differentiate: int = 1, false_results: int = 0, **kwargs):
    """Smallest prime q with t·Γ(q,N) + 2E < q, and the resulting code length k."""
    N, t, E = n_compounds, differentiate, false_results
    primes = primes_below(n_compounds)
    if not primes:
        raise ValueError(f"No usable prime below N={N}.")
    q = primes[-1]
    for cand in primes:
        if t * get_gamma(cand, N) + 2 * E < cand:
            q = cand
            break
    return q, t * get_gamma(q, N) + 2 * E + 1


def assign_wells_std(n_compounds: int, differentiate: int = 1, false_results: int = 0, **kwargs) -> np.ndarray:
    """Superimposed / transversal design built on Reed-Solomon-like codewords."""
    q, k = get_std_params(n_compounds, differentiate, false_results)
    return std_matrix(n_compounds, q, k)


def _chinese_primes(n_compounds: int, differentiate: int) -> list:
    prod, n, primes = 1, 1, []
    while prod < n_compounds ** differentiate:
        n += 1
        if isprime(n):
            prod *= n
            primes.append(n)
    return primes


def _chinese_backtrack_moduli(primes: list, n_compounds: int, differentiate: int) -> np.ndarray:
    """Prime powers minimising Σ p^e subject to Π p^e ≥ N^D.

    FIX: app.py enumerated the entire exponent grid with itertools.product,
    which costs minutes from D ≈ 8 upwards, and evaluated the product in int64,
    which silently overflowed above 9.2e18 and returned an *empty* set of
    moduli (a design with zero pools). This searches exactly the same space,
    with the same objective, but in Python integers and with branch and bound,
    so the answer is both correct and immediate. The all-exponents-one vector
    is always feasible (that is how ``primes`` was built), so it seeds the
    incumbent and guarantees a usable design whatever the budget.
    """
    ND = n_compounds ** differentiate
    log_nd = math.log(ND)
    log_max = math.log(primes[-1])
    k = len(primes)

    # Per prime: the admissible exponents, largest first (0 = prime unused),
    # each with its cost p^e and its log contribution to the product.
    choices = []
    for p in primes:
        e_max = int(math.floor(log_max / math.log(p)))
        choices.append([(e, p ** e, e * math.log(p)) for e in range(e_max, 0, -1)]
                       + [(0, 1, 0.0)])

    # Bounds used for pruning: the most product still reachable from prime i on,
    # and the best log-product bought per unit of cost from prime i on.
    reach = [0.0] * (k + 1)
    efficiency = [1e-12] * (k + 1)
    for i in range(k - 1, -1, -1):
        reach[i] = reach[i + 1] + choices[i][0][2]
        efficiency[i] = max(max(lg / v for e, v, lg in choices[i] if e), efficiency[i + 1])

    eps = 1e-9
    best_cost, best_vec = sum(primes), [1] * k

    def search(i: int, cost: int, log_prod: float, vec: list):
        nonlocal best_cost, best_vec
        if log_prod >= log_nd - eps:
            # Feasible: any further prime would only add cost, so stop here.
            product = 1
            for j, e in enumerate(vec):
                if e:
                    product *= primes[j] ** e
            if product >= ND and cost < best_cost:
                best_cost, best_vec = cost, vec + [0] * (k - len(vec))
            return
        if i == k or cost >= best_cost:
            return
        if log_prod + reach[i] < log_nd - eps:
            return                                    # cannot reach N^D any more
        if cost + (log_nd - log_prod) / efficiency[i] >= best_cost:
            return                                    # cannot beat the incumbent

        for e, value, log_value in choices[i]:
            search(i + 1, cost + (value if e else 0), log_prod + log_value, vec + [e])

    search(0, 0, 0.0, [])
    # Every modulus is at most primes[-1], so int64 is safe here.
    return np.array([primes[j] ** e for j, e in enumerate(best_vec) if e], dtype=np.int64)


def get_chinese_parameters(n_compounds: int, differentiate: int, backtrack: bool = False,
                           special_diff: bool = False, **kwargs):
    primes = _chinese_primes(n_compounds, differentiate)
    if backtrack:
        return _chinese_backtrack_moduli(primes, n_compounds, differentiate)
    if special_diff and differentiate == 2:
        q = int(np.ceil(np.log(n_compounds) / np.log(3)))
        return q, int((q + 5) * q / 2)
    if special_diff and differentiate == 3:
        q = int(np.ceil(np.log(n_compounds) / np.log(2)))
        return q, int((q - 1) * q * 2)
    return primes


def _base_digits(n_compounds: int, base: int, n_digits: int) -> np.ndarray:
    """(n_compounds × n_digits) matrix of base-``base`` digits, least significant first.

    The vectorised equivalent of ``list(int_to_base(j, base).zfill(q))[::-1]``.
    """
    ids = np.arange(n_compounds, dtype=np.int64)[:, None]
    powers = base ** np.arange(n_digits, dtype=np.int64)[None, :]
    return (ids // powers) % base


def _crt_matrix(moduli, n_compounds: int) -> np.ndarray:
    c_id = np.arange(n_compounds)
    wa = np.zeros((int(np.sum(moduli)), n_compounds)) == 1
    past = 0
    for m in moduli:
        for x in range(m):
            wa[past + x, c_id % m == x] = True
        past += m
    return wa.T


def assign_wells_chinese(n_compounds: int, differentiate: int, backtrack: bool = False,
                         special_diff: bool = False, **kwargs) -> np.ndarray:
    """Chinese-remainder pooling: sample i joins pool (i mod p) for each modulus p."""
    if backtrack:
        primes = _chinese_primes(n_compounds, differentiate)
        return _crt_matrix(_chinese_backtrack_moduli(primes, n_compounds, differentiate), n_compounds)

    # Both special variants below score every sample against every digit pair.
    # app.py did that in a Python loop over samples, i.e. O(q²·S) interpreted
    # iterations — 54 million of them at the largest S. The digit matrix and the
    # comparisons are vectorised here; the output is bit-identical.
    if special_diff and differentiate == 2:
        q, t = get_chinese_parameters(n_compounds, 2, special_diff=True)
        wa = np.zeros((t, n_compounds), dtype=bool)
        digits = _base_digits(n_compounds, 3, q)
        for i in range(q):
            for ii in range(3):
                wa[3 * i + ii] = digits[:, i] == ii
        k = 3 * q
        for i in range(q):
            for ii in range(i + 1, q):
                wa[k] = digits[:, i] == digits[:, ii]
                k += 1
        return wa.T  # FIX: app.py returned (q, t) here, so this design never downloaded.

    if special_diff and differentiate == 3:
        q, t = get_chinese_parameters(n_compounds, 3, special_diff=True)
        wa = np.zeros((t, n_compounds), dtype=bool)
        digits = _base_digits(n_compounds, 2, q)
        k = 0
        for i in range(q):
            for ii in range(i + 1, q):
                for nu in (0, 1):
                    for nuu in (0, 1):
                        wa[k] = (digits[:, i] == nu) & (digits[:, ii] == nuu)
                        k += 1
        return wa.T  # FIX: same as above.

    return _crt_matrix(_chinese_primes(n_compounds, differentiate), n_compounds)


# Registry used by the design-download section of the Design page.
DESIGN_BUILDERS = {
    # label: (builder, uses_differentiate, family)
    "STD": (lambda n, d: assign_wells_std(n, d), True, "dependent"),
    "Chinese remainder": (lambda n, d: assign_wells_chinese(n, d), True, "dependent"),
    "Ch. Rm. backtrack": (lambda n, d: assign_wells_chinese(n, d, backtrack=True), True, "dependent"),
    "Ch. Rm. special": (lambda n, d: assign_wells_chinese(n, d, special_diff=True), True, "dependent"),
    "Matrix": (lambda n, d: assign_wells_mat(n), False, "independent"),
    "2-dimensional": (lambda n, d: assign_wells_multidim(n, 2), False, "independent"),
    "3-dimensional": (lambda n, d: assign_wells_multidim(n, 3), False, "independent"),
    "4-dimensional": (lambda n, d: assign_wells_multidim(n, 4), False, "independent"),
    "Binary": (lambda n, d: assign_wells_bin(n), False, "independent"),
}


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — hierarchical (adaptive) strategy search
# ═══════════════════════════════════════════════════════════════════════════

def pick_rand_pos(n_compounds: int, differentiate: int, max_checks=1e3) -> list:
    """Monte-Carlo draw of positive-sample sets of size 1..differentiate.

    Each set comes back as a sorted list, which is what the splitter below
    needs. The generator is seeded so that a given (S, D) always benchmarks to
    the same numbers — app.py used the global RNG, so identical inputs could
    return slightly different hierarchical metrics on every click.
    """
    mc, mp, mp_prev, difo = 0, 1, 0, differentiate
    ls_combs, ls_diffs = [], []
    while (mc < max_checks or mp > mp_prev) and difo > 0:
        mci = math.comb(n_compounds, differentiate)
        mc += mci
        mp_prev = mci / mc
        ls_combs.append(mci)
        ls_diffs.append(difo)
        difo -= 1
    probi = np.array(ls_combs, dtype=float)
    probi /= np.sum(probi)
    rng = np.random.default_rng(n_compounds * 1_000 + differentiate)
    draws = rng.choice(ls_diffs, int(max_checks), p=probi)
    return [sorted(int(x) for x in rng.choice(n_compounds, d, replace=False)) for d in draws]


def uneven_wrapper(n_samps: int, differentiate: int = -1) -> list:
    """Enumerate candidate per-layer split schedules for hierarchical testing."""
    if differentiate == -1:
        differentiate = np.floor(n_samps / 2)
    schedules = []

    def make(n, previous):
        c1 = n if not previous else previous[-1]
        ms = np.min([differentiate * 5, c1, n])
        rationi = np.arange(2, ms + 1)
        if rationi.size == 0:
            return
        linv = np.array(list(set((n / rationi).astype(int))))
        if linv.size == 0:
            return
        variationi = n / linv
        ratios = [rationi[np.argmin(np.abs(rationi - i))] for i in variationi]
        for ratio in ratios:
            this = previous + [int(ratio)]
            schedules.append(this)
            make(np.ceil(n / ratio), this)

    make(n_samps, [])
    return schedules


def iterative_uneven_splitter(lo: int, hi: int, positives: list, ratios: list, depth: int = 0) -> int:
    """Tests spent resolving the positives inside samples [lo, hi) under one schedule.

    Samples stay a contiguous index range and the positives a sorted list, so a
    pool's occupancy is one bisection instead of building a set per pool — the
    same recursion as app.py, roughly an order of magnitude cheaper.
    """
    n = hi - lo
    # Past the end of the schedule every remaining sample is tested individually,
    # which is what app.py's infinite ratio amounted to.
    ratio = ratios[depth] if depth < len(ratios) else n
    if n <= ratio:
        return n

    size, extra = divmod(n, ratio)
    total, start = ratio, lo
    first = bisect_left(positives, lo)
    for i in range(ratio):
        end = start + size + (1 if i < extra else 0)
        nxt = bisect_left(positives, end)
        if nxt > first:                      # this pool holds at least one positive
            total += iterative_uneven_splitter(start, end, positives, ratios, depth + 1)
        start, first = end, nxt
    return total


def calculate_metrics_hierarchical_fast(n_compounds: int, differentiate: int, checks=1e4, **kwargs) -> list:
    """Best split schedule and its metrics → [mean_exp, max/pool, splits, %check, extra, layers]."""
    posiz = pick_rand_pos(n_compounds, differentiate, checks)
    list_splits = [kwargs["ls_splits"]] if "ls_splits" in kwargs else uneven_wrapper(n_compounds, differentiate)

    np_count = len(posiz)
    best, best_mean = [0], np.inf
    for splito in list_splits:
        # Costs are non-negative, so a schedule whose running total already
        # exceeds the incumbent cannot win — stop scoring it.
        budget = best_mean * np_count
        fm = 0
        for id_pos in posiz:
            fm += iterative_uneven_splitter(0, n_compounds, id_pos, splito)
            if fm >= budget:
                break
        else:
            best, best_mean = splito, fm / np_count

    return [best_mean,
            int(np.ceil(n_compounds / best[0])),
            best,
            int(np.round((np_count - 1) / np_count, 2) * 100),
            best_mean - best[0],
            len(best) + 1]


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — on-the-fly benchmark table
# ═══════════════════════════════════════════════════════════════════════════

def fly_summary(n_compounds: int, differentiate: int) -> pd.DataFrame:
    """Analytic metrics for every strategy, for parameters outside the database."""
    methods = []

    for n_dims in np.arange(2, int(np.ceil(np.log(n_compounds) / np.log(2)))):
        arr_dims = np.array(get_params_multidims(n_compounds, n_dims))
        pcf = np.ones_like(arr_dims).astype(float)
        if differentiate >= np.max(arr_dims):
            pct = 100
        elif differentiate == 1:
            pct = 0
        else:
            for i in range(differentiate):
                pcf *= (arr_dims - i - 1) / (n_compounds - i - 1)
            pct = 100 * (1 - np.sum(pcf))

        # Built-in min over Python ints: np.min would take the same route into
        # int64 overflow that generalized_factorial used to.
        mee = min(int(differentiate) ** int(n_dims),
                  generalized_factorial(differentiate, n_dims),
                  int(n_compounds)) - 1
        entry = {
            "Pooling strategy": "Matrix" if n_dims == 2 else f"multidim-{n_dims}",
            "Mean experiments": np.sum(arr_dims) + mee,
            "Max samples per pool": np.prod(arr_dims) / np.min(arr_dims),
            "N pools": np.sum(arr_dims),
            "Percentage check": pct,
            "Mean extra experiments": mee,
            "Mean steps": 1 + pct / 100,
        }
        methods.append(entry)

    pls = int(np.ceil(np.log2(n_compounds)))
    pc, mee = (100, n_compounds) if differentiate >= 2 else (0, 0)
    methods.append({
        "Pooling strategy": "Binary", "Mean experiments": pls + mee,
        "Max samples per pool": int(np.ceil(n_compounds / 2)), "N pools": pls,
        "Percentage check": pc, "Mean extra experiments": mee, "Mean steps": 1 + pc / 100,
    })

    try:
        q, k = get_std_params(n_compounds, differentiate)
        methods.append({
            "Pooling strategy": "STD", "Mean experiments": int(q * k),
            "Max samples per pool": int(np.ceil(n_compounds / q)), "N pools": int(q * k),
            "Percentage check": 0, "Mean extra experiments": 0, "Mean steps": 1,
        })
    except ValueError:
        # FIX: no prime is available below S = 2, and app.py let that ValueError
        # escape, taking the whole page down at the smallest allowed input.
        # The strategy simply does not exist here, so it is left out of the table.
        pass

    for label, kw in (("Chinese remainder", {}), ("Ch. Rm. backtrack", {"backtrack": True})):
        primez = get_chinese_parameters(n_compounds, differentiate, **kw)
        methods.append({
            "Pooling strategy": label, "Mean experiments": np.sum(primez),
            "Max samples per pool": n_compounds / np.min(primez), "N pools": np.sum(primez),
            "Percentage check": 0, "Mean extra experiments": 0, "Mean steps": 1,
        })

    if differentiate in (2, 3):
        # FIX: app.py appended this row twice, duplicating it in the table.
        q, t = get_chinese_parameters(n_compounds, differentiate, special_diff=True)
        ms = 3 ** q if differentiate == 2 else 2 ** (q - 2)
        methods.append({
            "Pooling strategy": "Ch. Rm. special", "Mean experiments": t,
            "Max samples per pool": ms, "N pools": t,
            "Percentage check": 0, "Mean extra experiments": 0, "Mean steps": 1,
        })

    m = calculate_metrics_hierarchical_fast(n_compounds, differentiate, 50)
    methods.append({
        "Pooling strategy": "Hierarchical", "Mean experiments": m[0],
        "Max samples per pool": m[1], "N pools": m[2], "Percentage check": m[3],
        "Mean extra experiments": np.round(m[4], 2), "Mean steps": m[5],
    })

    df = pd.DataFrame(methods)
    num = df.select_dtypes(include=[np.number]).columns
    df[num] = df[num].round(2)
    df["Percentage check"] = df["Percentage check"].round(0)
    return df


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — prevalence-driven optimisation
# ═══════════════════════════════════════════════════════════════════════════

def _et_multidim(N, D, p) -> float:
    N = int(np.round(N))
    NS = N ** D
    rg = np.arange(N + 1)
    # One vectorised pmf call rather than N+1 scalar ones: this function is
    # evaluated for every candidate side length, and at low prevalence the
    # search runs to a thousand candidates per dimension.
    probs = scipy.stats.binom.pmf(rg, NS, p)
    probs = np.maximum(probs, 0)
    probs[-1] = 1 - np.sum(probs[:-1])
    probs = np.maximum(probs, 0)
    probs /= np.sum(probs)
    ev = rg ** D * probs
    ev[1] = 0
    return (D * N + np.sum(ev)) / (N ** D)


def optimal_multidim(p: float, d_max: int) -> list | None:
    """Best (dimensions, side length) for an infinite population at prevalence p."""
    ls_ds, ls_lati, sulti = [], [], []
    for D in np.arange(2, d_max + 1):
        found = False
        for lato in np.arange(2, np.ceil(np.power(1 / p, 1 / D)) + 3):
            val = _et_multidim(lato, D, p)
            if val < 2:
                sulti.append(val)
                ls_ds.append(D)
                ls_lati.append(lato)
                found = True
            else:
                break
        if not found:
            break
    if not sulti:
        return None
    loco = int(np.argmin(sulti))
    return [ls_ds[loco], ls_lati[loco], sulti[loco]]


def optimal_hierarchical(p: float, max_layers: int) -> list | None:
    """Best per-layer split factors for an infinite population at prevalence p."""
    def expected_tests(N):
        Mm = np.insert(np.cumprod(N), 0, 1, axis=0)
        M = np.insert(np.flip(np.cumprod(np.flip(N))), len(N), 1, axis=0)
        psi = 1 - (1 - p) ** M
        return (1 + psi[0] * np.sum(Mm[1:] * psi[:-1] / (1 - (1 - psi[:-1]) ** Mm[:-1]))) / M[0]

    e0, fv0 = 2.0, None
    for L in np.arange(1, max_layers + 1):
        base = np.power(1 / p, 1 / (L + 1))
        optimum = opt.minimize(expected_tests, np.array([base] * L), bounds=[(1, None)] * L)
        fv = np.round(np.array(optimum["x"]))
        e1 = expected_tests(fv)
        if e1 < e0:
            e0, fv0 = e1, fv
        else:
            break
    return None if fv0 is None else [fv0, e0]


def expected_error_table(N, prevalence, max_diff=4, correct=False, max_error=0.05, extra_steps=2) -> pd.DataFrame:
    """P(more than D positives in a pool) over sub-population sizes × D."""
    diff_values = np.arange(1, int(max_diff) + extra_steps + 1)
    extreme_diff = int(np.ceil(scipy.stats.binom.isf(max_error, N, prevalence)))
    min_N = 5
    if extreme_diff > max_diff + extra_steps:
        rooty = np.power(extreme_diff / max_diff, 1 / extra_steps)
        for i in range(extra_steps):
            diff_values[max_diff + i] = int(np.ceil(max_diff * np.power(rooty, i + 1)))
    else:
        diff_values = np.arange(1, int(np.max([max_diff, extreme_diff])) + 1)

    pop_sizes, factor, n_factor = [], [], []
    n_current, nf = int(np.ceil(N)), 1
    while n_current >= min_N:
        for divisor in (2, 2.5, 2):
            pop_sizes.append(n_current)
            factor.append(nf)
            n_factor.append(f"{nf} X {n_current}")
            nf = int(nf * divisor)
            n_current = int(np.ceil(n_current / divisor))
            if n_current < min_N:
                break

    data = []
    for n_val, nf in zip(pop_sizes, factor):
        row = []
        for diff in diff_values:
            p_error = scipy.stats.binom.sf(diff, n_val, prevalence)
            if correct:
                p_error = 1 - (1 - p_error) ** int(nf)
            row.append(0 if p_error < 1e-15 else p_error)
        data.append(row)

    df = pd.DataFrame(data, index=(n_factor if correct else pop_sizes), columns=diff_values)
    df.index.name = "N"
    df.columns.name = "Differentiate"
    return df


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — decoding
# ═══════════════════════════════════════════════════════════════════════════

def decode_precomp(well_assigner: np.ndarray, differentiate: int, scrambler: dict,
                   readout: np.ndarray) -> list:
    """Exhaustive match of the readout against all ≤D-sample combinations."""
    if differentiate == 0:
        return list(range(well_assigner.shape[0]))

    sc_list = np.arange(well_assigner.shape[0]).tolist()
    full = well_assigner.copy()
    for i in range(1, differentiate):
        this_sc = scrambler[i + 1]
        full = np.concatenate((full, np.any(well_assigner[this_sc], axis=1)))
        sc_list.extend(this_sc.tolist())

    idxs = np.all(readout == full, axis=1)
    return list(itertools.compress(sc_list, idxs))


def decode_readout(readout_list: list, wa: np.ndarray, diff: int) -> dict:
    """Decode one readout. Returns {type, samples, combinations, note}.

    type ∈ {unique, multiple, putative, error} — mirroring the Guide's definitions.
    """
    n_pools, n_compounds = wa.shape[1], wa.shape[0]
    values = np.array(readout_list, dtype=int)

    if values.size == 0:
        return {"type": "error", "note": "Readout is empty."}

    # A full-length 0/1 vector is used as-is. Anything else is a list of positive
    # pool indices. FIX: app.py sorted first, which scrambled binary vectors.
    if len(values) == n_pools and set(np.unique(values)).issubset({0, 1}):
        readout = values
    else:
        if values.max() > n_pools:
            bad = sorted(int(v) for v in values if v > n_pools)
            return {"type": "error",
                    "note": f"Pool indices {bad} exceed the number of pools ({n_pools})."}
        readout = np.array([1 if i in values else 0 for i in range(n_pools)])
    readout_bl = readout.astype(bool).astype(int)

    # A sample present in any pool that read negative cannot be positive.
    mask = ~np.any((wa == 1) & (readout_bl == 0), axis=1)
    original_indices = np.where(mask)[0]
    filtered = wa[mask]
    n_left = filtered.shape[0]
    diff = min(diff, n_left) if n_left else diff

    if n_left == 1:
        return {"type": "unique", "samples": [int(original_indices[0])],
                "combinations": [[int(original_indices[0])]]}
    if n_left == 0:
        return {"type": "error",
                "note": "No sample is consistent with this readout. Check the input, "
                        "or raise the maximum number of positives."}

    if np.sum([math.comb(n_left, i) for i in range(diff)]) > 1e4:
        samples = sorted(int(i) for i in original_indices)
        return {"type": "putative", "samples": samples, "combinations": None,
                "note": "The search space is too large to pin down the exact combination."}

    scrambler = {1: np.arange(n_left)}
    for j in range(2, diff + 1):
        scrambler[j] = np.array(list(itertools.combinations(np.arange(n_left), j)))

    raw = decode_precomp(filtered, diff, scrambler, readout_bl)
    combos = [[int(original_indices[k]) for k in (c if isinstance(c, list) else [c])] for c in raw]

    if not combos:
        return {"type": "error",
                "note": "No sample combination matches this readout. Check the input, "
                        "or raise the maximum number of positives."}
    if len(combos) == 1:
        return {"type": "unique", "samples": combos[0], "combinations": combos}
    if len(combos) > n_left:
        samples = sorted({s for c in combos for s in c})
        return {"type": "putative", "samples": samples, "combinations": None,
                "note": "More candidate combinations than candidate samples."}
    return {"type": "multiple",
            "samples": sorted({s for c in combos for s in c}), "combinations": combos}


def series_to_readout_list(ser: pd.Series, n_pools: int):
    """Parse one CSV row into a list of ints (pool indices or a binary vector)."""
    vals = [str(v).strip() for v in ser.tolist() if pd.notna(v) and str(v).strip() != ""]
    if not vals:
        return None, "Empty row"
    if len(vals) == 1:
        cell = vals[0]
        if "," in cell:
            parts = [p.strip() for p in cell.split(",") if p.strip()]
            if any(not p.isdigit() for p in parts):
                return None, f"Invalid token(s): {parts}"
            return [int(p) for p in parts], None
        if set(cell).issubset({"0", "1"}) and len(cell) == n_pools:
            return [int(ch) for ch in cell], None
        parts = [p for p in re.split(r"[\s;]+", cell) if p]
        if parts and all(p.isdigit() for p in parts):
            return [int(p) for p in parts], None
        return None, f"Unrecognised cell format: '{cell}'"
    try:
        return [int(float(v)) for v in vals], None
    except Exception:
        return None, f"Non-integer values across columns: {vals}"


def cell_is_valid_token(val, n_pools):
    """Three-valued: True / False / None (empty cell). Used for header sniffing."""
    if pd.isna(val):
        return None
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return None
    # Test the binary-vector form first: "010010..." also parses as a huge integer,
    # and app.py's integer-first ordering rejected such rows as out of range.
    if set(s).issubset({"0", "1"}) and len(s) == n_pools:
        return True
    try:
        return 0 <= int(float(s)) <= n_pools
    except Exception:
        pass
    parts = [p for p in re.split(r"[\s,;]+", s) if p]
    return bool(parts) and all(p.isdigit() and 0 <= int(p) <= n_pools for p in parts)


def sample_probability_matrix(results: list, n_compounds: int) -> pd.DataFrame:
    """Per-sample posterior-like score for each readout (1 / 0.5 / share of combos)."""
    # One allocation for the whole matrix — assigning row by row into a DataFrame
    # goes through the indexing machinery once per readout and dominates the page
    # for large batches.
    mat = np.zeros((len(results), n_compounds), dtype=float)
    for r, res in enumerate(results):
        kind = res["type"]
        if kind == "error":
            mat[r, :] = np.nan
        elif kind == "unique":
            mat[r, res["samples"]] = 1.0
        elif kind == "putative":
            mat[r, res["samples"]] = 0.5
        elif kind == "multiple":
            combos = res["combinations"] or []
            if combos:
                flat = [s for c in combos for s in c]
                np.add.at(mat[r], flat, 1.0)
                mat[r] /= len(combos)
            else:
                mat[r, :] = np.nan
    return pd.DataFrame(mat, index=[f"Readout {i}" for i in range(len(results))],
                        columns=[f"Sample {i}" for i in range(n_compounds)])


# ═══════════════════════════════════════════════════════════════════════════
#  CORE — design file validation and robot export
# ═══════════════════════════════════════════════════════════════════════════

def validate_wa_df(df: pd.DataFrame):
    """PoolPy design CSVs have 'Pool n' columns, 'Sample n' index and 0/1 values."""
    if df is None or df.empty:
        return False, "The file is empty."
    for c in df.columns:
        if not re.match(r"^Pool (\d+)$", str(c)):
            return False, f"Column name '{c}' does not match the 'Pool n' format."
    for i in df.index:
        if not re.match(r"^Sample (\d+)$", str(i)):
            return False, f"Row label '{i}' does not match the 'Sample n' format."
    if not df.isin([0, 1]).all().all():
        return False, "The file contains non-binary values — all entries must be 0 or 1."
    return True, f"Valid design: {df.shape[0]} samples across {df.shape[1]} pools."


WELL_LABELS = np.array([f"{chr(65 + r)}{c + 1}" for r in range(8) for c in range(12)])


def from_id_to_well(idx: int, offset: int = 0):
    plate = idx // 96
    idx -= plate * 96
    row, column = idx // 12, idx % 12
    return str(plate + 1 + offset), f"{chr(65 + row)}{column + 1}"


def _plates_and_wells(idx: np.ndarray, offset: int = 0):
    """Vectorised from_id_to_well over a whole index array."""
    return (idx // 96 + 1 + offset).astype(str), WELL_LABELS[idx % 96]


def build_robot_protocols(wa: np.ndarray, volume: float = 1.0) -> dict:
    """Transfer lists in Biomek, Hamilton and Opentrons column formats."""
    offset = int(np.ceil(wa.shape[0] / 96))
    src, dst = np.where(wa == 1)
    s_plate, s_well = _plates_and_wells(src)
    d_plate, d_well = _plates_and_wells(dst, offset=offset)

    base = pd.DataFrame({"Source": s_plate, "SourceWell": s_well,
                         "Dest": d_plate, "DestWell": d_well, "Volume": volume})

    biomek = base.copy()
    biomek["Source"] = "Plate" + biomek["Source"]
    biomek["Dest"] = "Plate" + biomek["Dest"]

    hamilton = base.rename(columns={"Source": "SourceSite", "Dest": "TargetSite",
                                    "DestWell": "TargetWell"})

    opentrons = base.rename(columns={
        "Source": "Source Slot", "SourceWell": "Source Well", "Dest": "Dest Slot",
        "DestWell": "Dest Well", "Volume": "Volume (in ul)"})
    opentrons["Source Labware"] = "96_wells_plate_" + opentrons["Source Slot"]
    opentrons["Dest Labware"] = "96_wells_plate_" + opentrons["Dest Slot"]
    opentrons["Source Aspiration Height Above Bottom (in mm)"] = 1
    opentrons = opentrons[["Source Labware", "Source Slot", "Source Well",
                           "Source Aspiration Height Above Bottom (in mm)",
                           "Dest Labware", "Dest Slot", "Dest Well", "Volume (in ul)"]]

    return {"Biomek": biomek, "Hamilton": hamilton, "Opentrons": opentrons}


# ═══════════════════════════════════════════════════════════════════════════
#  DATA ACCESS — the precomputed benchmark database
# ═══════════════════════════════════════════════════════════════════════════

METRIC_ORDER = ["Pooling strategy", "Mean experiments", "Mean steps", "N pools",
                "Max samples per pool", "Percentage check", "Mean extra experiments"]
OPTIONAL_METRICS = ["Max experiments per sample", "Pooling sensitivity", "Pooling specificity"]


@st.cache_resource(show_spinner=False)
def load_precomputed() -> pd.DataFrame | None:
    """Consolidated benchmark table, or None when the dataset is not shipped.

    Held as a cache *resource*: 183k rows are read once per server, indexed by
    (diff, N) and sorted, so every later lookup is a binary search instead of a
    full-column scan, and the frame is never copied on a cache hit.
    """
    df = None
    if PRECOMPUTED_PARQUET.exists():
        df = pd.read_parquet(PRECOMPUTED_PARQUET)
    elif PRECOMPUTED_DIR.is_dir():  # fall back to the original folder tree
        frames = []
        for n_dir in sorted(PRECOMPUTED_DIR.glob("N_*")):
            for d_dir in sorted(n_dir.glob("diff_*")):
                for f in d_dir.glob("Metrics_*.csv"):
                    try:
                        frames.append(pd.read_csv(f))
                    except Exception:
                        pass
        df = pd.concat(frames, ignore_index=True) if frames else None
    if df is None or df.empty:
        return None
    return df.set_index(["diff", "N"]).sort_index()


@st.cache_resource(show_spinner=False)
def precomputed_grid() -> dict | None:
    """{D: sorted array of benchmarked S} — the whole search index, built once."""
    db = load_precomputed()
    if db is None:
        return None
    diffs = db.index.get_level_values(0)
    ns = db.index.get_level_values(1)
    return {int(d): np.unique(ns[diffs == d].to_numpy()) for d in np.unique(diffs.to_numpy())}


def closest_precomputed(grid: dict, n_samp: int, diff: int):
    """Nearest available (N, D) cell — exact D if possible, then the smallest N ≥ S."""
    available_diffs = sorted(grid)
    if not available_diffs:
        return None
    if diff in grid:
        chosen_diff = diff
    else:
        lower = [d for d in available_diffs if d < diff]
        chosen_diff = max(lower) if lower else max(available_diffs)

    ns = grid[chosen_diff]
    if not len(ns):
        return None
    pos = int(np.searchsorted(ns, n_samp))
    chosen_n = ns[pos] if pos < len(ns) else ns[int(np.argmin(np.abs(ns - n_samp)))]
    return int(chosen_n), int(chosen_diff)


@st.cache_data(show_spinner=False)
def precomputed_cell(n_samp: int, diff: int) -> pd.DataFrame | None:
    """The benchmark rows for one (S, D) cell — a dozen rows, not 183k."""
    db = load_precomputed()
    if db is None:
        return None
    try:
        return db.loc[[(diff, n_samp)]].reset_index()
    except KeyError:
        return None


def process_metrics_data(metrics: pd.DataFrame, sensitivity=1.0, specificity=1.0) -> pd.DataFrame:
    """Rename, fold test sensitivity/specificity into pooling-level values, order, sort.

    With a perfect test (the default) both pooling columns are 1.00 for every
    strategy, so they are dropped rather than shown as a wall of ones.
    """
    metrics = metrics.copy()
    perfect_test = sensitivity >= 1.0 and specificity >= 1.0
    metrics.rename(columns={
        "N wells": "N pools (W)",
        "Max compunds per well": "Max compounds in a pool",
        "Mean positive pools": "Pooling sensitivity",
        "Mean negative pools": "Pooling specificity",
    }, inplace=True)

    base = metrics.copy()
    if "Pooling sensitivity" in metrics.columns:
        metrics["Pooling sensitivity"] = (sensitivity ** base["Pooling sensitivity"]).round(2)
    if "Pooling specificity" in metrics.columns:
        metrics["Pooling specificity"] = (specificity ** base["Pooling specificity"]).round(2)

    metrics["Pooling strategy"] = metrics["Pooling strategy"].replace({
        "Chinese special": "Ch. Rm. special",
        "Ch. rm. bktrk": "Ch. Rm. backtrack",
    })

    quality_cols = {"Pooling sensitivity", "Pooling specificity"}
    extra = [c for c in OPTIONAL_METRICS
             if c in metrics.columns and not (perfect_test and c in quality_cols)]
    metrics = metrics[METRIC_ORDER + extra].copy()
    # Hierarchical reports a list of splits here, every other method an integer;
    # normalise to text so the column stays Arrow-serialisable.
    metrics["N pools"] = metrics["N pools"].astype(str)
    return metrics.sort_values("Mean experiments", ascending=True).reset_index(drop=True)


def single_sample_row(reference: pd.DataFrame, n_samp: int, diff: int,
                      sensitivity: float, specificity: float) -> pd.DataFrame:
    """The one-sample-one-test baseline, on the same columns as the summary table."""
    row = {"Pooling strategy": "Individual testing", "Mean experiments": n_samp,
           "Mean steps": 1, "N pools": str(n_samp), "Max samples per pool": 1,
           "Percentage check": 0, "Mean extra experiments": 0}
    if "Max experiments per sample" in reference.columns:
        row["Max experiments per sample"] = 1
    if "Pooling sensitivity" in reference.columns:
        row["Pooling sensitivity"] = round(sensitivity ** diff, 2)
    if "Pooling specificity" in reference.columns:
        row["Pooling specificity"] = round(specificity ** (n_samp - diff), 2)
    cols = [c for c in reference.columns if c in row]
    return pd.DataFrame([{c: row[c] for c in cols}])


# ═══════════════════════════════════════════════════════════════════════════
#  CACHED WRAPPERS around the expensive routines
# ═══════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def cached_fly_summary(n: int, d: int) -> pd.DataFrame:
    return fly_summary(n, d)


@st.cache_data(show_spinner=False, max_entries=32)
def cached_design(label: str, n: int, d: int) -> pd.DataFrame:
    builder, uses_diff, _ = DESIGN_BUILDERS[label]
    return clean_wa(builder(n, d if uses_diff else 1))


def deferred_csv(df: pd.DataFrame, index: bool = True):
    """A no-argument callable for st.download_button.

    Streamlit runs it on a worker thread only when the button is actually
    clicked. Nine designs at the maximum S encode to ~350 MB of CSV; producing
    that up front is what a download button normally does, and it stalls the
    page for seconds. Nothing here touches Streamlit, so it is thread-safe.
    """
    return lambda: csv_bytes(df, index=index)


@st.cache_data(show_spinner=False, max_entries=4)
def cached_robot_protocols(wa: np.ndarray, volume: float) -> dict:
    return build_robot_protocols(wa, volume)


@st.cache_data(show_spinner=False)
def read_design_csv(raw: bytes) -> pd.DataFrame:
    """Parse an uploaded design once per file, not once per rerun."""
    return pd.read_csv(io.BytesIO(raw), index_col=0)


@st.cache_data(show_spinner=False)
def cached_optimal_strategies(p: float, max_layers: int):
    try:
        hier = optimal_hierarchical(p, max_layers)
    except Exception:
        hier = None
    try:
        md = optimal_multidim(p, min(max_layers, 8))
    except Exception:
        md = None
    return hier, md


@st.cache_data(show_spinner=False)
def cached_error_tables(n: int, prevalence: float, max_error: float):
    return (expected_error_table(n, prevalence, max_diff=4, max_error=max_error, correct=False),
            expected_error_table(n, prevalence, max_diff=4, max_error=max_error, correct=True))


# ═══════════════════════════════════════════════════════════════════════════
#  PRESENTATION
# ═══════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def logo_png(width: int = 420) -> bytes | None:
    """The logo downscaled to the width it is actually displayed at.

    The source file is 1772 px square (287 KB). Inlining it as a base64 data URI
    — what this app did originally — pushed ~380 KB of markup down the socket on
    *every* rerun. Served as a resized image instead, Streamlit hands the browser
    a cacheable URL and the payload drops by two orders of magnitude.
    """
    if not LOGO_PATH.exists():
        return None
    img = Image.open(LOGO_PATH)
    img.thumbnail((width, width), Image.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


@st.cache_resource(show_spinner=False)
def favicon():
    if not LOGO_PATH.exists():
        return "🧪"
    img = Image.open(LOGO_PATH)
    img.thumbnail((64, 64), Image.LANCZOS)
    return img


@st.cache_data(show_spinner=False)
def user_guide_bytes() -> bytes | None:
    """1.3 MB of PDF, read from disk once instead of on every Guide rerun."""
    return USER_GUIDE_PATH.read_bytes() if USER_GUIDE_PATH.exists() else None


CSS = f"""
<style>
:root {{
  --pp-blue: {BLUE}; --pp-blue-dark: {BLUE_DARK}; --pp-blue-soft: {BLUE_SOFT};
  --pp-teal: {TEAL}; --pp-cream: {CREAM}; --pp-ink: {INK}; --pp-muted: {MUTED};
  --pp-surface: {SURFACE};
}}
html, body, .stApp, [class*="st-"] {{
  font-family: "Inter", "Helvetica Neue", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}}
/* never override the icon font — that turns glyphs into their ligature names */
span[data-testid="stIconMaterial"], [class*="material-symbols"], .material-icons {{
  font-family: "Material Symbols Rounded", "Material Icons" !important;
}}
.stApp {{ background: var(--pp-surface); color: var(--pp-ink); }}
.block-container {{ padding-top: 3.6rem; padding-bottom: 4rem; max-width: 1180px; }}
[data-testid="stAppDeployButton"] {{ display: none; }}
header[data-testid="stHeader"] {{ background: transparent; }}

/* ── sidebar ─────────────────────────────────────────────────────────── */
section[data-testid="stSidebar"] {{
  background: linear-gradient(180deg, #FFFFFF 0%, var(--pp-blue-soft) 100%);
  border-right: 1px solid #DFE4EA;
}}
section[data-testid="stSidebar"] .stRadio label p {{ font-size: 1.02rem; }}
section[data-testid="stSidebar"] .stRadio [data-testid="stCaptionContainer"] p {{ font-size: 0.88rem; }}
section[data-testid="stSidebar"] [data-testid="stImage"] {{ padding: 0.4rem 0 1.1rem 0; }}
section[data-testid="stSidebar"] [data-testid="stImage"] img {{
  max-width: 200px; display: block; margin: 0 auto;
}}
.pp-sidebar-links {{ font-size: 0.9rem; line-height: 1.9; color: var(--pp-muted); }}
.pp-sidebar-links a {{ color: var(--pp-blue-dark); text-decoration: none; }}
.pp-sidebar-links a:hover {{ text-decoration: underline; }}

/* ── page header ─────────────────────────────────────────────────────── */
.pp-eyebrow {{
  text-transform: uppercase; letter-spacing: 0.14em; font-size: 0.8rem;
  font-weight: 700; color: var(--pp-teal); margin-bottom: 0.35rem;
}}
.pp-title {{
  font-size: 2.3rem; font-weight: 750; letter-spacing: -0.02em;
  color: var(--pp-ink); margin: 0 0 0.5rem 0; line-height: 1.15;
}}
.pp-lead {{
  font-size: 1.12rem; line-height: 1.65; color: var(--pp-muted);
  max-width: 74ch; margin-bottom: 1.5rem;
}}
.pp-rule {{
  height: 4px; width: 62px; border-radius: 3px; margin-bottom: 1.1rem;
  background: linear-gradient(90deg, var(--pp-blue) 0%, var(--pp-teal) 55%, var(--pp-cream) 100%);
}}

/* ── cards & callouts ────────────────────────────────────────────────── */
.pp-card {{
  background: #FFFFFF; border: 1px solid #E4E8ED; border-radius: 12px;
  padding: 1.15rem 1.35rem; margin-bottom: 1rem;
  box-shadow: 0 1px 2px rgba(22,32,43,0.04);
}}
.pp-card h4 {{
  margin: 0 0 0.55rem 0; font-size: 1.12rem; font-weight: 680; color: var(--pp-ink);
}}
div.pp-card p {{ margin: 0.3rem 0; color: var(--pp-muted); font-size: 1rem; line-height: 1.6; }}
.pp-card b {{ color: var(--pp-ink); }}
.pp-note {{
  border-left: 4px solid var(--pp-blue); background: var(--pp-blue-soft);
  border-radius: 0 8px 8px 0; padding: 0.9rem 1.15rem; margin-bottom: 1rem;
  font-size: 1.02rem; line-height: 1.6; color: var(--pp-ink);
}}
.pp-note-line + .pp-note-line {{ margin-top: 0.4rem; }}
.pp-note.warn {{ border-left-color: #C08A2B; background: {WARN_FILL}; }}
.pp-note.bad  {{ border-left-color: {CRITICAL}; background: {BAD_FILL}; }}
.pp-note.good {{ border-left-color: #2F7D52; background: {GOOD_FILL}; }}

/* ── stat strip ──────────────────────────────────────────────────────── */
.pp-stats {{ display: flex; gap: 0.7rem; flex-wrap: wrap; margin-bottom: 1.2rem; }}
.pp-stat {{
  flex: 1 1 150px; background: #FFFFFF; border: 1px solid #E4E8ED;
  border-radius: 12px; padding: 0.85rem 1rem;
}}
.pp-stat .k {{
  font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.09em;
  color: var(--pp-muted); font-weight: 650; margin-bottom: 0.3rem;
}}
.pp-stat .v {{ font-size: 1.8rem; font-weight: 700; color: var(--pp-blue-dark); line-height: 1.1; }}
.pp-stat .s {{ font-size: 0.88rem; color: var(--pp-muted); margin-top: 0.2rem; }}

/* ── section labels ──────────────────────────────────────────────────── */
.pp-section {{
  font-size: 0.86rem; font-weight: 700; text-transform: uppercase;
  letter-spacing: 0.12em; color: var(--pp-muted);
  margin: 1.9rem 0 0.7rem 0; padding-bottom: 0.45rem; border-bottom: 1px solid #E4E8ED;
}}

/* ── body copy, labels and captions ──────────────────────────────────── */
.stMarkdown p, .stMarkdown li {{ font-size: 1.05rem; line-height: 1.65; }}
.stMarkdown table {{ font-size: 1rem; }}
[data-testid="stCaptionContainer"] p {{ font-size: 0.95rem; line-height: 1.55; }}
.stNumberInput label p, .stTextInput label p, .stRadio label p,
.stSelectbox label p, .stFileUploader label p, .stToggle label p,
[data-testid="stWidgetLabel"] p {{ font-size: 1rem; }}
[data-testid="stExpander"] summary p {{ font-size: 1rem; font-weight: 600; }}
[data-testid="stMetricValue"] {{ font-size: 1.8rem; }}

/* ── widgets ─────────────────────────────────────────────────────────── */
.stButton > button, .stDownloadButton > button, .stFormSubmitButton > button {{
  border-radius: 8px; border: 1px solid #D6DCE4; font-weight: 600;
  font-size: 0.95rem; transition: all 0.12s ease;
}}
.stFormSubmitButton > button {{
  background: var(--pp-blue); color: #FFFFFF; border-color: var(--pp-blue);
  padding: 0.45rem 1.9rem;
}}
.stFormSubmitButton > button:hover {{ background: var(--pp-blue-dark); border-color: var(--pp-blue-dark); }}
.stDownloadButton > button:hover {{ border-color: var(--pp-blue); color: var(--pp-blue); }}
div[data-testid="stExpander"] details {{ border: 1px solid #E4E8ED; border-radius: 10px; background: #FFF; }}
code {{ background: #EFF2F6; color: var(--pp-blue-dark); border-radius: 4px; padding: 0.1rem 0.35rem; }}
.pp-cmd {{
  background: #FFFFFF; border: 1px solid #E4E8ED; border-left: 4px solid var(--pp-teal);
  border-radius: 0 8px 8px 0; padding: 0.75rem 1rem; font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 0.92rem; color: var(--pp-ink); overflow-x: auto; white-space: pre; margin-bottom: 0.9rem;
}}
.pp-legend {{ font-size: 0.95rem; color: var(--pp-muted); line-height: 2.1; }}
.pp-chip {{ padding: 0.16rem 0.6rem; border-radius: 5px; border: 1px solid #E0E4EA; margin-right: 0.45rem; }}
.pp-foot {{
  margin-top: 3rem; padding-top: 1.1rem; border-top: 1px solid #E4E8ED;
  font-size: 0.9rem; color: var(--pp-muted);
}}

/* ── static tables ───────────────────────────────────────────────────── */
/* Rendered as plain HTML rather than st.dataframe: every row and column is
   visible at once, headers wrap instead of being clipped, and there is no
   Arrow round-trip or data-grid to mount on each rerun. */
.pp-table {{ width: 100%; overflow-x: auto; margin-bottom: 0.6rem; }}
.pp-table table {{
  width: 100%; border-collapse: collapse; background: #FFFFFF;
  border: 1px solid #E4E8ED; border-radius: 10px; overflow: hidden;
  font-size: 0.95rem; font-variant-numeric: tabular-nums;
  table-layout: auto;
}}
.pp-table th {{
  background: #F2F5F9; color: var(--pp-ink); font-weight: 650; font-size: 0.9rem;
  text-align: right; padding: 0.55rem 0.6rem; border-bottom: 1px solid #DDE3EA;
  line-height: 1.25; vertical-align: bottom; hyphens: auto;
}}
.pp-table td {{
  padding: 0.5rem 0.6rem; text-align: right; color: var(--pp-ink);
  border-bottom: 1px solid #EEF1F5; white-space: nowrap;
}}
/* Dense: many narrow numeric columns (the prevalence grids, design previews)
   have to fit half the page width without a scrollbar. */
.pp-table.dense table {{ font-size: 0.85rem; }}
.pp-table.dense th {{ font-size: 0.82rem; padding: 0.4rem 0.4rem; }}
.pp-table.dense td {{ padding: 0.32rem 0.4rem; }}
/* Long text is allowed to wrap rather than push the table sideways. */
.pp-table.wrap-label td.pp-left {{ white-space: normal; }}
.pp-table.wrap-cells td {{ white-space: normal; word-break: break-word; }}
.pp-table th.pp-left, .pp-table td.pp-left {{ text-align: left; }}
.pp-table tbody tr:last-child td {{ border-bottom: none; }}
/* Footer band: the reference row a table is measured against, set apart from
   the ranking above it by a gap, a rule and a tint of its own. */
.pp-table tfoot tr.pp-gap td {{
  height: 10px; padding: 0; border: none; background: var(--pp-surface);
}}
.pp-table tfoot tr.pp-foot-row td {{
  border-top: 2px solid #C7D2DE; border-bottom: none; background: #F2F5F9;
  color: var(--pp-muted);
}}
.pp-table tfoot tr.pp-foot-row td.pp-left {{ font-style: italic; color: var(--pp-ink); }}
.pp-table caption {{
  caption-side: bottom; text-align: left; padding-top: 0.5rem;
  font-size: 0.9rem; color: var(--pp-muted);
}}
</style>
"""


def page_header(eyebrow: str, title: str, lead: str):
    st.markdown(
        f'<div class="pp-eyebrow">{eyebrow}</div>'
        f'<div class="pp-title">{title}</div>'
        f'<div class="pp-rule"></div>'
        f'<div class="pp-lead">{lead}</div>',
        unsafe_allow_html=True,
    )


def section(label: str):
    st.markdown(f'<div class="pp-section">{label}</div>', unsafe_allow_html=True)


def note(text: str, kind: str = ""):
    st.markdown(f'<div class="pp-note {kind}">{text}</div>', unsafe_allow_html=True)


def stat_strip(items: list):
    """items: list of (label, value, sub)."""
    html = '<div class="pp-stats">' + "".join(
        f'<div class="pp-stat"><div class="k">{k}</div><div class="v">{v}</div>'
        f'<div class="s">{s}</div></div>' for k, v, s in items
    ) + "</div>"
    st.markdown(html, unsafe_allow_html=True)


def csv_bytes(df: pd.DataFrame, index: bool = True) -> bytes:
    return df.to_csv(index=index).encode("utf-8")


def plural(n: int, word: str) -> str:
    return word if n == 1 else word + "s"


def _fmt_cell(value) -> str:
    """Numbers to at most two decimals, integers without a decimal point."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        return f"{value:,.0f}" if float(value).is_integer() else f"{value:,.2f}"
    return str(value)


def html_table(df: pd.DataFrame, cell_style=None, index_label: str | None = None,
               formatter=None, caption: str | None = None, dense: bool = False,
               wrap: str = "", footer: pd.DataFrame | None = None):
    """Render a DataFrame as a complete static table — no scrolling, no clipping.

    ``cell_style(row_position, column, value) -> css`` colours individual cells.
    ``dense`` tightens type and padding for wide numeric grids; ``wrap`` is
    "label" (first column may wrap) or "cells" (every cell may wrap), which is
    what keeps a long text column from forcing a horizontal scrollbar.
    ``footer`` rows are set below a rule, in their own band: they share the
    columns but are not part of the ranking above them.
    """
    formatter = formatter or _fmt_cell
    classes = "pp-table" + (" dense" if dense else "") + (f" wrap-{wrap}" if wrap else "")
    left = ' class="pp-left"'
    head = "".join(f'<th{left if i == 0 and index_label is None else ""}>{c}</th>'
                   for i, c in enumerate(df.columns))
    if index_label is not None:
        head = f'<th{left}>{index_label}</th>' + head

    def rows(frame: pd.DataFrame, style=None, row_class: str = "") -> str:
        out = []
        cls_attr = f' class="{row_class}"' if row_class else ""
        for pos, (label, row) in enumerate(frame.iterrows()):
            cells = [f"<td{left}>{label}</td>"] if index_label is not None else []
            for i, col in enumerate(df.columns):
                css = style(pos, col, row[col]) if style else ""
                cls = left if i == 0 and index_label is None else ""
                attr = f' style="{css}"' if css else ""
                cells.append(f"<td{cls}{attr}>{formatter(row[col])}</td>")
            out.append(f"<tr{cls_attr}>" + "".join(cells) + "</tr>")
        return "".join(out)

    cap = f"<caption>{caption}</caption>" if caption else ""
    foot = ""
    if footer is not None and not footer.empty:
        # An empty row in the page colour opens a real gap; the rule and the
        # tint then read as a separate band rather than one more ranked row.
        span = len(df.columns) + (1 if index_label is not None else 0)
        foot = (f'<tfoot><tr class="pp-gap"><td colspan="{span}"></td></tr>'
                f'{rows(footer, row_class="pp-foot-row")}</tfoot>')
    st.markdown(
        f'<div class="{classes}"><table>{cap}<thead><tr>{head}</tr></thead>'
        f'<tbody>{rows(df, cell_style)}</tbody>{foot}</table></div>',
        unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════
#  CHARTS
# ═══════════════════════════════════════════════════════════════════════════

CHART_FONT = "Inter, Helvetica Neue, sans-serif"
STATUS_DOMAIN = ["Fewer tests than individual", "No better than individual"]
STATUS_RANGE = [BLUE, CRITICAL]


def _chart_base(df: alt.Chart):
    return df.configure_view(strokeWidth=0).configure_axis(
        labelFont=CHART_FONT, titleFont=CHART_FONT, labelColor=MUTED, titleColor=INK,
        labelFontSize=14, titleFontSize=15, titleFontWeight=600, titlePadding=10,
        grid=True, gridColor="#ECEFF3", domainColor="#DDE2E8", tickColor="#DDE2E8",
    ).configure_legend(
        labelFont=CHART_FONT, titleFont=CHART_FONT, labelColor=INK, titleColor=MUTED,
        labelFontSize=14, titleFontSize=14, orient="top", direction="horizontal",
        symbolType="square", symbolSize=140, offset=6, title=None,
        labelLimit=0,   # the default 160px clipped "Fewer tests than individual"
    ).configure_title(font=CHART_FONT, fontSize=18, color=INK, anchor="start",
                      fontWeight=600, offset=14)


def _chart_frame(summary: pd.DataFrame, n_samp: int) -> pd.DataFrame:
    df = summary.copy()
    df["Tests"] = pd.to_numeric(df["Mean experiments"], errors="coerce")
    df["Pool size"] = pd.to_numeric(df["Max samples per pool"], errors="coerce")
    df["Rounds"] = pd.to_numeric(df["Mean steps"], errors="coerce")
    df["Saving"] = (100 * (1 - df["Tests"] / n_samp)).round(0)
    df["Status"] = np.where(df["Tests"] < n_samp, STATUS_DOMAIN[0], STATUS_DOMAIN[1])
    return df.dropna(subset=["Tests"])


def _status_colour(df: pd.DataFrame):
    """Fixed status→hue mapping, restricted to the states actually present.

    A single state needs no legend box — the chart title and caption name it.
    """
    present = [s for s in STATUS_DOMAIN if (df["Status"] == s).any()]
    colours = [STATUS_RANGE[STATUS_DOMAIN.index(s)] for s in present]
    return alt.Color("Status:N",
                     scale=alt.Scale(domain=present, range=colours),
                     legend=alt.Legend() if len(present) > 1 else None)


def _tooltip():
    return [alt.Tooltip("Pooling strategy:N", title="Strategy"),
            alt.Tooltip("Tests:Q", title="Mean tests", format=".1f"),
            alt.Tooltip("Saving:Q", title="Saving vs individual (%)", format=".0f"),
            alt.Tooltip("Pool size:Q", title="Max samples per pool", format=".0f"),
            alt.Tooltip("Rounds:Q", title="Mean steps", format=".2f")]


def tests_chart(summary: pd.DataFrame, n_samp: int):
    df = _chart_frame(summary, n_samp)
    if df.empty:
        return None

    # Explicit category order — a layered chart cannot resolve sort="x" reliably.
    order = df.sort_values("Tests")["Pooling strategy"].tolist()
    y = alt.Y("Pooling strategy:N", sort=order, title=None,
              axis=alt.Axis(labelOverlap=False, labelLimit=200, labelPadding=8,
                            labelFontSize=14, labelColor=INK))

    bars = alt.Chart(df).mark_bar(cornerRadiusEnd=4, size=18).encode(
        x=alt.X("Tests:Q", title="Mean number of tests"),
        y=y, color=_status_colour(df), tooltip=_tooltip(),
    )
    values = alt.Chart(df).mark_text(
        align="left", dx=7, font=CHART_FONT, fontSize=13, color=MUTED).encode(
        x="Tests:Q", y=y, text=alt.Text("Tests:Q", format=".1f"))
    baseline = alt.Chart(pd.DataFrame({"x": [n_samp]})).mark_rule(
        color=INK, strokeDash=[5, 4], strokeWidth=1.5, opacity=0.6).encode(x="x:Q")
    label = alt.Chart(
        pd.DataFrame({"x": [n_samp], "t": [f"Individual testing — {n_samp} tests"]})).mark_text(
        align="right", dx=-7, dy=-4, baseline="top", color=MUTED, font=CHART_FONT,
        fontSize=13).encode(x="x:Q", text="t:N")

    chart = (bars + values + baseline + label).properties(
        height=alt.Step(34), title="Tests required per strategy")
    return _chart_base(chart)


def _pareto_front(df: pd.DataFrame) -> pd.DataFrame:
    """Designs not beaten on both axes at once — the usable trade-off frontier."""
    keep = []
    for i, row in df.iterrows():
        dominated = ((df["Tests"] <= row["Tests"]) & (df["Pool size"] <= row["Pool size"]) &
                     ((df["Tests"] < row["Tests"]) | (df["Pool size"] < row["Pool size"])))
        if not dominated.any():
            keep.append(i)
    return df.loc[keep].sort_values("Tests")


TRADEOFF_HEIGHT = 430
LABEL_GAP = 24          # px of vertical room one label needs to stay readable
LABEL_WIDTH_HINT = 760  # px assumed when placing labels; the chart itself is responsive
LABEL_CHAR_PX = 7.6     # mean advance of the 13px semibold label face


def _dodge_labels(values: np.ndarray, lo: float, hi: float, height: int,
                  gap: int = LABEL_GAP) -> np.ndarray:
    """Push labels apart vertically until none of them overlaps.

    Works in pixels against a known y domain, then hands the positions back in
    data units. Alternating fixed offsets — what this chart used before — only
    helps when neighbours alternate; a cluster of three or more frontier points
    at a similar pool size still collided.
    """
    span = (hi - lo) or 1.0
    px = height * (1 - (values - lo) / span)          # data → pixels, y grows down
    order = np.argsort(px, kind="stable")

    placed = px.copy()
    last = -np.inf
    for i in order:                                    # top to bottom, keep the gap
        last = placed[i] = max(placed[i], last + gap)

    overflow = placed[order[-1]] - height              # ran off the bottom: shift up
    if overflow > 0:
        placed -= overflow
        last = np.inf
        for i in order[::-1]:                          # bottom to top, keep the gap
            last = placed[i] = min(placed[i], last - gap)
        placed = np.maximum(placed, 0)

    return lo + (1 - placed / height) * span


def _label_alignments(front: pd.DataFrame, df: pd.DataFrame, x_lo, x_hi, y_lo, y_hi) -> list:
    """Put each label on whichever side of its point is clear of other markers.

    Dodging keeps labels off each other; this keeps them off the *points*, which
    is the other half of the clutter on a scatter this dense.
    """
    w, h = LABEL_WIDTH_HINT, TRADEOFF_HEIGHT
    to_px_x = lambda v: (v - x_lo) / ((x_hi - x_lo) or 1) * w
    to_px_y = lambda v: (1 - (v - y_lo) / ((y_hi - y_lo) or 1)) * h

    px = to_px_x(df["Tests"].to_numpy(dtype=float))
    py = to_px_y(df["Pool size"].to_numpy(dtype=float))

    alignments = []
    for _, row in front.iterrows():
        lx, ly = to_px_x(row["Tests"]), to_px_y(row["LabelY"])
        reach = len(str(row["Pooling strategy"])) * LABEL_CHAR_PX + 16
        same_band = np.abs(py - ly) < 11              # close enough vertically to clash
        after = np.any(same_band & (px > lx + 6) & (px < lx + reach)) or lx + reach > w
        before = np.any(same_band & (px < lx - 6) & (px > lx - reach)) or lx - reach < 0
        alignments.append("right" if after and not before else "left")
    return alignments


def tradeoff_chart(summary: pd.DataFrame, n_samp: int):
    df = _chart_frame(summary, n_samp).dropna(subset=["Pool size"])
    if df.empty:
        return None
    front = _pareto_front(df).copy()

    # Explicit padded domains rather than nice/padding: the label layout below
    # needs to know exactly which pixel a value lands on.
    def domain(series, pad=0.08):
        lo, hi = float(series.min()), float(series.max())
        margin = (hi - lo) * pad or max(abs(hi) * pad, 1.0)
        return lo - margin, hi + margin

    x_lo, x_hi = domain(df["Tests"])
    y_lo, y_hi = domain(df["Pool size"])
    y_title = "Max samples per pool"
    x = alt.X("Tests:Q", title="Mean number of tests",
              scale=alt.Scale(domain=[x_lo, x_hi], nice=False, clamp=True))
    y = alt.Y("Pool size:Q", title=y_title,
              scale=alt.Scale(domain=[y_lo, y_hi], nice=False, clamp=True))

    front["LabelY"] = _dodge_labels(front["Pool size"].to_numpy(dtype=float),
                                    y_lo, y_hi, TRADEOFF_HEIGHT)
    front["Side"] = _label_alignments(front, df, x_lo, x_hi, y_lo, y_hi)
    # A connector is only drawn where dodging actually moved the label away.
    front["Moved"] = (np.abs(front["LabelY"] - front["Pool size"])
                      > 0.02 * (y_hi - y_lo))

    frontier = alt.Chart(front).mark_line(
        color=MUTED, strokeWidth=1.5, strokeDash=[4, 3], opacity=0.55,
        interpolate="step-after").encode(x=x, y=y)
    leaders = alt.Chart(front[front["Moved"]]).mark_rule(
        color=MUTED, strokeWidth=0.9, opacity=0.4).encode(
        x=x, y=y, y2=alt.Y2("LabelY:Q"))
    # 2px surface ring so overlapping marks stay separable
    halo = alt.Chart(df).mark_point(filled=True, size=230, opacity=1, color=SURFACE).encode(x=x, y=y)
    dots = alt.Chart(df).mark_point(filled=True, size=150, opacity=0.95).encode(
        x=x, y=y, color=_status_colour(df), tooltip=_tooltip())
    # Direct-label the frontier only; everything else is one hover away. Each
    # label is drawn twice: a surface-coloured outline first, so a name that has
    # to pass close to a marker still reads cleanly, then the text itself.
    # The label layer plots a second field on the shared y scale, so it has to
    # repeat the axis title verbatim: layered axes are merged, and either None or
    # a different string here loses the axis for the whole chart.
    label_y = alt.Y("LabelY:Q", title=y_title,
                    scale=alt.Scale(domain=[y_lo, y_hi], nice=False, clamp=True))
    labels = []
    for side, part in front.groupby("Side"):
        base = alt.Chart(part)
        common = dict(align=side, dx=(13 if side == "left" else -13), baseline="middle",
                      font=CHART_FONT, fontSize=13, fontWeight=600)
        encode = dict(x=x, y=label_y, text="Pooling strategy:N")
        labels.append(base.mark_text(**common, color=SURFACE, stroke=SURFACE,
                                     strokeWidth=4, opacity=0.95).encode(**encode))
        labels.append(base.mark_text(**common, color=INK).encode(**encode))

    chart = alt.layer(frontier, leaders, halo, dots, *labels).properties(
        height=TRADEOFF_HEIGHT,
        title="Efficiency against pool size — fewer tests means larger pools")
    return _chart_base(chart)


def summary_cell_style(df: pd.DataFrame, n_samp: int):
    """Green = among the three most efficient; red row = no better than individual testing."""
    numeric = pd.to_numeric(df["Mean experiments"], errors="coerce").to_numpy(dtype=float)
    ranked = np.argsort(numeric, kind="stable")[:3]
    best = {int(i) for i in ranked if numeric[i] < n_samp}
    over = {int(i) for i in np.flatnonzero(numeric > n_samp)}

    def _style(row_position, column, _value):
        if column == "Mean experiments" and row_position in best:
            return f"background-color: {GOOD_FILL}; font-weight: 600"
        if row_position in over:
            return f"background-color: {BAD_FILL}"
        return ""

    return _style


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — DESIGN
# ═══════════════════════════════════════════════════════════════════════════

def render_design():
    page_header(
        "Method comparison",
        "Find the right pooling design",
        "Enter how many samples you need to screen and how many positives you are willing to "
        "resolve in one batch. PoolPy benchmarks ten conceptually different pooling algorithms "
        "against one-sample-one-test and returns the designs ready to run.",
    )

    with st.form("design_form"):
        c1, c2 = st.columns(2)
        n_samp = c1.number_input("Number of samples (S)", min_value=2, max_value=100_000,
                                 value=100, step=1)
        diff = c2.number_input("Max. number of positives (D)", min_value=1, max_value=1_000,
                               value=1, step=1)

        # Keep the panel open once imperfect values are in play, so the settings
        # behind the current table stay visible.
        previous = st.session_state.get("design", {})
        imperfect = previous.get("sens", 1.0) < 1.0 or previous.get("spec", 1.0) < 1.0
        with st.expander("Advanced options", expanded=imperfect):
            st.caption("A perfect test is assumed by default. Lower either value to see how "
                       "test errors propagate through each design; two extra columns then "
                       "appear in the benchmark table.")
            c3, c4 = st.columns(2)
            sens = c3.number_input("Test sensitivity", min_value=0.0, max_value=1.0,
                                   value=1.0, step=0.01,
                                   help="Per-test probability of detecting a true positive.")
            spec = c4.number_input("Test specificity", min_value=0.0, max_value=1.0,
                                   value=1.0, step=0.01,
                                   help="Per-test probability of correctly calling a negative.")
        submitted = st.form_submit_button("Benchmark designs")

    if submitted:
        st.session_state["design"] = dict(n=int(n_samp), d=int(diff),
                                          sens=float(sens), spec=float(spec))

    params = st.session_state.get("design")
    if not params:
        note("Set your parameters above and select <b>Benchmark designs</b>. "
             f"Results come from a database of benchmarked designs for S ≤ {MAX_N} and "
             f"D ≤ {MAX_DIFFERENTIATE}; anything outside that range is evaluated analytically "
             "on the fly.")
        return

    n_samp, diff = params["n"], params["d"]
    sens, spec = params["sens"], params["spec"]

    if diff >= n_samp:
        note(f"The maximum number of positives (<b>{diff}</b>) must be smaller than the "
             f"number of samples (<b>{n_samp}</b>).", "bad")
        return

    summary, source, precomp_cell = None, None, None
    prevalence = diff / n_samp

    if diff <= MAX_DIFFERENTIATE and n_samp <= MAX_N and prevalence <= MAX_PREVALENCE:
        grid = precomputed_grid()
        if grid is None:
            note("The benchmark database was not found. Expected "
                 f"<code>{PRECOMPUTED_PARQUET.name}</code> under <code>data/</code> or the "
                 "original <code>precomputed/</code> folder tree. Falling back to on-the-fly "
                 "evaluation.", "warn")
        else:
            cell = closest_precomputed(grid, n_samp, diff)
            raw = precomputed_cell(*cell) if cell else None
            if raw is not None:
                precomp_cell = cell
                summary = process_metrics_data(raw, sens, spec)
                source = "database"

    if summary is None:
        if n_samp * diff * diff > MAX_PROD:
            reason = (f"D = {diff} is above the benchmarked maximum of {MAX_DIFFERENTIATE}"
                      if diff > MAX_DIFFERENTIATE else
                      f"S = {n_samp} is above the benchmarked maximum of {MAX_N}"
                      if n_samp > MAX_N else
                      f"the prevalence D/S = {prevalence:.1%} is above the benchmarked "
                      f"maximum of {MAX_PREVALENCE:.0%}")
            note(f"No benchmark table for these parameters, because {reason}. "
                 f"Evaluating them on the fly would also be too expensive "
                 f"(S·D² = {n_samp * diff * diff:,} exceeds the {MAX_PROD:,} limit). "
                 "You can still download the designs below, or generate everything locally "
                 "with the command at the bottom of this page.", "warn")
        else:
            with st.spinner("Evaluating every strategy analytically…"):
                summary = process_metrics_data(cached_fly_summary(n_samp, diff), sens, spec)
            source = "onthefly"

    # ── result banner ────────────────────────────────────────────────────
    if source == "database":
        new_n, new_d = precomp_cell
        if (new_n, new_d) == (n_samp, diff):
            note(f"Benchmarked results for <b>{n_samp} samples</b> with up to "
                 f"<b>{diff} {plural(diff, 'positive')}</b>.", "good")
        else:
            note(f"No benchmark table exists for exactly {n_samp} samples with up to "
                 f"{diff} {plural(diff, 'positive')}. Showing the closest benchmarked cell: "
                 f"<b>S = {new_n}, D = {new_d}</b>. The designs offered below are still "
                 "generated for your exact parameters.", "warn")
    elif source == "onthefly":
        note(f"These parameters lie outside the benchmark database, so the metrics below were "
             f"derived analytically for <b>S = {n_samp}</b>, <b>D = {diff}</b>. They are "
             "approximations — the hierarchical row in particular comes from a short "
             "Monte-Carlo search.", "warn")

    if summary is not None and not summary.empty:
        _render_summary(summary, n_samp, diff, sens, spec)

    _render_design_downloads(n_samp, diff)
    _render_local_commands(n_samp, diff)


@st.fragment
def _render_summary(summary: pd.DataFrame, n_samp: int, diff: int, sens: float, spec: float):
    numeric = pd.to_numeric(summary["Mean experiments"], errors="coerce")
    viable = summary[numeric < n_samp]

    if not viable.empty:
        fewest = viable.sort_values(["Mean experiments", "Mean steps"]).iloc[0]
        quickest = viable.sort_values(["Mean steps", "Mean experiments"]).iloc[0]
        smallest = viable.sort_values(["Max samples per pool", "Mean experiments"]).iloc[0]
        best_tests = float(fewest["Mean experiments"])
        stat_strip([
            ("Fewest tests", f"{best_tests:,.0f}", f"{fewest['Pooling strategy']}"),
            ("Saving vs individual", f"{100 * (1 - best_tests / n_samp):.0f}%",
             f"{n_samp - best_tests:,.0f} tests avoided"),
            ("Fastest to answer", f"{float(quickest['Mean steps']):.2f}",
             f"mean steps — {quickest['Pooling strategy']}"),
            ("Smallest pools", f"{float(smallest['Max samples per pool']):,.0f}",
             f"samples/pool — {smallest['Pooling strategy']}"),
        ])
        note("".join(f'<div class="pp-note-line">{line}</div>' for line in (
            f"If you mainly need to <b>reduce the number of tests</b>, consider the "
            f"<b>{fewest['Pooling strategy']}</b> design.",
            f"If you need results <b>quickly</b>, consider the "
            f"<b>{quickest['Pooling strategy']}</b> design.",
            f"If you <b>cannot pool many samples</b> together, consider the "
            f"<b>{smallest['Pooling strategy']}</b> design.",
        )))
    else:
        note("No strategy beats individual testing at these parameters. That usually means the "
             "prevalence is too high for pooling to pay off — try a smaller D, or split the "
             "samples into smaller batches.", "warn")

    section("Benchmark table")
    # Individual testing is the reference the pooling designs are judged against,
    # not one of them, so it sits below the rule in its own band.
    baseline = single_sample_row(summary, n_samp, diff, sens, spec)
    table = pd.concat([summary, baseline], ignore_index=True)      # the CSV keeps both
    # Strategy names stay on one line: measured in the browser, letting them wrap
    # only made rows ragged without buying width the table needed.
    html_table(summary, cell_style=summary_cell_style(summary, n_samp),
               dense=len(summary.columns) > 8, footer=baseline)
    st.markdown(
        f'<div class="pp-legend">'
        f'<span class="pp-chip" style="background:{GOOD_FILL}">Among the three most efficient</span>'
        f'<span class="pp-chip" style="background:{BAD_FILL}">No better than individual testing</span>'
        f'Below the rule is the one-sample-one-test baseline every design is compared '
        f'against. Column definitions are in the Guide.'
        f'</div>', unsafe_allow_html=True)

    st.download_button("Download benchmark table (CSV)", csv_bytes(table, index=False),
                       file_name=f"poolpy_summary_S{n_samp}_D{diff}.csv", mime="text/csv",
                       key="dl_summary")

    section("Trade-offs")
    c1 = tests_chart(summary, n_samp)
    if c1 is not None:
        st.altair_chart(c1, width="stretch")
    c2 = tradeoff_chart(summary, n_samp)
    if c2 is not None:
        st.altair_chart(c2, width="stretch")
        st.caption("There is no universally best pooling strategy — designs that need the "
                   "fewest tests generally require the largest pools, which leads to "
                   "higher signal dilution. Named points sit on the trade-off frontier (dashed "
                   "line): nothing else does better on both axes at once. Hover any point for "
                   "its full metrics.")


@st.fragment
def _render_design_downloads(n_samp: int, diff: int):
    section("Download designs for your exact parameters")
    st.caption("Designs are generated on demand for S = %d and D = %d, as a samples × pools "
               "binary matrix. These same CSVs feed the Decoder and Automation pages."
               % (n_samp, diff))

    families = [
        ("Positive-count dependent", "dependent",
         "These designs adapt to the number of expected positives."),
        ("Positive-count independent", "independent",
         "These designs do not adapt to D — use them with care above one positive."),
    ]
    for title, family, blurb in families:
        st.markdown(f"**{title}**  \n<span style='color:{MUTED};font-size:0.95rem'>{blurb}</span>",
                    unsafe_allow_html=True)
        labels = [k for k, v in DESIGN_BUILDERS.items() if v[2] == family]
        cols = st.columns(len(labels))
        for col, label in zip(cols, labels):
            with col:
                if label == "Ch. Rm. special" and diff not in (2, 3):
                    st.button(label, disabled=True, key=f"dis_{label}",
                              help="Available only for D = 2 or D = 3.", width="stretch")
                    continue
                try:
                    design = cached_design(label, n_samp, diff)
                    st.download_button(
                        label, deferred_csv(design),
                        file_name=f"{label.replace('. ', '_').replace(' ', '_')}"
                                  f"_S{n_samp}_D{diff}.csv",
                        mime="text/csv", key=f"dl_{label}", width="stretch",
                        help=f"{design.shape[0]} samples × {design.shape[1]} pools")
                except Exception as exc:
                    st.button(label, disabled=True, key=f"err_{label}",
                              help=f"Not available: {exc}", width="stretch")
        st.write("")

    note("Every benchmarked design is archived on Zenodo — "
         f'<a href="{ZENODO_URL}" target="_blank">DOI: 10.5281/zenodo.18660061</a>.')


def _render_local_commands(n_samp: int, diff: int):
    section("Run locally")
    st.caption("For large screens, or to reproduce the benchmark yourself, run the same code "
               "from the repository. Use the notebooks `pool_interface.ipynb` and "
               "`decode_interface.ipynb`, or the command line:")
    st.markdown(
        f'<div class="pp-cmd">python pool_N.py --n_samp {n_samp} --differentiate {diff} '
        f'--directory your/path</div>'
        f'<div class="pp-cmd">python decode_N.py --differentiate {diff} '
        f'--path_to_WA your/design.csv --readout your/readout.csv</div>',
        unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — PREVALENCE
# ═══════════════════════════════════════════════════════════════════════════

def render_prevalence():
    page_header(
        "Parameter choice",
        "Turn a prevalence estimate into pooling parameters",
        "Two complementary tools. The first sizes pools for sampling from a large population of "
        "known prevalence. The second reports how often a chosen maximum number of positives (D) "
        "would be exceeded, so you can pick D and the number of batches with a stated error rate.",
    )

    _prevalence_optimiser()
    _prevalence_risk()


@st.fragment
def _prevalence_optimiser():
    section("Optimal pool size at a known prevalence")
    with st.form("bg2_form"):
        c1, c2, c3 = st.columns([2, 2, 1])
        p = c1.number_input("Prevalence (as a fraction)", min_value=1e-6, max_value=0.5,
                            value=0.01, step=0.001, format="%.5f")
        steps = c2.number_input("Max. hierarchical steps", min_value=1, max_value=100,
                                value=10, step=1,
                                help="Upper bound on the number of successive splitting rounds.")
        c3.markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
        go = c3.form_submit_button("Optimise")

    if go:
        st.session_state["prev_opt"] = (float(p), int(steps))

    if "prev_opt" in st.session_state:
        pp, ss = st.session_state["prev_opt"]
        with st.spinner("Searching the strategy space…"):
            hier, md = cached_optimal_strategies(pp, ss)

        c1, c2 = st.columns(2)
        with c1:
            if hier is None:
                st.markdown('<div class="pp-card"><h4>Optimal hierarchical strategy</h4>'
                            '<p>No strategy beats individual testing at this prevalence.</p>'
                            '</div>', unsafe_allow_html=True)
            else:
                fv0, e0 = hier
                splits = [int(x) for x in fv0]
                st.markdown(
                    f'<div class="pp-card"><h4>Optimal hierarchical strategy</h4>'
                    f'<p>Split factors per level: <b>{splits}</b></p>'
                    f'<p>Initial pool size: <b>{int(np.prod(fv0))}</b> samples</p>'
                    f'<p>Expected tests per sample: <b>{e0:.4f}</b> '
                    f'(individual testing = 1.0000)</p>'
                    f'<p style="margin-top:0.6rem">At each step, positive pools are subdivided '
                    f'by the next factor until samples are tested individually.</p></div>',
                    unsafe_allow_html=True)
        with c2:
            if md is None:
                st.markdown('<div class="pp-card"><h4>Optimal multidimensional strategy</h4>'
                            '<p>No strategy beats individual testing at this prevalence.</p>'
                            '</div>', unsafe_allow_html=True)
            else:
                D, N, ET = int(md[0]), int(md[1]), md[2]
                st.markdown(
                    f'<div class="pp-card"><h4>Optimal multidimensional strategy</h4>'
                    f'<p>Dimensions: <b>{D}</b> &nbsp;·&nbsp; side length: <b>{N}</b></p>'
                    f'<p>Total pool size: <b>{N ** D}</b> samples</p>'
                    f'<p>Expected tests per sample: <b>{ET:.4f}</b> '
                    f'(individual testing = 1.0000)</p>'
                    f'<p style="margin-top:0.6rem">Arrange {N ** D} samples in a {D}-dimensional '
                    f'{N}×… grid; every coordinate on every axis defines one pool.</p></div>',
                    unsafe_allow_html=True)


@st.fragment
def _prevalence_risk():
    section("Misparametrisation risk for a chosen D")
    st.caption("The tables give the probability that at least one pool contains more than D "
               "positives — per pool on the left, and across all pools of a combinatorial batch "
               "on the right (family-wise error rate). Rows are sample counts, columns are D.")

    with st.form("prev_form"):
        c1, c2, c3, c4 = st.columns([2, 2, 2, 1])
        n = c1.number_input("Number of samples", min_value=5, max_value=100_000,
                            value=100, step=1)
        prev = c2.number_input("Prevalence (as a fraction)", min_value=1e-6, max_value=0.5,
                               value=0.01, step=0.001, format="%.5f")
        max_err = c3.number_input("Max. misparametrisation probability", min_value=1e-6,
                                  max_value=1.0, value=0.05, step=0.01, format="%.4f")
        c4.markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
        go2 = c4.form_submit_button("Evaluate")

    if go2:
        st.session_state["prev_tbl"] = (int(n), float(prev), float(max_err))

    if "prev_tbl" in st.session_state:
        nn, pv, me = st.session_state["prev_tbl"]
        df_orig, df_corr = cached_error_tables(nn, pv, me)

        def colour(v_orig, v_corr):
            if v_orig > me and v_corr > me:
                return f"background-color: {BAD_FILL}"
            if v_corr > me:
                return f"background-color: {WARN_FILL}"
            return f"background-color: {GOOD_FILL}"

        # Both tables share one colour scale: every cell is judged on its own value
        # and on its counterpart in the other table, which sit at the same position.
        def cell_colour(row_position, column, _value):
            j = df_orig.columns.get_loc(column)
            return colour(df_orig.iloc[row_position, j], df_corr.iloc[row_position, j])

        def probability(value) -> str:
            # 1.23e-05 → 1.23e-5: three significant digits in as few characters
            # as possible, so six columns fit half the page without scrolling.
            return "0" if value == 0 else f"{value:.3g}".replace("e-0", "e-")

        # The batch multiplier reads "16 X 25" in the source table; a thin space
        # around × keeps it one unbreakable token.
        corrected = df_corr.set_axis([str(i).replace(" X ", "×") for i in df_corr.index])

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Single batch**")
            html_table(df_orig, cell_style=cell_colour, index_label="N", formatter=probability,
                       dense=True, caption="Columns are the maximum number of positives D.")
        with c2:
            st.markdown("**Corrected across batches (FWER)**")
            html_table(corrected, cell_style=cell_colour, index_label="Batches×N",
                       formatter=probability, dense=True,
                       caption="Columns are the maximum number of positives D.")

        st.markdown(
            f'<div class="pp-legend">'
            f'<span class="pp-chip" style="background:{GOOD_FILL}">Within the target rate '
            f'both per pool and across batches</span><br>'
            f'<span class="pp-chip" style="background:{WARN_FILL}">Within the target rate per '
            f'pool, but not across batches</span><br>'
            f'<span class="pp-chip" style="background:{BAD_FILL}">Above the target rate in both '
            f'cases</span></div>', unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — DECODER
# ═══════════════════════════════════════════════════════════════════════════

TYPE_COPY = {
    "unique": ("good", "Unique solution"),
    "multiple": ("warn", "Multiple solutions"),
    "putative": ("warn", "Putative positives"),
    "error": ("bad", "No solution"),
}


def render_decoder():
    page_header(
        "Result interpretation",
        "Decode a pooled experiment",
        "Upload the design you ran, then tell PoolPy which pools came back positive. "
        "The decoder returns the individual samples consistent with that readout, either as a "
        "single solution, a short list of candidate combinations, or a putative set to confirm.",
    )

    design_file = st.file_uploader("Pooling design (CSV generated by PoolPy)", type=["csv"],
                                   key="dec_design")
    wa_df = None
    if design_file is not None:
        try:
            wa_df = read_design_csv(design_file.getvalue())
        except Exception as exc:
            note(f"Could not read the file: {exc}", "bad")
            wa_df = None
        if wa_df is not None:
            ok, msg = validate_wa_df(wa_df)
            note(msg, "good" if ok else "bad")
            if not ok:
                wa_df = None

    if wa_df is None:
        note("Only design CSVs produced by PoolPy are accepted: columns named "
             "<code>Pool 0…</code>, rows named <code>Sample 0…</code>, and 0/1 entries. "
             "Download one from the Design page to try this out.")
        return

    wa = wa_df.values
    n_compounds, n_pools = wa.shape

    with st.expander(f"Preview the design — {n_compounds} samples × {n_pools} pools"):
        # A corner of the matrix is enough; sending a 1000 × 800 grid to the
        # browser just to preview it is what made this expander sluggish.
        preview = wa_df.iloc[:20, :20]
        html_table(preview, index_label="", dense=True, caption=(
            f"First {preview.shape[0]} of {n_compounds} samples × "
            f"first {preview.shape[1]} of {n_pools} pools."))

    section("Readout")
    c1, c2 = st.columns([3, 2])
    with c1:
        mode = st.radio("Readout mode", ["Single readout", "Multiple readouts (CSV)"],
                        horizontal=True, label_visibility="collapsed")
    with c2:
        diff_txt = st.text_input(
            "Max. number of positives", value="", placeholder="leave empty for no limit",
            help="Bounds the combination search. Left empty the decoder considers any number "
                 "of positives, which usually exhausts the search budget and falls back to a "
                 "putative set. Setting it to the D you designed for gives sharper calls.")
    try:
        diff = int(diff_txt) if diff_txt.strip() else n_compounds
    except ValueError:
        note("The maximum number of positives must be an integer. Ignoring it.", "warn")
        diff = n_compounds
    if diff < 1:
        note("The maximum number of positives must be at least 1. Ignoring it.", "warn")
        diff = n_compounds

    if mode == "Single readout":
        _decode_single(wa, n_pools, n_compounds, diff)
    else:
        _decode_multi(wa, n_pools, n_compounds, diff)


@st.fragment
def _decode_single(wa, n_pools, n_compounds, diff):
    readout_str = st.text_input(
        "Positive pools", value="", placeholder="e.g. 1,3,5",
        help="Comma-separated indices of the pools that tested positive.")
    signature = (wa.shape, readout_str, diff)

    if st.button("Decode", type="primary"):
        if not readout_str.strip():
            note("Enter at least one positive pool.", "bad")
            return
        tokens = [t.strip() for t in readout_str.split(",") if t.strip()]
        if any(not t.isdigit() for t in tokens):
            note("Only non-negative integers separated by commas are accepted.", "bad")
            return
        with st.spinner("Searching sample combinations…"):
            st.session_state["dec_single"] = (signature,
                                              decode_readout([int(t) for t in tokens], wa, diff))

    stored = st.session_state.get("dec_single")
    # Results are kept across the reruns that download buttons trigger, but only
    # while the design, readout and D are unchanged.
    if not stored or stored[0] != signature:
        return
    result = stored[1]
    kind, title = TYPE_COPY[result["type"]]

    if result["type"] == "error":
        note(f"<b>{title}.</b> {result['note']}", "bad")
        return

    if result["type"] == "unique":
        samples = ", ".join(str(s) for s in result["samples"])
        body = (f"<b>{title}.</b> The readout is explained by exactly one set of positive "
                f"samples: <b>Sample {samples}</b>.")
    elif result["type"] == "multiple":
        combos = "<br>".join("&nbsp;&nbsp;· Samples " + ", ".join(map(str, c))
                             for c in result["combinations"])
        body = (f"<b>{title}.</b> {len(result['combinations'])} sets of positive samples are "
                f"consistent with this readout:<br>{combos}")
    else:
        body = (f"<b>{title}.</b> {result['note']} There are up to "
                f"{min(diff, len(result['samples']))} positives among samples "
                f"<b>{', '.join(map(str, result['samples']))}</b>. Test them individually, or "
                f"use a design that resolves more positives.")
    note(body, kind)

    text = (f"PoolPy decoder\nReadout (positive pools): {readout_str}\n"
            f"Design: {n_compounds} samples x {n_pools} pools\n"
            f"Max. positives: {diff}\nResult type: {result['type']}\n"
            f"Positive samples: {result['samples']}\n")
    if result["type"] == "multiple":
        text += "Candidate combinations:\n" + "\n".join(
            "  " + ", ".join(map(str, c)) for c in result["combinations"]) + "\n"
    st.download_button("Download decoder output (TXT)", text.encode(),
                       file_name="decoder_output.txt", mime="text/plain")


@st.fragment
def _decode_multi(wa, n_pools, n_compounds, diff):
    readout_file = st.file_uploader(
        "Readout CSV — one row per readout, each listing the positive pools",
        type=["csv"], key="dec_readout")
    if readout_file is None:
        note("Each row may be a comma-separated list of positive pools (<code>1,3,5</code>) or a "
             f"binary vector of length {n_pools}. Header rows and ID columns are detected and "
             "skipped automatically.")
        return

    raw = readout_file.getvalue().decode("utf-8", errors="replace")
    signature = (wa.shape, hash(raw), diff)

    if st.button("Decode all readouts", type="primary"):
        # csv.reader, not str.split — a cell like "1,3,5" is one quoted field.
        rows = [r for r in csv.reader(io.StringIO(raw)) if any(c.strip() for c in r)]
        if not rows:
            note("The readout file appears to be empty.", "bad")
            return
        width = max(len(r) for r in rows)
        df = pd.DataFrame([r + [""] * (width - len(r)) for r in rows])

        dropped_row = dropped_col = False
        if df.shape[0] > 0:
            checks = [v for v in (cell_is_valid_token(x, n_pools) for x in df.iloc[0])
                      if v is not None]
            if checks and not all(checks):
                df, dropped_row = df.iloc[1:].reset_index(drop=True), True
        if df.shape[1] > 0:
            checks = [v for v in (cell_is_valid_token(x, n_pools) for x in df.iloc[:, 0])
                      if v is not None]
            if checks and not all(checks):
                df, dropped_col = df.iloc[:, 1:].reset_index(drop=True), True

        results, table = [], []
        total = max(len(df), 1)
        # One socket message per readout would cost more than the decoding for
        # large files; 40 updates are enough to look continuous.
        step = max(1, total // 40)
        progress = st.progress(0.0, text="Decoding readouts…")
        for i, (_, row) in enumerate(df.iterrows()):
            parsed, err = series_to_readout_list(row, n_pools)
            res = ({"type": "error", "note": err} if parsed is None
                   else decode_readout(parsed, wa, diff))
            results.append(res)
            if res["type"] == "error":
                output = res.get("note", "error")
            elif res["type"] == "multiple":
                output = " | ".join(", ".join(map(str, c)) for c in res["combinations"])
            else:
                output = ", ".join(map(str, res["samples"]))
            table.append({"decoded_type": res["type"], "decoder_output": output})
            if (i + 1) % step == 0 or i + 1 == total:
                progress.progress((i + 1) / total)
        progress.empty()

        out_df = pd.DataFrame(table)
        out_df.index = [f"Readout {i}" for i in range(len(out_df))]
        st.session_state["dec_multi"] = (signature, out_df, results,
                                         dropped_row, dropped_col)

    stored = st.session_state.get("dec_multi")
    if not stored or stored[0] != signature:
        return
    _, out_df, results, dropped_row, dropped_col = stored

    skipped = []
    if dropped_row:
        skipped.append("the first row (header)")
    if dropped_col:
        skipped.append("the first column (identifiers)")
    note(f"Decoded <b>{len(out_df)}</b> readouts"
         + (f", skipping {' and '.join(skipped)}." if skipped else "."), "good")

    counts = out_df["decoded_type"].value_counts()
    stat_strip([(TYPE_COPY[t][1], f"{counts.get(t, 0)}", f"of {len(out_df)} readouts")
                for t in ("unique", "multiple", "putative", "error") if t in counts.index])

    if counts.get("putative", 0) > len(out_df) / 2 and diff >= n_compounds:
        note("Most readouts came back putative because the maximum number of positives is "
             "unbounded, so the combination search runs out of budget. Set it to the D you "
             "designed for and decode again.", "warn")

    # A decoded row can name dozens of samples, so that column wraps instead of
    # widening the table. Long batches are capped rather than dumped in full —
    # the CSV below always holds every row.
    limit = 200
    display = out_df.rename(columns={"decoded_type": "Type",
                                     "decoder_output": "Positive samples"}).head(limit)
    display.index = [i.replace("Readout ", "") for i in display.index]
    html_table(display, index_label="Readout", wrap="cells",
               caption=(f"Showing the first {limit} of {len(out_df):,} readouts — "
                        "download the CSV for all of them."
                        if len(out_df) > limit else None))

    probs = sample_probability_matrix(results, n_compounds)
    c1, c2 = st.columns(2)
    c1.download_button("Download decoded readouts (CSV)", csv_bytes(out_df),
                       file_name="decoded_readouts.csv", mime="text/csv", width="stretch")
    c2.download_button("Download sample probabilities (CSV)", csv_bytes(probs),
                       file_name="sample_probabilities.csv", mime="text/csv", width="stretch")
    st.caption("Sample probabilities give, per readout, the fraction of consistent combinations "
               "containing each sample: 1 for a unique call, 0.5 for putative sets.")


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — AUTOMATION
# ═══════════════════════════════════════════════════════════════════════════

def render_automation():
    page_header(
        "Liquid handling",
        "Turn a design into a robot protocol",
        "Upload a PoolPy design and get the full list of pipetting steps in the transfer-list "
        "format used by Beckman Biomek, Hamilton and Opentrons. Source plates hold the samples, "
        "destination plates hold the pools; both are addressed as 96-well plates.",
    )

    design_file = st.file_uploader("Pooling design (CSV generated by PoolPy)", type=["csv"],
                                   key="auto_design")
    if design_file is None:
        note("Volumes and plate names will very likely need adjusting for your deck layout. "
             "If your instrument is not covered here, get in touch and we will add it.")
        return

    try:
        wa_df = read_design_csv(design_file.getvalue())
    except Exception as exc:
        note(f"Could not read the file: {exc}", "bad")
        return

    ok, msg = validate_wa_df(wa_df)
    note(msg, "good" if ok else "bad")
    if not ok:
        return

    _render_protocols(wa_df.values)


@st.fragment
def _render_protocols(wa: np.ndarray):
    volume = st.number_input("Transfer volume per step (µL)", min_value=0.1, max_value=1000.0,
                             value=1.0, step=0.5,
                             help="Applied to every transfer. Adjust to your assay.")

    protocols = cached_robot_protocols(wa, volume)
    n_transfers = len(protocols["Biomek"])
    source_plates = int(np.ceil(wa.shape[0] / 96))
    dest_plates = int(np.ceil(wa.shape[1] / 96))

    stat_strip([
        ("Samples", f"{wa.shape[0]:,}", f"{source_plates} source plate(s)"),
        ("Pools", f"{wa.shape[1]:,}", f"{dest_plates} destination plate(s)"),
        ("Transfers", f"{n_transfers:,}", "individual pipetting steps"),
        ("Total volume", f"{n_transfers * volume:,.0f} µL", f"at {volume:g} µL per step"),
    ])

    section("Protocols")
    tabs = st.tabs(list(protocols))
    for tab, (name, df) in zip(tabs, protocols.items()):
        with tab:
            html_table(df.head(20).reset_index(drop=True), dense=True, wrap="cells",
                       caption=f"Showing the first 20 of {len(df):,} transfers.")
            st.download_button(f"Download {name} transfer list (CSV)",
                               deferred_csv(df, index=False),
                               file_name=f"{name.lower()}_protocol.csv", mime="text/csv",
                               key=f"robot_{name}")


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — GUIDE
# ═══════════════════════════════════════════════════════════════════════════

def render_guide():
    page_header(
        "Reference",
        "User guide",
        "What each pooling family does, how to read every column of the benchmark table, and "
        "what the decoder's verdicts mean.",
    )

    guide = user_guide_bytes()
    if guide is not None:
        st.download_button("Download the full user guide (PDF)", guide,
                           file_name="PoolPy_user_guide.pdf", mime="application/pdf")

    section("Positive-count dependent methods")
    st.markdown(
        "These change their design according to the number of expected positives. They are "
        "usually non-adaptive and suit high values of D.\n\n"
        "- **Random** — semi-adaptive; pools are formed by assigning samples at random.\n"
        "- **STD** — non-adaptive; built from prime numbers and modulus operations.\n"
        "- **Chinese remainder** — non-adaptive; built on the Chinese Remainder Theorem, with "
        "*backtrack* (prime powers minimising the pool count) and *special* (D = 2 or 3) variants.\n"
        "- **Hierarchical** — strictly adaptive; iteratively splits the positive pools. Its "
        "*N pools* entry is the list of splits per stage: `[3, 3]` means split into 3 pools, test, "
        "split each positive pool into 3 again, test, then test the remaining samples individually."
    )

    section("Positive-count independent methods")
    st.markdown(
        "These keep the same design whatever D is, so use them with care when more than one "
        "positive is expected.\n\n"
        "- **Matrix** — semi-adaptive; each sample sits in one row pool and one column pool.\n"
        "- **Multidimensional (3D, 4D, …)** — samples on a higher-dimensional grid, one pool per "
        "coordinate per axis.\n"
        "- **Binary** — samples are assigned to pools by a binary code, maximising information "
        "per test."
    )

    section("Benchmark table columns")
    st.markdown(
        "| Column | Meaning |\n|---|---|\n"
        "| Pooling strategy | Name of the method. |\n"
        "| Mean experiments | Average number of tests needed to identify the positives. |\n"
        "| Mean steps | Average number of successive testing rounds; 1 means fully non-adaptive. |\n"
        "| N pools (W) | Pools used in the first step. For Hierarchical, the list of splits. |\n"
        "| Max samples per pool | Largest number of samples combined in any one pool — the "
        "practical limit set by assay dilution. |\n"
        "| Percentage check | Share of cases needing a confirmation round beyond the first step. |\n"
        "| Mean extra experiments | Average number of tests beyond the first step. |\n"
        "| Max experiments per sample | Most times any individual sample is tested. |\n"
        "| Pooling sensitivity / specificity | Test sensitivity and specificity propagated "
        "through the pooling design, when you supply them. |"
    )

    section("Decoder verdicts")
    st.markdown(
        "- **Unique** — exactly one set of positive samples explains the readout.\n"
        "- **Multiple** — several sets of at most D samples explain the readout; all are listed.\n"
        "- **Putative** — more candidate sets than candidate samples, so only the pooled set of "
        "possible positives is reported. Confirm them individually.\n"
        "- **Error** — nothing is consistent with the readout. Check the input, or raise D."
    )

    section("Notation")
    st.markdown(
        "| Symbol | Meaning |\n|---|---|\n"
        "| **S** | Number of samples to test. |\n"
        "| **D** | Differentiate — maximum number of positives to resolve. |\n"
        "| **ρ** | Prevalence of positives in the population. |\n"
        "| **W** | Total number of pools a method needs. |"
    )


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE — ABOUT
# ═══════════════════════════════════════════════════════════════════════════

def render_about():
    page_header(
        "About",
        "PoolPy",
        "Combinatorial pooling designs that increase throughput and cut costs across biological "
        "experiments.",
    )

    st.markdown(
        "PoolPy replaces one-sample-one-test workflows with structured pooling designs that "
        "still identify positives while using far fewer measurements. It brings design "
        "selection, benchmarking, automation and decoding together in a single end-to-end "
        "platform.\n\n"
        "Ten conceptually different pooling algorithms are benchmarked across more than "
        "100,000 in-silico conditions, so designs can be chosen on the trade-offs that matter "
        "for a given application: test efficiency, pool size, number of steps, and robustness "
        "to dilution.\n\n"
        "A core idea is that **there is no universally best pooling strategy**. PoolPy tailors "
        "combinatorial testing to each use case by balancing throughput, prevalence, robustness "
        "and experimental constraints.\n\n"
        "Both standard binary screening assays and more complex multi-readout profiling "
        "experiments are supported. Across protein–ligand screening, RT-qPCR viral testing and "
        "genome-wide protein–DNA interaction profiling, PoolPy reduced the number of required "
        "measurements by roughly 60% to 93% while keeping results interpretable."
    )

    c1, c2, c3 = st.columns(3)
    c1.markdown(f'<div class="pp-card"><h4>Code</h4>'
                f'<p><a href="{GITHUB_URL}" target="_blank">github.com/trouillon-lab/PoolPy</a></p>'
                f'<p>Issues and contributions welcome.</p></div>', unsafe_allow_html=True)
    c2.markdown(f'<div class="pp-card"><h4>Paper</h4>'
                f'<p><a href="{ARXIV_URL}" target="_blank">PoolPy: Flexible Group Testing '
                f'Design for Large-Scale Screening</a></p><p>arXiv:2509.03481</p></div>',
                unsafe_allow_html=True)
    c3.markdown(f'<div class="pp-card"><h4>Data</h4>'
                f'<p><a href="{ZENODO_URL}" target="_blank">DOI: 10.5281/zenodo.18660061</a></p>'
                f'<p>All precomputed pooling designs.</p></div>', unsafe_allow_html=True)

    st.markdown(f'<div class="pp-lead" style="margin-top:1.4rem">For enquiries and feedback, '
                f'contact {CONTACT}.</div>', unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════
#  APP SHELL
# ═══════════════════════════════════════════════════════════════════════════

RENDERERS = {
    "Design": render_design,
    "Prevalence": render_prevalence,
    "Decoder": render_decoder,
    "Automation": render_automation,
    "Guide": render_guide,
    "About": render_about,
}

PAGE_BLURB = {
    "Design": "Compare methods and get designs",
    "Prevalence": "Choose D from a prevalence estimate",
    "Decoder": "Read out a pooled experiment",
    "Automation": "Generate robot protocols",
    "Guide": "Methods, metrics and notation",
    "About": "The project, paper and data",
}


def main():
    st.set_page_config(page_title="PoolPy", page_icon=favicon(),
                       layout="wide", initial_sidebar_state="expanded")
    st.markdown(CSS, unsafe_allow_html=True)

    with st.sidebar:
        logo = logo_png()
        if logo:
            st.image(logo, width="stretch")
        else:
            st.markdown('<div class="pp-title">PoolPy</div>', unsafe_allow_html=True)

        page = st.radio("Navigation", PAGES, label_visibility="collapsed",
                        captions=[PAGE_BLURB[p] for p in PAGES])

        st.markdown(
            f'<div class="pp-sidebar-links" style="margin-top:1.6rem">'
            f'<a href="{GITHUB_URL}" target="_blank">GitHub repository</a><br>'
            f'<a href="{ARXIV_URL}" target="_blank">Read the paper</a><br>'
            f'<a href="{ZENODO_URL}" target="_blank">Design archive (Zenodo)</a><br>'
            f'<span>{CONTACT}</span></div>', unsafe_allow_html=True)

    RENDERERS[page]()

    st.markdown(
        f'<div class="pp-foot">PoolPy · Trouillon lab, ETH Zürich · '
        f'<a href="{GITHUB_URL}" target="_blank" style="color:{BLUE_DARK}">source</a> · '
        f'<a href="{ARXIV_URL}" target="_blank" style="color:{BLUE_DARK}">arXiv:2509.03481</a>'
        f'</div>', unsafe_allow_html=True)


if __name__ == "__main__":
    main()
