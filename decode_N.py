import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import os

from Fast_functions import *



def _cell_is_valid_token(val, n_pools_local: int):
    try:
        # Skip NaN/None/empty
        if pd.isna(val):
            return None
        s = str(val).strip()
        if s == "" or s.lower() == "nan":
            return None
        # integer scalar
        try:
            iv = int(float(s))
            return 0 <= iv <= n_pools_local
        except Exception:
            pass
        # binary string of length n_pools
        if set(s).issubset({"0", "1"}) and len(s) == n_pools_local:
            return True
        # delimited list of ints
        parts = [p for p in re.split(r"[\s,;]+", s) if p]
        if parts and all(p.isdigit() and 0 <= int(p) <= n_pools_local for p in parts):
            return True
        return False
    except Exception:
        return False


def _parse_text_readout(readout_str: str):
    if not readout_str or not readout_str.strip():
        return None, "Error: Please enter a readout as comma-separated text."
    readout_list = []
    try:
        for x in readout_str.split(','):
            x_clean = x.strip()
            if x_clean == '':
                continue
            if not x_clean.isdigit():
                return None, f"Error: Invalid entry '{x_clean}' in readout. Please enter only integers separated by commas."
            readout_list.append(int(x_clean))
    except Exception:
        return None, "Error: Could not parse readout. Please enter only integers separated by commas."
    if len(readout_list) == 0:
        return None, "Error: Readout appears to be empty."
    return readout_list, None


def _decode_with_filter(readout_arr: np.ndarray, WA: np.ndarray, diff_deco: int):
    n_pools_local = WA.shape[1]
    n_compounds_local = WA.shape[0]

    if np.max(readout_arr) > 1 or len(readout_arr) != n_pools_local:
        readout_bin_ls = [1 if i in readout_arr else 0 for i in range(n_pools_local)]
        readout_use = np.array(readout_bin_ls)
    else:
        readout_use = readout_arr

    readout_bl = np.array(readout_use.astype(bool).astype(int))
    mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
    original_indices = np.where(mask)[0]
    filtered_WA = WA[mask]
    n_compounds_f = filtered_WA.shape[0]
    dval = min(diff_deco, n_compounds_f)

    if n_compounds_f < 2:
        if n_compounds_f == 1:
            decoded = [int(original_indices[0])]
            return "unique", decoded, decoded, None, n_compounds_f, None
        return "error", "No matches found. Check input or increase the differentiate value.", [], None, n_compounds_f, None

    ls_combs = [math.comb(n_compounds_f, i) for i in range(dval)]
    max_combs = np.sum(ls_combs)
    if max_combs > 1e4:
        decoded = [int(original_indices[j]) for j in range(len(original_indices))]
        decoded = sorted(decoded)
        return "putative", decoded, decoded, decoded, n_compounds_f, "max_combs"

    scrambler = {1: np.arange(n_compounds_f)}
    for j in range(2, dval + 1):
        scrambler[j] = np.array(list(itertools.combinations(np.arange(n_compounds_f), j)))

    decoded_pre = decode_precomp(
        well_assigner=filtered_WA,
        differentiate=dval,
        scrambler=scrambler,
        readout=readout_bl,
    )
    decoders = [list(c) if isinstance(c, (list, tuple, np.ndarray)) else [c] for c in decoded_pre]
    decoded = [[int(original_indices[k]) for k in comb] for comb in decoders]

    if len(decoded) == 0:
        return "error", "No matches found. Check input or increase the differentiate value.", decoded, None, n_compounds_f, None
    if len(decoded) == 1:
        return "unique", decoded[0], decoded, None, n_compounds_f, None
    if len(decoded) > n_compounds_f:
        decoded_set = sorted(list(set([x for combo in decoded for x in combo])))
        return "putative", decoded_set, decoded, decoded_set, n_compounds_f, "too_many_combos"
    return "multiple", decoded, decoded, None, n_compounds_f, None


def _build_text_message(file_name: str, readout_list: list, diff_deco: int, n_compounds: int, n_pools: int,
                        decoded_type: str, decoded: list, decoded_set: list, n_compounds_f: int, putative_reason: str):
    lines = []
    lines.append(f"Processing file {file_name} with max. {diff_deco} positive samples and readout:")
    lines.append(f"{readout_list}")
    lines.append(f"The uploaded pooling strategy comprizes {n_compounds} samples in {n_pools} pools.")

    if n_compounds_f < 2:
        if n_compounds_f == 1:
            lines.append("A single positive sample was found:")
            lines.append(f"Sample: {decoded[0]}")
        else:
            lines.append("We found no matches for the given parameters, check your input or try increasing the differentiate value.")
        return "\n".join(lines)

    if decoded_type == "putative" and decoded_set is not None and len(decoded_set) > 0:
        if putative_reason == "max_combs":
            lines.append("Putative positive samples were identified, but the app does not have the computational power to attempt to decode the exact combination.")
            lines.append("Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.")
        else:
            lines.append("Putative positive samples were identified, but the exact combination could not be pinpointed.")
            lines.append("Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.")
        lines.append(f"There are up to {min([diff_deco, len(decoded_set)])} positive samples among the following samples: {decoded_set}.")
        return "\n".join(lines)

    if decoded_type == "error":
        lines.append("We found no matches for the given parameters, check your input or try increasing the differentiate value.")
    elif decoded_type == "unique":
        if isinstance(decoded, list) and len(decoded) == 1 and isinstance(decoded[0], list) and len(decoded[0]) == 1:
            lines.append("A single positive sample was found:")
            lines.append(f"Sample: {decoded[0][0]}")
        else:
            lines.append("A single possible combination of positive samples was found. The positive samples are:")
            combo = decoded[0] if isinstance(decoded, list) and len(decoded) == 1 else decoded
            lines.append(f"Samples: {', '.join(map(str, combo))}")
    else:
        lines.append(f"{len(decoded)} possible combinations of positive samples were found. The possible combinations are:")
        for i, deco in enumerate(decoded):
            if i != 0:
                lines.append("or")
            lines.append(f"Samples: {deco}")

    return "\n".join(lines)




parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=-1, help='The number of maximum number of positives in the pooling strategy.')
parser.add_argument('--path_to_WA', type=str, help="A string argument containing the path to the pooling strategy file.")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout (readout in the form 0,1,0,1,0,0 or 3,6,14)")
parser.add_argument('--extensive_search', type=str, default='False', help="weather to search all possibilities (True) or use a faster exclusion based implementation (False, default)")

args = parser.parse_args()
args_dict = vars(args)

path_to_wa = args_dict['path_to_WA']
diff = args_dict['differentiate']
readout_in = args_dict['readout']

if not path_to_wa or not readout_in:
    raise SystemExit('Please provide both --path_to_WA and --readout.')

WA_df = pd.read_csv(path_to_wa, index_col=0)
WA = WA_df.values
n_pools = WA.shape[1]
n_compounds = WA.shape[0]

if diff == -1:
    diff = n_compounds

if readout_in.lower().endswith('csv'):
    readout_df = pd.read_csv(readout_in, header=None)
    if readout_df is None or readout_df.empty:
        raise SystemExit('Error: Readout CSV appears to be empty.')

    # Header handling: decide whether to include first row/first column as data
    if readout_df.shape[0] > 0:
        try:
            top = readout_df.iloc[0, :].tolist()
            validations = [_cell_is_valid_token(x, n_pools) for x in top]
            non_empty_validations = [v for v in validations if v is not None]
            if len(non_empty_validations) > 0 and not all(non_empty_validations):
                readout_df = readout_df.iloc[1:, :].reset_index(drop=True)
        except Exception:
            pass

    if readout_df.shape[1] > 0:
        try:
            first = readout_df.iloc[:, 0].tolist()
            validations = [_cell_is_valid_token(x, n_pools) for x in first]
            non_empty_validations = [v for v in validations if v is not None]
            if len(non_empty_validations) > 0 and not all(non_empty_validations):
                readout_df = readout_df.iloc[:, 1:].reset_index(drop=True)
        except Exception:
            pass

    results_rows = []
    for idx, row in readout_df.iterrows():
        rl, err = series_to_readout_list(row, n_pools)
        if rl is None:
            results_rows.append({
                "decoded_type": "error",
                "decoder_output": f"Error: {err}",
            })
            continue

        try:
            rl_sorted = sorted([int(x) for x in rl])
        except Exception:
            results_rows.append({
                "decoded_type": "error",
                "decoder_output": "Error: could not coerce values to integers",
            })
            continue

        readout_arr = np.array(rl_sorted, dtype=int)
        if len(readout_arr) == 0:
            results_rows.append({
                "decoded_type": "error",
                "decoder_output": "Error: empty readout after parsing",
            })
            continue

        try:
            if readout_arr.max() > n_pools:
                invalid_entries = [int(val) for val in rl_sorted if int(val) > n_pools]
                msg = (
                    f"Entries exceed number of pools ({n_pools}): {invalid_entries}. "
                    "Please correct your readout."
                )
                results_rows.append({
                    "decoded_type": "error",
                    "decoder_output": msg,
                })
                continue
        except Exception:
            pass

        decoded_type, decoder_output, _, _, _, _ = _decode_with_filter(readout_arr, WA, diff)
        results_rows.append({
            "decoded_type": decoded_type,
            "decoder_output": decoder_output,
        })

    df_out = pd.DataFrame(results_rows)
    if isinstance(df_out, pd.DataFrame) and not df_out.empty:
        df_with_index = df_out.copy()
        df_with_index.index = [f"Readout {i}" for i in range(len(df_with_index))]
        if "decoder_output" in df_with_index.columns:
            def format_output(val):
                if isinstance(val, list):
                    return ", ".join(map(str, val))
                if isinstance(val, str):
                    if val.startswith("[") and val.endswith("]"):
                        return val[1:-1]
                    return val
                return str(val)
            df_with_index["decoder_output"] = df_with_index["decoder_output"].apply(format_output)
        df_with_index.columns = [col.replace("_", " ").title() for col in df_with_index.columns]
    else:
        df_with_index = df_out

    decoded_csv = "decoded_readouts.csv"
    if isinstance(df_out, pd.DataFrame) and not df_out.empty:
        df_with_index.to_csv(decoded_csv, index=True)
    else:
        df_with_index.to_csv(decoded_csv, index=False)

else:
    readout_list, err = _parse_text_readout(readout_in)
    if err:
        raise SystemExit(err)

    readout_list.sort()
    readout_arr = np.array(readout_list, dtype=int)

    file_name = os.path.basename(path_to_wa)

    if max(readout_arr) > n_pools:
        invalid_entries = [val for val in readout_list if val > n_pools]
        msg = f"Processing file <b>{file_name}</b> with max. <b>{diff}</b> positive samples and readout: <br>"
        msg += f"{readout_list}<br>"
        msg += f"<br>The following entries in your readout are bigger than the number of pools ({n_pools}):<br> <b>{invalid_entries}</b>"
        msg += "<br>Please <b>correct your readout</b>.<br>"
        text_msg = _html_to_text(msg)
        decoded_txt = "decoder_output.txt"
        txt = f"Uploaded file: {file_name}\nReadout: {readout_in}\n\n{text_msg}"
        with open(decoded_txt, 'w+') as f:
            f.write(txt)
        print(text_msg)
        raise SystemExit(1)

    decoded_type, decoder_output, decoded, decoded_set, n_compounds_f, putative_reason = _decode_with_filter(readout_arr, WA, diff)
    text_msg = _build_text_message(file_name, readout_list, diff, n_compounds, n_pools, decoded_type, decoded, decoded_set, n_compounds_f, putative_reason)

    decoded_txt = "decoder_output.txt"
    txt = f"Uploaded file: {file_name}\nReadout: {readout_in}\n\n{text_msg}"
    with open(decoded_txt, 'w+') as f:
        f.write(txt)
    print(text_msg)
