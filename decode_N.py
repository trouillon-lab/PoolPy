import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import os

from Fast_functions import *


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
