import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import os
import json

from Fast_functions import decode_multi_readout_df, decode_single_readout_payload


parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=-1, help='The number of maximum number of positives in the pooling strategy.')
parser.add_argument('--path_to_WA', type=str, help="A string argument containing the path to the pooling strategy file.")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout (readout in the form 0,1,0,1,0,0 or 3,6,14)")
parser.add_argument('--extensive_search', type=str, default='False', help="weather to search all possibilities (True) or use a faster exclusion based implementation (False, default)")
parser.add_argument('--min_signal', type=float, default=None, help='Minimum signal threshold for constrained continuous decoding.')
parser.add_argument('--diluting', type=str, default='True', help='Whether to use dilution scaling in continuous decoding (True/False).')
parser.add_argument('--alpha', type=float, default=0.05, help='Elastic net regularization strength for continuous decoding.')
parser.add_argument('--l1_ratio', type=float, default=1, help='Elastic net l1 ratio for continuous decoding (0=ridge, 1=lasso).')
parser.add_argument('--force_continuous', type=str, default='False', help='Force decoder to treat readout as continuous (True/False).')
parser.add_argument('--grid_search', type=str, default='False', help='Grid search alpha and l1_ratio for continuous decoding (True/False).')
parser.add_argument('--return_all_fits', type=str, default='False', help='When grid_search=True for multi-readouts, include all grid-search fits in output (True/False).')
parser.add_argument('--save_grid_decoders', type=str, default='False', help='Save one decoder per grid element for each readout (requires grid_search=True).')
parser.add_argument('--grid_objective', type=str, default='mse', help='Objective used for grid search: mse, mae, rmse, r2.')
parser.add_argument('--decoded_output_base', type=str, default='decoded', help='Base name for saved decoded CSV outputs (suffixes are appended automatically).')

args = parser.parse_args()
args_dict = vars(args)

path_to_wa = args_dict['path_to_WA']
diff = args_dict['differentiate']
readout_in = args_dict['readout']
min_signal = args_dict['min_signal']
alpha = args_dict['alpha']
l1_ratio = args_dict['l1_ratio']
diluting_in = str(args_dict['diluting']).strip().lower()
diluting = diluting_in in ['1', 'true', 't', 'yes', 'y']
force_continuous_in = str(args_dict['force_continuous']).strip().lower()
force_continuous = force_continuous_in in ['1', 'true', 't', 'yes', 'y']
grid_search_in = str(args_dict['grid_search']).strip().lower()
grid_search = grid_search_in in ['1', 'true', 't', 'yes', 'y']
return_all_fits_in = str(args_dict['return_all_fits']).strip().lower()
return_all_fits = return_all_fits_in in ['1', 'true', 't', 'yes', 'y']
save_grid_decoders_in = str(args_dict['save_grid_decoders']).strip().lower()
save_grid_decoders = save_grid_decoders_in in ['1', 'true', 't', 'yes', 'y']
grid_objective = str(args_dict['grid_objective']).strip().lower()
decoded_output_base = str(args_dict['decoded_output_base']).strip() or 'decoded'
if grid_objective not in ['mse', 'mae', 'rmse', 'r2']:
    raise SystemExit("Error: --grid_objective must be one of: mse, mae, rmse, r2.")

if save_grid_decoders and not grid_search:
    print("Note: --save_grid_decoders was set but --grid_search is False. Disabling save_grid_decoders.")
    save_grid_decoders = False

if not path_to_wa or not readout_in:
    raise SystemExit('Please provide both --path_to_WA and --readout.')

WA_df = pd.read_csv(path_to_wa, index_col=0)
WA = WA_df.values
n_pools = WA.shape[1]
n_compounds = WA.shape[0]

if diff == -1:
    diff = n_compounds

if readout_in.lower().endswith('csv'):
    # Try reading with index_col=0 first (if CSV has row IDs)
    try:
        readout_df_with_index = pd.read_csv(readout_in, index_col=0, dtype=str)
        readout_ids = readout_df_with_index.index.tolist() if hasattr(readout_df_with_index, 'index') else None
        
        # Get data without index for processing
        readout_df = readout_df_with_index.reset_index(drop=True)
        readout_df = readout_df.astype(str)
    except (pd.errors.ParserError, pd.errors.EmptyDataError):
        # If that fails, try reading without index
        try:
            readout_df = pd.read_csv(readout_in, header=None, dtype=str)
            # First column should be the readout IDs
            if readout_df.shape[1] > 1:
                readout_ids = readout_df.iloc[:, 0].tolist()
                readout_df = readout_df.iloc[:, 1:].reset_index(drop=True)
            else:
                readout_ids = None
        except (pd.errors.ParserError, pd.errors.EmptyDataError):
            # Last resort: manually parse for inconsistent column counts or special delimiters
            readout_data = []
            readout_ids = []
            
            with open(readout_in, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    # First, check if line uses pipe delimiter (ID|data)
                    if '|' in line:
                        parts = line.split('|')
                        if len(parts) >= 2:
                            readout_ids.append(parts[0])
                            # Treat rest as comma-separated data
                            data_str = '|'.join(parts[1:])
                            remaining = [x.strip() for x in data_str.split(',') if x.strip()]
                            readout_data.append(remaining)
                            continue
                    
                    # Otherwise parse with quote handling for comma-separated
                    parts = []
                    current = ""
                    in_quotes = False
                    for char in line:
                        if char == '"':
                            in_quotes = not in_quotes
                        elif char == ',' and not in_quotes:
                            parts.append(current.strip(' "'))
                            current = ""
                        else:
                            current += char
                    if current:
                        parts.append(current.strip(' "'))
                    
                    if len(parts) > 0:
                        readout_ids.append(parts[0])
                        # For remaining parts, collect all commas into a single comma-separated string
                        remaining = parts[1:]
                        # If we have multiple remaining parts, they were split on commas - rejoin them
                        if len(remaining) > 1:
                            # Multiple parts mean they were comma-separated (e.g., 1,2,3)
                            remaining = [','.join(remaining)]
                        # Now split the comma-separated string into list of values
                        if remaining and remaining[0]:
                            values = remaining[0].split(',')
                            remaining = [x.strip() for x in values]
                        readout_data.append(remaining)
            
            # Create DataFrame from parsed data
            # Pad rows to same length if needed
            max_cols = max(len(row) for row in readout_data) if readout_data else 0
            for row in readout_data:
                while len(row) < max_cols:
                    row.append('')
            
            readout_df = pd.DataFrame(readout_data, dtype=str)
    
    if readout_df is None or readout_df.empty:
        raise SystemExit('Error: Readout CSV appears to be empty.')

    df_out = decode_multi_readout_df(
        readout_df,
        WA,
        diff,
        min_signal=min_signal,
        diluting=diluting,
        readout_ids=readout_ids,
        alpha=alpha,
        l1_ratio=l1_ratio,
        force_continuous=force_continuous,
        grid_search=grid_search,
        return_all_fits=return_all_fits,
        save_grid_decoders=save_grid_decoders,
        grid_objective=grid_objective,
    )

    if isinstance(df_out, pd.DataFrame) and not df_out.empty and save_grid_decoders and "grid_search_row_decoders" in df_out.columns:
        decoder_rows = []
        for i, row in df_out.iterrows():
            row_id = row.get("readout_id", f"Readout {i + 1}")
            row_decoders = row.get("grid_search_row_decoders", None)
            if isinstance(row_decoders, str):
                try:
                    row_decoders = json.loads(row_decoders)
                except Exception:
                    row_decoders = None
            if not isinstance(row_decoders, list):
                continue
            for fit in row_decoders:
                decoder_output = fit.get("decoder_output")
                if isinstance(decoder_output, np.ndarray):
                    decoder_output = decoder_output.tolist()
                decoder_rows.append({
                    "readout_id": row_id,
                    "alpha": fit.get("alpha"),
                    "l1_ratio": fit.get("l1_ratio"),
                    "alpha_scaled": fit.get("alpha_scaled"),
                    "objective": fit.get("objective", grid_objective),
                    "objective_value": fit.get("objective_value"),
                    "decoder_output": decoder_output,
                })

        if len(decoder_rows) > 0:
            pd.DataFrame(decoder_rows).to_csv(f"{decoded_output_base}_grid_decoders.csv", index=False)

    if isinstance(df_out, pd.DataFrame) and not df_out.empty:
        df_with_index = df_out.copy()
        if "grid_search_row_decoders" in df_with_index.columns:
            df_with_index = df_with_index.drop(columns=["grid_search_row_decoders"])
        if "grid_search_all_fits" in df_with_index.columns:
            df_with_index = df_with_index.drop(columns=["grid_search_all_fits"])
        if "readout_id" in df_with_index.columns:
            df_with_index = df_with_index.set_index("readout_id")
        else:
            df_with_index.index = [f"Readout {i + 1}" for i in range(len(df_with_index))]
        if "decoder_output" in df_with_index.columns:
            def format_output(val):
                if isinstance(val, np.ndarray):
                    # Format numpy array to single line with space-separated values
                    return " ".join([str(x) for x in val])
                if isinstance(val, list):
                    return " ".join([str(x) for x in val])
                if isinstance(val, str):
                    # Remove newlines and brackets
                    val = val.replace('\n', ' ').replace('  ', ' ').strip()
                    if val.startswith("[") and val.endswith("]"):
                        val = val[1:-1].strip()
                    return val
                return str(val)
            df_with_index["decoder_output"] = df_with_index["decoder_output"].apply(format_output)
        df_with_index.columns = [col.replace("_", " ").title() for col in df_with_index.columns]
    else:
        df_with_index = df_out

    # Keep standard decode output separate from grid-search decode output.
    decoded_csv = f"{decoded_output_base}_grid_optimum.csv" if grid_search else f"{decoded_output_base}_decoded_readouts.csv"
    if isinstance(df_out, pd.DataFrame) and not df_out.empty:
        df_with_index.to_csv(decoded_csv, index=True)
    else:
        df_with_index.to_csv(decoded_csv, index=False)

else:
    file_name = os.path.basename(path_to_wa)
    single_result = decode_single_readout_payload(
        readout_in,
        n_pools,
        WA,
        diff,
        min_signal=min_signal,
        diluting=diluting,
        alpha=alpha,
        l1_ratio=l1_ratio,
        force_continuous=force_continuous,
        grid_search=grid_search,
        grid_objective=grid_objective,
    )

    if single_result.get("decoded_type") == "error":
        text_msg = str(single_result.get("decoder_output"))
    else:
        text_msg = (
            f"Processing file {file_name} with max. {diff} positive samples.\n"
            f"Readout: {readout_in}\n"
            f"Decoded type: {single_result.get('decoded_type')}\n"
            f"Decoder output: {single_result.get('decoder_output')}"
        )

    decoded_txt = "decoder_output.txt"
    txt = f"Uploaded file: {file_name}\nReadout: {readout_in}\n\n{text_msg}"
    with open(decoded_txt, 'w+') as f:
        f.write(txt)
    print(text_msg)
