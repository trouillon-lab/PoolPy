import os
import re
import pandas as pd
from io import StringIO
from Functions import assign_wells_chinese
import numpy as np

def extract_min_tests(filename):
    """Extract float number after '_ME_' and before '.csv'."""
    match = re.search(r'_ME_([\d\.]+)\.csv$', filename)
    if match:
        return float(match.group(1))
    return None

def clean_wa_files(WApath):
    """In each folder under WApath, only keep the file with smallest _ME_ value."""
    for root, dirs, files in os.walk(WApath):
        wa_files = [f for f in files if f.startswith('WA_Random_N_') and '_ME_' in f and f.endswith('.csv')]
        if not wa_files:
            continue

        files_with_vals = [(f, extract_min_tests(f)) for f in wa_files if extract_min_tests(f) is not None]
        if not files_with_vals:
            continue

        # Get min ME value for this folder
        min_me = min(files_with_vals, key=lambda x: x[1])[1]

        # Remove all files with ME != min
        for fname, val in files_with_vals:
            if val != min_me:
                try:
                    os.remove(os.path.join(root, fname))
                    print("Removed:", os.path.join(root, fname))
                except Exception as e:
                    print("Failed to remove", fname, ":", e)


def replace_method_string_and_filter_metrics(dpath):
    """Process all Metrics_N_*_diff_*.csv files:
       - Replace 'method 1' → 'First method'
       - Drop duplicate Methods (keep min Mean experiments)
    """
    # File name pattern: Metrics_N_<num>_diff_<float>.csv
    metrics_filename_pattern = re.compile(r'^Metrics_N_\d+_diff_[\d\.]+\.csv$')

    for root, dirs, files in os.walk(dpath):
        for fname in files:
            if metrics_filename_pattern.match(fname):
                # Extract N and diff values from filename
                fpath = os.path.join(root, fname)
                match = re.match(r'^Metrics_N_(\d+)_diff_([\d\.]+)\.csv$', fname)
                if match:
                    N_value = int(match.group(1))
                    diff_value = float(match.group(2))
                else:
                    N_value = None
                    diff_value = None
                try:
                    # Load file as text
                    with open(fpath, 'r', encoding='utf-8') as f:
                        content = f.read()

                    # Replace 'method 1' → 'First method'
                    new_content = content.replace('method 1', 'First method')

                    # Try to read into DataFrame
                    df = pd.read_csv(StringIO(new_content))

                    if 'Method' in df.columns and 'Mean experiments' in df.columns:
                        # Drop duplicates keeping the one with the minimum 'Mean experiments'
                        df.sort_values('Mean experiments', inplace=True)
                        df = df.drop_duplicates(subset='Method', keep='first')

                        # Save back to CSV
                        df.to_csv(fpath, index=False)
                        print(f"Processed: {fpath} (N={N_value}, diff={diff_value})")
                    else:
                        # Save only the replaced version if columns not found
                        with open(fpath, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"Updated text only: {fpath} (N={N_value}, diff={diff_value})")

                except Exception as e:
                    print(f"Error processing {fpath}: {e}")



def replace_method_filter_metrics_add_CT(dpath):
    """Process all Metrics_N_*_diff_*.csv files:
       - Replace 'method 1' → 'First method'
       - Drop duplicate Methods (keep min Mean experiments)
    """
    # File name pattern: Metrics_N_<num>_diff_<float>.csv
    metrics_filename_pattern = re.compile(r'^Metrics_N_\d+_diff_[\d\.]+\.csv$')

    for root, dirs, files in os.walk(dpath):
        for fname in files:
            if metrics_filename_pattern.match(fname):
                fpath = os.path.join(root, fname)
                match = re.match(r'^Metrics_N_(\d+)_diff_([\d\.]+)\.csv$', fname)
                if match:
                    N_value = int(match.group(1))
                    diff_value = float(match.group(2))
                else:
                    N_value = None
                    diff_value = None
                try:
                    # Load file as text
                    with open(fpath, 'r', encoding='utf-8') as f:
                        content = f.read()

                    new_content = content.replace('Chinese trick', 'Chinese reminder')

                    # Try to read into DataFrame
                    df = pd.read_csv(StringIO(new_content))

                    # Prepare WAs subfolder path
                    was_dir = os.path.join(root, 'WAs')
                    if not os.path.exists(was_dir):
                        os.makedirs(was_dir)

                    # Call assign_wells_chinese with backtrack=True

                    WA_bktrk = assign_wells_chinese(n_compounds=N_value, differentiate=diff_value, backtrack=True)
                    bktrk_fname = f'WA_chinese_bktrk_N_{N_value}_diff_{diff_value}.csv'
                    np.savetxt(os.path.join(was_dir, bktrk_fname), WA_bktrk.astype(bool), delimiter=",")

                    # If diff_value is 2 or 3, call assign_wells_chinese with special_diff=True
                    if diff_value in [2, 3]:
                        WA_special = assign_wells_chinese(n_compounds=N_value, differentiate=diff_value, special_diff=True)
                        special_fname = f'WA_chinese_special_N_{N_value}_diff_{diff_value}.csv'
                        np.savetxt(os.path.join(was_dir, special_fname), WA_special.astype(bool), delimiter=",")

                    if 'Method' in df.columns and 'Mean experiments' in df.columns:
                        # Drop duplicates keeping the one with the minimum 'Mean experiments'
                        df.sort_values('Mean experiments', inplace=True)
                        df = df.drop_duplicates(subset='Method', keep='first')

                        # Save back to CSV
                        df.to_csv(fpath, index=False)
                        print(f"Processed: {fpath}")
                    else:
                        # Save only the replaced version if columns not found
                        with open(fpath, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"Updated text only: {fpath}")

                except Exception as e:
                    print(f"Error processing {fpath}: {e}")

# === Usage ===
# Set your root paths and variables
dpath = 'D:\precomputed\N_10'        # <-- Change this

clean_wa_files(dpath)
replace_method_string_and_filter_metrics(dpath)
