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
                    bktrk_fname = f'WA_Chinese_bktrk_N_{N_value}_diff_{diff_value}.csv'
                    np.savetxt(os.path.join(was_dir, bktrk_fname), WA_bktrk.astype(bool), delimiter=",")

                                        # Add row for Ch. Rm. Bktrk
                    bktrk_row = {
                        'Unnamed: 0': 'Ch. Rm. Bktrk',
                        'Method': 'Ch. Rm. Bktrk',
                        'Mean experiments': WA_bktrk.shape[1],
                        'Max compunds per well': int(np.max(np.sum(WA_bktrk, axis=0))),
                        'N wells': WA_bktrk.shape[1],
                        'Percentage check': 0,
                        'Mean extra experiments': 0,
                        'Mean steps': 1
                    }
                    if 'Method' in df.columns:
                        df = pd.concat([df, pd.DataFrame([bktrk_row])], ignore_index=True)

                    # If diff_value is 2 or 3, call assign_wells_chinese with special_diff=True
                    if diff_value in [2, 3]:
                        WA_special = assign_wells_chinese(n_compounds=N_value, differentiate=diff_value, special_diff=True)
                        special_fname = f'WA_Chinese_special_N_{N_value}_diff_{diff_value}.csv'
                        np.savetxt(os.path.join(was_dir, special_fname), WA_special.astype(bool), delimiter=",")
                        # Add row for Ch. Rm. Special
                        special_row = {
                            'Unnamed: 0': 'Ch. Rm. Special',
                            'Method': 'Ch. Rm. Special',
                            'Mean experiments': WA_special.shape[1],
                            'Max compunds per well': int(np.max(np.sum(WA_special, axis=0))),
                            'N wells': WA_special.shape[1],
                            'Percentage check': 0,
                            'Mean extra experiments': 0,
                            'Mean steps': 1
                        }
                        if 'Method' in df.columns:
                            df = pd.concat([df, pd.DataFrame([special_row])], ignore_index=True)

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

def sum_mean_experiments(row):
    if row['Method'] == 'Hierarchical':
        return row['Mean experiments']  # leave as is
    #print(type(row['Mean extra experiments']), type(row['N wells']))
    return row['Mean extra experiments'] + int(row['N wells'])



def adjust_mean_extra_experiments(df, N_value):
    """
    Ensure 'Mean extra experiments' <= N_value for all rows except 'Hierarchical' method.
    For 'Hierarchical', leave as is.
    """
    def adjust(row):
        if row['Method'] == 'Hierarchical':
            return row['Mean extra experiments']
        #print(row['Mean extra experiments'], N_value)
        return row['Mean extra experiments'] if row['Mean extra experiments'] <= N_value*row['Percentage check']/100 else N_value*row['Percentage check']/100
    df['Mean extra experiments'] = df.apply(adjust, axis=1)
    return df


def process_metrics_and_adjust_experiments(dpath):
    """
    Process all Metrics_N_*_diff_*.csv files:
    - Replace 'Chinese trick' → 'Chinese reminder'
    - Add rows for Ch. Rm. Bktrk and Ch. Rm. Special
    - After building the DataFrame, ensure 'Mean extra experiments' <= N_value.
      If not, set it to N_value.
    - Set 'Mean experiments' = 'Mean extra experiments' + 'N wells'
    - Drop duplicate Methods (keep min Mean experiments)
    """
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
                    with open(fpath, 'r', encoding='utf-8') as f:
                        content = f.read()

                    n_content = content.replace('Chinese trick', 'Chinese remainder')
                    new_content = n_content.replace('Chinese reminder', 'Chinese remainder')
                    df = pd.read_csv(StringIO(new_content))

                    was_dir = os.path.join(root, 'WAs')
                    if not os.path.exists(was_dir):
                        os.makedirs(was_dir)

                    WA_bktrk = assign_wells_chinese(n_compounds=N_value, differentiate=diff_value, backtrack=True)
                    bktrk_fname = f'WA_Chinese_bktrk_N_{N_value}_diff_{diff_value}.csv'
                    np.savetxt(os.path.join(was_dir, bktrk_fname), WA_bktrk.astype(bool), delimiter=",")

                    bktrk_row = {
                        'Unnamed: 0': 'Ch. Rm. Bktrk',
                        'Method': 'Ch. Rm. Bktrk',
                        'Mean experiments': WA_bktrk.shape[1],
                        'Max compunds per well': int(np.max(np.sum(WA_bktrk, axis=0))),
                        'N wells': WA_bktrk.shape[1],
                        'Percentage check': 0,
                        'Mean extra experiments': 0,
                        'Mean steps': 1
                    }
                    if 'Method' in df.columns:
                        df = pd.concat([df, pd.DataFrame([bktrk_row])], ignore_index=True)

                    if diff_value in [2, 3]:
                        WA_special = assign_wells_chinese(n_compounds=N_value, differentiate=diff_value, special_diff=True)
                        special_fname = f'WA_Chinese_special_N_{N_value}_diff_{diff_value}.csv'
                        np.savetxt(os.path.join(was_dir, special_fname), WA_special.astype(bool), delimiter=",")
                        special_row = {
                            'Unnamed: 0': 'Ch. Rm. Special',
                            'Method': 'Ch. Rm. Special',
                            'Mean experiments': WA_special.shape[1],
                            'Max compunds per well': int(np.max(np.sum(WA_special, axis=0))),
                            'N wells': WA_special.shape[1],
                            'Percentage check': 0,
                            'Mean extra experiments': 0,
                            'Mean steps': 1
                        }
                        if 'Method' in df.columns:
                            df = pd.concat([df, pd.DataFrame([special_row])], ignore_index=True)

                    # --- Adjust columns as requested ---
                    if 'Mean extra experiments' in df.columns and 'N wells' in df.columns:
                        df = adjust_mean_extra_experiments(df, N_value)
                        df['Mean experiments'] = df.apply(sum_mean_experiments, axis=1)

                        
                    for col in df.select_dtypes(include=['float', 'int']).columns:
                        df[col] = df[col].round(2)

                    if 'Method' in df.columns and 'Mean experiments' in df.columns:
                        df.sort_values('Mean experiments', inplace=True)
                        df = df.drop_duplicates(subset='Method', keep='first')
                        df.to_csv(fpath, index=False)
                        print(f"Processed: {fpath}")
                    else:
                        with open(fpath, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"Updated text only: {fpath}")

                except Exception as e:
                    print(f"Error processing {fpath}: {e}")


def remove_decoders_folders(root_path):
    """
    Recursively delete all folders named 'decoders' and their contents under root_path.
    """
    for current_root, dirs, files in os.walk(root_path, topdown=True):
        # Make a copy of dirs to avoid modifying while iterating
        for dir_name in list(dirs):
            if dir_name == 'decoders':
                decoders_path = os.path.join(current_root, dir_name)
                # Remove all files inside the decoders folder
                for root_dec, dirs_dec, files_dec in os.walk(decoders_path, topdown=False):
                    for file_dec in files_dec:
                        try:
                            os.remove(os.path.join(root_dec, file_dec))
                        except Exception as e:
                            print(f"Failed to remove file {file_dec}: {e}")
                    for dir_dec in dirs_dec:
                        try:
                            os.rmdir(os.path.join(root_dec, dir_dec))
                        except Exception as e:
                            print(f"Failed to remove subfolder {dir_dec}: {e}")
                # Remove the decoders folder itself
                try:
                    os.rmdir(decoders_path)
                    print(f"Removed folder: {decoders_path}")
                except Exception as e:
                    print(f"Failed to remove folder {decoders_path}: {e}")
                # Remove from dirs so os.walk doesn't visit it again
                dirs.remove(dir_name)
# === Usage ===
# Set your root paths and variables
dpath = 'D:\\precomputed'        # <-- Change this

clean_wa_files(dpath)
remove_decoders_folders(dpath)
#replace_method_filter_metrics_add_CT(dpath)
process_metrics_and_adjust_experiments(dpath)
