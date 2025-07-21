import os
import re
import glob

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

def replace_method_string_for_all_metrics(dpath):
    """Find all 'Metrics_N_*_diff_*.csv' files under dpath and replace 'method 1'."""
    # Regex pattern to match valid Metrics filenames
    metrics_filename_pattern = re.compile(r'^Metrics_N_\d+_diff_[\d\.]+\.csv$')

    for root, dirs, files in os.walk(dpath):
        for file in files:
            if metrics_filename_pattern.match(file):
                filepath = os.path.join(root, file)

                try:
                    with open(filepath, 'r', encoding='utf-8') as f:
                        content = f.read()

                    new_content = content.replace('method 1', 'First method')

                    if new_content != content:
                        with open(filepath, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"Replaced in: {filepath}")

                except Exception as e:
                    print(f"Failed to process {filepath}: {e}")

# === Usage ===
# Set your root paths and variables
dpath = '/your/main/folder'        # <-- Change this
WApath = os.path.join(dpath, 'WAs')
N = 10                             # <-- Set the value of N
diff = 0.1                         # <-- Set the value of diff

clean_wa_files(WApath)
replace_method_string_for_all_metrics(dpath)
