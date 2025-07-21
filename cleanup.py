import os
import re

def extract_min_tests(filename):
    """Extract the number after '_ME_' and before '.csv'."""
    match = re.search(r'_ME_([\\d\\.]+)\\.csv$', filename)
    if match:
        return float(match.group(1))
    return None

def clean_wa_files(start_path):
    for root, dirs, files in os.walk(start_path):
        # Find files matching the pattern in this folder
        target_files = [
            f for f in files
            if f.startswith('WA_Random_N_') and '_diff_' in f and '_ME_' in f and f.endswith('.csv')
        ]
        if not target_files:
            continue

        # Map files to their min_tests value
        files_with_mt = []
        for fname in target_files:
            mt = extract_min_tests(fname)
            if mt is not None:
                files_with_mt.append((fname, mt))

        if not files_with_mt:
            continue

        # Find the smallest min_tests value
        smallest_mt = min(files_with_mt, key=lambda x: x[1])[1]

        # Remove all except those with the smallest min_tests
        for fname, mt in files_with_mt:
            if mt != smallest_mt:
                try:
                    os.remove(os.path.join(root, fname))
                    print(f"Removed: {os.path.join(root, fname)}")
                except Exception as e:
                    print(f"Could not remove {fname}: {e}")

# Example usage:
# clean_wa_files('/your/path/here')
