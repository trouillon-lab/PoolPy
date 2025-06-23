import os
import re
import sys

def delete_files_by_regex(root_folder, regex_pattern, safe_folder=None):
    pattern = re.compile(regex_pattern)
    for dirpath, dirnames, filenames in os.walk(root_folder):
        # Exclude the safe folder if specified
        if safe_folder:
            dirnames[:] = [d for d in dirnames if d != safe_folder]
            # Skip current dir if it's the safe folder
            if safe_folder in os.path.normpath(dirpath).split(os.sep):
                continue

        for filename in filenames:
            if pattern.match(filename):
                file_path = os.path.join(dirpath, filename)
                print(f"Deleting {file_path}")
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")

if __name__ == "__main__":
    # Command-line arguments handling
    if len(sys.argv) >= 3:  # At least root and regex provided
        root = sys.argv[1]
        regex = sys.argv[2]
        safe = sys.argv[3] if len(sys.argv) >= 4 else None
    else:
        # Interactive mode
        root = input("Enter the path to the root folder: ").strip()
        regex = input("Enter the regex pattern for filenames to delete: ").strip()
        safe_input = input("Enter safe folder name (optional): ").strip()
        safe = safe_input if safe_input else None

    delete_files_by_regex(root, regex, safe)



'''
# Delete files starting with 'decode_' everywhere
python script.py /path/to/root "^decode_.*"

# Delete files ending with '.tmp', but skip 'backups' folders
python script.py /path/to/root ".*\.tmp$" backups

# Interactive mode
python script.py
# (Follow prompts, leave 'safe folder' blank if you want no exclusions)

'''
