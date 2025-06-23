import os
import re

def delete_files_by_regex(root_folder, regex_pattern):
    pattern = re.compile(regex_pattern)
    for dirpath, dirnames, filenames in os.walk(root_folder):
        # Do not descend into any 'decoders' folders
        dirnames[:] = [d for d in dirnames if d != 'decoders']

        for filename in filenames:
            if pattern.match(filename):
                file_path = os.path.join(dirpath, filename)
                # Skip files already in any 'decoders' folder
                if 'decoders' in os.path.normpath(dirpath).split(os.sep):
                    continue
                print(f"Deleting {file_path}")
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")

if __name__ == "__main__":
    root = input("Enter the path to the root folder: ").strip()
    regex = input("Enter the regex pattern for filenames to delete: ").strip()
    delete_files_by_regex(root, regex)
