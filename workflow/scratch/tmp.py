import os
import argparse
import sys
from pathlib import Path

def rename_mixture_folders(
    root_dir: str,
    old_name: str = "mixture_k=None",
    new_name: str = "mixture_k=50",
    dry_run: bool = True
):
    """
    Recursively finds and renames subfolders.

    Args:
        root_dir (str): The starting directory to search from.
        old_name (str): The name of the folder to be renamed.
        new_name (str): The new name for the folder.
        dry_run (bool): If True, only prints the changes that would be made.
                        If False, performs the actual renaming.
    """
    print("--- Folder Renaming Script ---")
    if dry_run:
        print("[DRY RUN MODE] No changes will be made to the filesystem.")
    else:
        print("[LIVE MODE] Changes will be written to the filesystem.")
    print(f"Searching in: '{Path(root_dir).resolve()}'")
    print(f"Target folder name: '{old_name}'")
    print(f"New folder name: '{new_name}'")
    print("-" * 30)

    # Convert to Path objects for robust path handling
    root_path = Path(root_dir)
    if not root_path.is_dir():
        print(f"Error: The specified root directory does not exist: '{root_dir}'")
        sys.exit(1)

    folders_to_rename = []

    # First pass: find all folders to rename.
    # We walk from the bottom up to handle nested renames safely, e.g.,
    # /path/to/mixture_k=None/another/mixture_k=None
    for dirpath, dirnames, _ in os.walk(root_dir, topdown=False):
        if old_name in dirnames:
            # Construct the full path to the old and new folder names
            old_folder_path = Path(dirpath) / old_name
            new_folder_path = Path(dirpath) / new_name
            folders_to_rename.append((old_folder_path, new_folder_path))

    if not folders_to_rename:
        print("No folders found with the target name. Exiting.")
        return

    print(f"Found {len(folders_to_rename)} folder(s) to rename:")
    for old, new in folders_to_rename:
        if dry_run:
            print(f"  [DRY RUN] Would rename '{old}'\n                  TO '{new}'")
        else:
            print(f"  Attempting to rename '{old}' TO '{new}'")
            try:
                # Check if the new name already exists
                if new.exists():
                    print(f"  [SKIPPED] Cannot rename because destination already exists: '{new}'")
                    continue
                
                # Perform the rename
                os.rename(old, new)
                print(f"  [SUCCESS] Renamed successfully.")
            except OSError as e:
                print(f"  [ERROR] Failed to rename. Reason: {e}")

    print("-" * 30)
    print("Script finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Recursively find and rename subfolders named 'mixture_k=None' to 'mixture_k=50'.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "root_dir",
        type=str,
        help="The root directory where the search should begin."
    )
    
    parser.add_argument(
        "--live",
        action="store_true",
        help="Perform the actual renaming. If this flag is not set, the script will run in dry-run mode."
    )
    
    args = parser.parse_args()

    # The script runs in dry-run mode unless --live is specified
    is_dry_run = not args.live

    rename_mixture_folders(
        root_dir=args.root_dir,
        old_name="mixture_k=None",
        new_name="mixture_k=50",
        dry_run=is_dry_run
    )