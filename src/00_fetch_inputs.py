"""
Step 0 - Download every input the pipeline needs.
Files already present are skipped
Sources are listed in config.py
"""

import sys
import requests
import config


def download_url(url, path):
    print(f"  downloading {path.name} from the original source ...")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)


def download_drive(file_id, path):
    """Download a file from Google Drive by its file ID
    """
    try:
        import gdown
    except ImportError:
        sys.exit("gdown is not installed. Run: pip install -r requirements.txt")
    print(f"  downloading {path.name} from Google Drive ...")
    gdown.download(id=file_id, output=str(path), quiet=False)


def main():
    downloaded = 0
    skipped = 0
    missing_ids = []

    print("Public sources:")
    for path, url in config.PUBLIC_URLS.items():
        if path.exists():
            print(f"  already here: {path.name}")
            skipped += 1
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        download_url(url, path)
        downloaded += 1

    print("\nGoogle Drive:")
    for path, file_id in config.DRIVE_IDS.items():
        if path.exists():
            print(f"  already here: {path.name}")
            skipped += 1
            continue
        if not file_id:
            print(f"  NO FILE ID    {path.name}")
            missing_ids.append(path)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        download_drive(file_id, path)
        downloaded += 1

    print(f"\nDownloaded {downloaded}, already present {skipped}.")

    if missing_ids:
        print(f"\n{len(missing_ids)} file(s) could not be fetched because no Drive ID "
              f"is set for them in config.py DRIVE_IDS:")
        for path in missing_ids:
            print(f"  {path.relative_to(config.ROOT)}")
        sys.exit(1)

    print("All inputs are in place. Next: python src/01_build_mutation_matrix.py")


if __name__ == "__main__":
    main()
