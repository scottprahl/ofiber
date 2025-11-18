"""Update CITATION.cff and README.rst citation block from latest GitHub release."""

from __future__ import annotations

import re
from pathlib import Path
import requests
import yaml

# Replace with your username and repo name
USERNAME = "scottprahl"
REPO = "pygrin"

GITHUB_API_URL = f"https://api.github.com/repos/{USERNAME}/{REPO}/releases/latest"

# Add a header so GitHub does not apply restrictive rate limits
HEADERS = {
    "User-Agent": f"{REPO}-citation-updater",
    "Accept": "application/vnd.github+json",
}

# --------------------------------------------------------------------
# Fetch latest release info from GitHub
# --------------------------------------------------------------------
response = requests.get(
    GITHUB_API_URL,
    timeout=10,
    headers=HEADERS,
)
response.raise_for_status()

release_info = response.json()
release_date = release_info["published_at"].split("T")[0]  # e.g. "2025-11-17"
tag_version = release_info["tag_name"]  # e.g. "v0.6.0" or "0.6.0"

# Normalize version: strip leading "v" if present
version = tag_version.lstrip("v")
year = release_date.split("-")[0]

# --------------------------------------------------------------------
# Update CITATION.cff
# --------------------------------------------------------------------
citation_path = Path("CITATION.cff")
if citation_path.exists():
    with citation_path.open("r", encoding="utf-8") as f:
        cff_data = yaml.safe_load(f)

    cff_changed = False

    if cff_data.get("date-released") != release_date:
        cff_data["date-released"] = release_date
        cff_changed = True

    if cff_data.get("version") != version:
        cff_data["version"] = version
        cff_changed = True

    if cff_changed:
        with citation_path.open("w", encoding="utf-8") as f:
            yaml.dump(cff_data, f, sort_keys=False)
        print(f"CITATION.cff updated → version: {version}, date: {release_date}")
    else:
        print("CITATION.cff: no change in release date or version.")
else:
    print("CITATION.cff not found; skipping CITATION.cff update.")

# --------------------------------------------------------------------
# Update citation block in README.rst
# --------------------------------------------------------------------
readme_path = Path("README.rst")
if not readme_path.exists():
    print("README.rst not found; skipping README citation update.")
else:
    text = readme_path.read_text(encoding="utf-8")
    original_text = text

    # 1. Prose citation line:
    #    Prahl, S. (2023). *pygrin: ...* (Version 0.5.1) [Computer software]. Zenodo. https://doi.org/...
    # Update the year in "Prahl, S. (YYYY)."
    text = re.sub(
        r"(Prahl,\s*S\.\s*\()(\d{4})(\)\.)",
        lambda m: f"{m.group(1)}{year}{m.group(3)}",
        text,
    )

    # Update "(Version X.Y.Z)" in the prose citation
    text = re.sub(
        r"\(Version [^)]+\)",
        f"(Version {version})",
        text,
    )

    # 2. BibTeX key:
    #    @software{pygrin_prahl_2023,
    text = re.sub(
        r"(@software\{pygrin_prahl_)(\d{4})(\s*,)",
        lambda m: f"{m.group(1)}{year}{m.group(3)}",
        text,
    )

    # 3. BibTeX year field:
    #    year      = {2023},
    text = re.sub(
        r"(year\s*=\s*\{)(\d{4})(\s*\},)",
        lambda m: f"{m.group(1)}{year}{m.group(3)}",
        text,
    )

    # 4. BibTeX version field:
    #    version   = {0.5.1},
    text = re.sub(
        r"(version\s*=\s*\{)([^}]+)(\s*\},)",
        lambda m: f"{m.group(1)}{version}{m.group(3)}",
        text,
    )

    if text != original_text:
        readme_path.write_text(text, encoding="utf-8")
        print(f"README.rst citation block updated → version: {version}, year: {year}")
    else:
        print("README.rst citation block already up to date.")
