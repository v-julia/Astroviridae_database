import os
import re
import calendar
import pandas as pd

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def load_mapping(mapping_file):
    """
    Load a TSV mapping file with columns: regex_pattern, short_code.
    Returns a list of (compiled_regex, short_code) tuples.
    """
    if not mapping_file or not os.path.exists(mapping_file):
        return []
    df = pd.read_csv(mapping_file, sep='\t', header=None, names=['pattern', 'code'])
    compiled = []
    for _, row in df.iterrows():
        try:
            pat = re.compile(str(row['pattern']), re.IGNORECASE)
            compiled.append((pat, row['code']))
        except re.error as e:
            print(f"Warning: invalid regex '{row['pattern']}': {e}")
    return compiled


def apply_mapping(value, compiled_patterns, default='NA'):
    """
    Apply regex patterns to a value and return the first matching short code.
    If no pattern matches, return a cleaned version of the original value.
    """
    if pd.isna(value):
        return default
    for pat, code in compiled_patterns:
        if pat.search(str(value)):
            return code
    # Fallback: clean the original value
    return str(value).replace(' ', '-').replace('/', '-')


def clean_field(value):
    """Replace problematic characters for safe use in filenames."""
    if pd.isna(value):
        return 'NA'
    return str(value).replace('/', '-').replace(' ', '-').replace(',', '-').replace('(', '-').replace(')', '-')


def standardize_date(date_str):
    """
    Convert various date formats to ISO 8601-like string:
      - YYYY-MM-DD if day available
      - YYYY-MM if only month available
      - YYYY if only year available
    Handles ranges (e.g., 08-Jul-2014/11-Oct-2016) by joining with '_'.
    """
    if pd.isna(date_str):
        return 'NA'
    s = str(date_str).strip()

    # Handle date range: split on '/' and standardize each part, then join with '_'
    if '/' in s and not re.match(r'^\d{4}$', s) and not re.match(r'^\d{4}-\d{2}$', s):
        parts = s.split('/')
        std_parts = [standardize_date(p) for p in parts]
        return '_'.join(std_parts)

    # DD-MMM-YYYY (e.g., 13-Jun-2024)
    match = re.match(r'^(\d{1,2})-([A-Za-z]{3})-(\d{4})$', s)
    if match:
        day, mon, year = match.groups()
        month_num = list(calendar.month_abbr).index(mon.capitalize())
        return f"{year}-{month_num:02d}-{int(day):02d}"

    # DD.MM.YYYY (e.g., 02.10.2019)
    match = re.match(r'^(\d{1,2})\.(\d{1,2})\.(\d{4})$', s)
    if match:
        day, mon, year = match.groups()
        return f"{year}-{int(mon):02d}-{int(day):02d}"

    # YYYY-MM (e.g., 2019-10)
    match = re.match(r'^(\d{4})-(\d{1,2})$', s)
    if match:
        year, mon = match.groups()
        return f"{year}-{int(mon):02d}"

    # MM-YYYY or MM.YYYY (e.g., 11-2025 or 11.2025)
    match = re.match(r'^(\d{1,2})[-.](\d{4})$', s)
    if match:
        mon, year = match.groups()
        return f"{year}-{int(mon):02d}"

    # MMM-YYYY (e.g., Nov-2025)
    match = re.match(r'^([A-Za-z]{3})-(\d{4})$', s)
    if match:
        mon, year = match.groups()
        month_num = list(calendar.month_abbr).index(mon.capitalize())
        return f"{year}-{month_num:02d}"

    # Just YYYY
    match = re.match(r'^(\d{4})$', s)
    if match:
        return s

    # If nothing matches, clean and return as is
    return clean_field(s)
