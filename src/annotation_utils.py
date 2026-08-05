import csv
import os

"""
Utility functions for mapping ORF names, reading CSV,
and handling keywords / exception lists.
"""


def read_exception_file(file_name):
    with open(file_name) as f:
        return [line.split(":")[0].strip() for line in f]

def read_csv(file_name, strip_it=True):
    """Read two‑column CSV mapping (base, new)."""
    if not os.path.exists(file_name):
        return {}
    with open(file_name) as f:
        reader = csv.DictReader(f, delimiter=",", fieldnames=["base", "new"])
        return {row["base"].strip(): row["new"].strip() for row in reader}

def map_feature(feature, feature_map):
    for k, v in feature_map.items():
        if feature.lower() == k.lower():
            return v
    return feature

def map_keywords(text, keywords_dict):
    codes = []
    for word, code in keywords_dict.items():
        if word.lower() in text.lower():
            codes.append(code)
    return list(set(codes))

