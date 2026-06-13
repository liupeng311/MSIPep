"""
Shared filtering logic for database search result post-processing.

- Decoy: recognize common target-decoy naming (DECOY, REV__, reverse DB, etc.) to remove pure decoy matches.
- Score thresholds are passed in by each engine script; this module provides tool-agnostic predicate functions.
"""

from __future__ import annotations

import re
from typing import List

# Strip modification tags: remove [xxx] mass offsets and non-letter characters
_MOD_STRIP = re.compile(r'\[.*?\]|[^A-Za-z]')
# Common MSFragger / Comet N/C-terminal modification prefix lowercase letters (n, c)
_TERM_MOD = re.compile(r'^[nc]', re.IGNORECASE)


# Split multi-protein fields by | ; ,
_PROTEIN_SPLIT = re.compile(r"[|;,]\s*")


def looks_like_decoy(accession: str) -> bool:
    """
    Whether a single accession / protein identifier should be treated as decoy.

    Covers: DECOY, decoy_, REV__, rev_, Reverse, shuffle, etc.;
    UniProt-style sp|DECOY_xxx| also matches within the string.
    """
    if not accession or not str(accession).strip():
        return False
    s = str(accession).strip()
    u = s.upper()

    if u.startswith("DECOY"):
        return True
    if u.startswith("REV__") or u.startswith("REV_"):
        return True
    if "SHUFFLE" in u or "SHUFFLED" in u:
        return True
    if u.startswith("REVERSE") or u.startswith("RANDOM"):
        return True
    # Decoy marker near start of string (e.g. sp|DECOY_...)
    if re.search(r"(?<![A-Z0-9])DECOY[_-]", u):
        return True
    return False


def split_protein_identifiers(protein_field: str) -> List[str]:
    if not protein_field or not str(protein_field).strip():
        return []
    parts = _PROTEIN_SPLIT.split(str(protein_field).strip())
    return [p.strip() for p in parts if p.strip()]


def target_proteins_only(protein_field: str) -> List[str]:
    """Return identifiers in the protein field judged as target (decoy parts removed)."""
    return [p for p in split_protein_identifiers(protein_field) if not looks_like_decoy(p)]


def row_has_target_hit(protein_field: str) -> bool:
    """Keep row if the PSM protein list has at least one target (otherwise discard entire row)."""
    return bool(target_proteins_only(protein_field))


def peptide_aa_length(peptide: str) -> int:
    """
    Sequence length for MHC-I peptide range checks.

    Handles common MSFragger / Comet modification formats, e.g.:
      M[147]PEPTIDE  →  8 aa
      n[43]ACDEFGHIK →  9 aa (N-terminal mod prefix n removed)
      C[160]PEPTIDER →  9 aa

    Steps:
      1. Remove [...] mass offset annotations
      2. Remove N/C-terminal modification lowercase prefix (n, c)
      3. Remove remaining non-letter characters
      4. Letter count is amino acid count
    """
    if not peptide:
        return 0
    seq = _MOD_STRIP.sub('', peptide)      # Remove [mod] and non-letters
    seq = _TERM_MOD.sub('', seq)           # Remove leading n/c terminal mod prefix
    return len(seq)


def in_class1_length_range(
    peptide: str, min_len: int, max_len: int
) -> bool:
    ln = peptide_aa_length(peptide.strip())
    return min_len <= ln <= max_len
