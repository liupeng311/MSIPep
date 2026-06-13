import argparse
import csv
import io
import os
import sys
from pathlib import Path

_MS_DIR = Path(__file__).resolve().parent
if str(_MS_DIR) not in sys.path:
    sys.path.insert(0, str(_MS_DIR))

from peptide_filters import (
    in_class1_length_range,
    row_has_target_hit,
    target_proteins_only,
)


def _norm_key(k: str) -> str:
    return (k or "").strip().lower().replace(" ", "_")


def _get_float(row: dict, *keys, default=None):
    for key in keys:
        v = row.get(key)
        if v is None or v == "":
            continue
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return default


def _normalize_reader_row(row: dict) -> dict:
    return {_norm_key(k): v for k, v in row.items()}


def _comet_txt_dict_reader(filepath: str):
    """
    Comet PEP export often has metadata lines before the table; the real header starts with scan\\t.
    Returns a csv.DictReader for that table.
    """
    with open(filepath, "r", encoding="utf-8") as infile:
        while True:
            line = infile.readline()
            if not line:
                raise ValueError(
                    f"{os.path.basename(filepath)}: Comet header not found (expected first column scan)"
                )
            if line.startswith("\ufeff"):
                line = line[1:]
            parts = line.split("\t")
            first = parts[0].strip().lower() if parts else ""
            if first == "scan":
                tbl = line + infile.read()
                return csv.DictReader(
                    io.StringIO(tbl),
                    delimiter="\t",
                    skipinitialspace=True,
                )
            # Skip banner lines until scan header


def comet_to_fasta(
    input_file,
    output_fasta,
    output_combined_fasta,
    xcorr_threshold=2.0,
    evalue_threshold=0.01,
    min_deltacn=0.05,
    min_length=8,
    max_length=11,
    exclude_decoy=True,
    append=False,
    strict_deltacn=False,
    seq_count_start=1,
):
    """
    Comet PSM export.

    Default thresholds: XCorr≥2.0, E-value≤0.01; ΔCn≥0.05 is applied only when the field exists and has a valid value
    (missing ΔCn does not discard the row unless strict_deltacn=True).

    seq_count_start: starting index for Class I FASTA headers; passed across files by caller to avoid duplicate headers.
    Returns: next available index after this file (for caller to pass to the next file).
    """
    mode = "a" if append else "w"
    try:
        reader = _comet_txt_dict_reader(input_file)

        seq_count = seq_count_start

        n_psm = 0
        n_main = 0
        n_cls = 0

        with open(output_fasta, mode, encoding="utf-8") as main_fasta, open(
            output_combined_fasta, mode, encoding="utf-8"
        ) as combined_fasta:

            for row_raw in reader:
                n_psm += 1
                try:
                    row = _normalize_reader_row(row_raw)
                    peptide = (
                        row.get("plain_peptide")
                        or row.get("peptide")
                        or row.get("plainpeptide")
                    )
                    protein_raw = (
                        row.get("protein")
                        or row.get("protein_id")
                        or row.get("proteins")
                    )
                    if not peptide or not protein_raw:
                        continue

                    peptide = str(peptide).strip()
                    if exclude_decoy and not row_has_target_hit(protein_raw):
                        continue

                    xcorr = _get_float(row, "xcorr", "x_corr", default=0.0)
                    evalue = _get_float(
                        row,
                        "e-value",
                        "e_value",
                        "expect",
                        "expectation",
                        "spectral_e-value",
                        default=9.0,
                    )

                    if xcorr < xcorr_threshold or evalue > evalue_threshold:
                        continue

                    dcn = _get_float(
                        row,
                        "delta_cn",
                        "deltacn",
                        "delta_c",
                        "deltacorr",
                        "delta_corr",
                        default=None,
                    )
                    if min_deltacn is not None and min_deltacn > 0:
                        if dcn is not None:
                            if dcn < min_deltacn:
                                continue
                        elif strict_deltacn:
                            continue

                    if exclude_decoy:
                        tprot = target_proteins_only(protein_raw)
                        if not tprot:
                            continue
                        label_parts = sorted(set(tprot))
                    else:
                        from peptide_filters import split_protein_identifiers
                        label_parts = sorted(set(split_protein_identifiers(protein_raw))) or [protein_raw]

                    main_label = label_parts[0] if len(label_parts) == 1 else "|".join(label_parts)

                    main_fasta.write(f">{main_label}\n{peptide}\n")
                    n_main += 1

                    if in_class1_length_range(
                        peptide, min_length, max_length
                    ):
                        combined_fasta.write(
                            f">Sequence_{seq_count}\n{peptide}\n"
                        )
                        seq_count += 1
                        n_cls += 1
                except Exception:
                    continue

        extra = ""
        if n_psm > 0 and n_main == 0:
            extra = (
                " | Warning: rows>0 but no peptides written: check column names (plain_peptide/peptide, xcorr, "
                "expect/e-value, delta_cn, protein) and thresholds match your Comet version."
            )
        elif n_main > 0 and n_cls == 0:
            extra = (
                f" | Note: no class I peptides in {min_length}-{max_length} aa range."
            )
        print(
            f"✅ Processed: {os.path.basename(input_file)} | "
            f"PSM rows {n_psm}; comet_all {n_main}; comet_classI {n_cls}{extra}"
        )
        return seq_count   # Return next available index
    except Exception as e:
        print(f"❌ Failed to process {input_file}: {e}")
        return seq_count_start  # On error, return original start index so later file numbering is unaffected


def main():
    p = argparse.ArgumentParser(
        description="Comet result post-processing: merge multiple files to FASTA; decoy removal, XCorr/E-value/ΔCn and MHC-I length filtering."
    )
    p.add_argument("--input_files", nargs="+", required=True)
    p.add_argument("--output_fasta", required=True)
    p.add_argument("--output_combined_fasta", required=True)
    p.add_argument(
        "--xcorr_threshold",
        type=float,
        default=2.0,
        help="Minimum XCorr (typical peptide ID 1.5–2.5, default 2.0)",
    )
    p.add_argument(
        "--evalue_threshold",
        type=float,
        default=0.01,
        help="Maximum E-value / expect (default 0.01)",
    )
    p.add_argument(
        "--min_deltacn",
        type=float,
        default=0.05,
        help="Minimum ΔCn; separates top from second-best hit (typical 0.05–0.1). 0 disables ΔCn filtering.",
    )
    p.add_argument("--min_length", type=int, default=8, help="Class I minimum aa length (default 8)")
    p.add_argument("--max_length", type=int, default=11, help="Class I maximum aa length (default 11)")
    p.add_argument(
        "--keep_decoy",
        action="store_true",
        help="Keep decoy protein matches (default: remove typical DECOY/REV identifiers)",
    )
    p.add_argument(
        "--strict_deltacn",
        action="store_true",
        help="Discard when ΔCn column is missing or row has no value (default: discard only when ΔCn is present and below threshold)",
    )

    args = p.parse_args()
    exclude_decoy = not args.keep_decoy
    min_dcn = args.min_deltacn if args.min_deltacn > 0 else None

    if os.path.exists(args.output_fasta):
        os.remove(args.output_fasta)
    if os.path.exists(args.output_combined_fasta):
        os.remove(args.output_combined_fasta)

    print(f"Processing {len(args.input_files)} Comet file(s)...")
    if exclude_decoy:
        print("  Filters: remove decoy; XCorr/E-value; ΔCn (strict only when column present); Class I length.")
    else:
        print("  Warning: keeping decoy matches.")

    seq_count_start = 1   # Accumulate across files so Class I header indices are globally unique
    for i, infile in enumerate(args.input_files):
        seq_count_start = comet_to_fasta(
            infile,
            args.output_fasta,
            args.output_combined_fasta,
            xcorr_threshold=args.xcorr_threshold,
            evalue_threshold=args.evalue_threshold,
            min_deltacn=min_dcn,
            min_length=args.min_length,
            max_length=args.max_length,
            exclude_decoy=exclude_decoy,
            append=(i > 0),
            strict_deltacn=args.strict_deltacn,
            seq_count_start=seq_count_start,
        )

    print("🎉 Comet post-processing complete!")
    print(f"   → {args.output_fasta}")
    print(f"   → {args.output_combined_fasta}")


if __name__ == "__main__":
    main()
