import os

import sys

import shutil

import argparse

import subprocess

from collections import OrderedDict

from pathlib import Path



# Post-processing scripts (resolved relative to this file to avoid cwd issues)

_SCRIPT_DIR = Path(__file__).resolve().parent

COMET_PEP_SCRIPT = _SCRIPT_DIR / "comet-pep.py"

MSFRAGGER_PEP_SCRIPT = _SCRIPT_DIR / "msfragger-pep-I.py"



def merge_fasta_files(fasta_list, output_path):

    with open(output_path, 'w', encoding='utf-8') as outfile:

        for fasta in fasta_list:

            if not os.path.exists(fasta):

                print(f"⚠️ FASTA file not found: {fasta}")

                continue

            with open(fasta, 'r', encoding='utf-8', errors='replace') as infile:

                shutil.copyfileobj(infile, outfile)

    print(f"✅ FASTA merge complete: {output_path}")



def _iter_fasta_records(path: str):

    """Parse FASTA records one at a time."""

    header = None

    seq_chunks = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:

        for line in f:

            line = line.rstrip("\n\r")

            if not line:

                continue

            if line.startswith(">"):

                if header is not None:

                    yield header, "".join(seq_chunks).replace(" ", "")

                header = line[1:].strip()

                seq_chunks = []

            else:

                seq_chunks.append(line.strip())

        if header is not None:

            yield header, "".join(seq_chunks).replace(" ", "")



def merge_classi_fastas_dedupe(input_paths, output_path, prefer_first_header: bool = False):

    """Merge multiple Class I FASTAs with peptide sequence deduplication; join headers with | by default."""

    seq_to_headers: "OrderedDict[str, list[str]]" = OrderedDict()



    def add_record(h: str, seq: str):

        if not seq:

            return

        key = seq.upper()

        if prefer_first_header:

            if key not in seq_to_headers:

                seq_to_headers[key] = [h]

            return

        if key not in seq_to_headers:

            seq_to_headers[key] = []

        if h not in set(seq_to_headers[key]):

            seq_to_headers[key].append(h)



    for path in input_paths:

        if not os.path.isfile(path):

            print(f"⚠️ Skipping missing file during Class I merge: {path}")

            continue

        for h, seq in _iter_fasta_records(path):

            add_record(h, seq)



    n = 0

    with open(output_path, "w", encoding="utf-8", newline="\n") as out:

        for seq_key, headers in seq_to_headers.items():

            if prefer_first_header:

                title = headers[0] if headers else seq_key[:20]

            else:

                seen_t = set()

                uniq = []

                for t in headers:

                    if t and t not in seen_t:

                        seen_t.add(t)

                        uniq.append(t)

                title = "|".join(uniq) if uniq else seq_key[:20]

            out.write(f">{title}\n{seq_key}\n")

            n += 1

    print(f"✅ Class I deduplicated merge wrote {n} peptides → {output_path}")





def _dedupe_abs_dirs(paths):

    """Convert paths to directories and deduplicate while preserving search order."""

    seen = set()

    dirs = []

    for p in paths:

        if not p:

            continue

        p = os.path.abspath(p)

        d = os.path.dirname(p) if os.path.isfile(p) else p

        if not os.path.isdir(d):

            continue

        if d not in seen:

            seen.add(d)

            dirs.append(d)

    return dirs





def _find_named_file(filename: str, directories):

    for d in directories:

        cand = os.path.join(d, filename)

        if os.path.isfile(cand):

            return cand

    return None





def update_param_file(param_path, db_path, new_param_path):

    """Update database path in parameter file."""

    if not os.path.exists(param_path):

        raise FileNotFoundError(f"Parameter template not found: {param_path}")



    with open(param_path, 'r', encoding='utf-8') as f:

        lines = f.readlines()



    updated = False

    with open(new_param_path, 'w', encoding='utf-8') as f:

        for line in lines:

            if line.strip().startswith('database_name'):

                comment = "  " + line[line.find('#'):] if '#' in line else ""

                f.write(f"database_name = {db_path}{comment}\n")

                updated = True

            else:

                f.write(line)



    if not updated:

        with open(new_param_path, 'a', encoding='utf-8') as f:

            f.write(f"\ndatabase_name = {db_path}\n")



    print(f"✅ Parameter file updated: {new_param_path}")





def run_msfragger(jar_path, param_file, raw_files, java_memory="64G"):

    for raw in raw_files:

        if not os.path.exists(raw):

            raise FileNotFoundError(f"RAW file not found: {raw}")



    raw_str = " ".join(raw_files)

    cmd = f"java -Xmx{java_memory} -jar {jar_path} {param_file} {raw_str}"

    print(f"🚀 Running MSFragger ({len(raw_files)} RAW file(s))")

    subprocess.run(cmd, shell=True, check=True)



def run_comet(comet_exe, param_file, mgf_file, output_dir):

    """Run Comet using absolute paths."""

    if not os.path.exists(comet_exe):

        raise FileNotFoundError(f"Comet executable not found: {comet_exe}")

    if not os.path.exists(mgf_file):

        raise FileNotFoundError(f"MGF file not found: {mgf_file}")



    abs_exe = os.path.abspath(comet_exe)

    abs_param = os.path.abspath(param_file)

    abs_mgf = os.path.abspath(mgf_file)



    cmd = f"{abs_exe} -P{abs_param} {abs_mgf}"

    print(f"🚀 Running Comet: {cmd}")

    subprocess.run(cmd, shell=True, check=True, cwd=output_dir)



def main():

    parser = argparse.ArgumentParser(description="MS database search workflow — merge DB, search, post-process, Class I merge")



    parser.add_argument(

        "--fasta_list",

        nargs="+",

        default=None,

        help="One or more search FASTA paths (concatenated in order as combined_protein_db.fasta)",

    )

    parser.add_argument("--output_dir", required=True, help="Output directory")



    parser.add_argument("--raw", nargs='+', help="One or more .raw files (for MSFragger)")

    parser.add_argument("--mgf", help="MGF file path (for Comet)")



    parser.add_argument("--java_memory", default="64G", help="Java memory limit")

    parser.add_argument(

        "--msfragger_dir",

        default="./software/MSFragger-3.8",

        help="MSFragger install directory (contains jar and closed_fragger.params)",

    )

    parser.add_argument(

        "--comet_dir",

        default="./software/comet",

        help="Comet config directory (contains comet.params template)",

    )

    parser.add_argument(

        "--comet_exe",

        default="/data/liup/MSIPep/software/comet/comet.2021010.linux.exe",

        help="Absolute path to Comet executable",

    )

    parser.add_argument(

        "--search-only",

        action="store_true",

        help="Merge DB + search only; skip post-processing and Class I merge",

    )

    parser.add_argument(

        "--postprocess-only",

        action="store_true",

        help="Skip search; collect existing results and post-process (requires search output in output_dir)",

    )

    parser.add_argument(

        "--keep_decoy",

        action="store_true",

        help="Keep decoy matches in post-processing (default: remove; for troubleshooting only)",

    )

    parser.add_argument(

        "--comet_strict_deltacn",

        action="store_true",

        help="Comet post-processing: discard rows when ΔCn is missing or has no value (default: discard only when ΔCn is present and below threshold)",

    )

    parser.add_argument(

        "--no-merge-class1",

        action="store_true",

        help="Skip generating merged_classI_dedup.fasta (merge runs by default even if only one Class I source exists)",

    )

    parser.add_argument(

        "--class1-merge-headers",

        choices=("join", "first"),

        default="join",

        help="Class I merge: join multiple headers for same peptide with |; first keeps first header only",

    )

    # MSFragger post-processing thresholds

    parser.add_argument("--msf_min_length", type=int, default=8)

    parser.add_argument("--msf_max_length", type=int, default=11)

    parser.add_argument("--msf_max_expectation", type=float, default=0.01)

    parser.add_argument("--msf_min_hyperscore", type=float, default=15)

    # Comet post-processing thresholds

    parser.add_argument("--comet_xcorr_threshold", type=float, default=2.0)

    parser.add_argument("--comet_evalue_threshold", type=float, default=0.01)

    parser.add_argument("--comet_min_deltacn", type=float, default=0.05)

    parser.add_argument("--comet_min_length", type=int, default=8)

    parser.add_argument("--comet_max_length", type=int, default=11)



    args = parser.parse_args()



    if not args.postprocess_only and not args.raw and not args.mgf:

        parser.error("Must provide --raw or --mgf (or use --postprocess-only)")



    if not args.postprocess_only and not args.fasta_list:

        parser.error("Non postprocess-only mode requires --fasta_list")



    output_dir = os.path.abspath(args.output_dir)



    # Validate all inputs before creating output directory (avoid empty dirs on failure)

    fasta_paths = [os.path.abspath(p) for p in (args.fasta_list or [])]

    if fasta_paths:

        missing = [p for p in fasta_paths if not os.path.isfile(p)]

        if missing:

            raise FileNotFoundError(f"The following FASTA files do not exist or are not files: {missing}")



    os.makedirs(output_dir, exist_ok=True)



    merged_fasta = os.path.join(output_dir, "combined_protein_db.fasta")

    msfragger_dir = os.path.abspath(args.msfragger_dir)

    comet_dir = os.path.abspath(args.comet_dir)



    if not args.postprocess_only:

        # 1. Merge FASTAs (all DB files from --fasta_list)

        merge_fasta_files(fasta_paths, merged_fasta)



        # 2. Generate parameter files

        msfragger_params = os.path.join(output_dir, "msfragger_run.params")

        comet_params = os.path.join(output_dir, "comet_run.params")



        update_param_file(

            os.path.join(msfragger_dir, "closed_fragger.params"), merged_fasta, msfragger_params

        )

        update_param_file(os.path.join(comet_dir, "comet.params"), merged_fasta, comet_params)



        # 3. Run search

        if args.raw:

            run_msfragger(

                jar_path=os.path.join(msfragger_dir, "MSFragger-3.8.jar"),

                param_file=msfragger_params,

                raw_files=args.raw,

                java_memory=args.java_memory,

            )



        if args.mgf:

            run_comet(

                comet_exe=args.comet_exe,

                param_file=comet_params,

                mgf_file=args.mgf,

                output_dir=output_dir,

            )

    else:

        print("⏭ Skipping search (--postprocess-only)")



    if args.search_only:

        print(f"\n🎉 Search complete (postprocess skipped)! Results directory: {output_dir}")

        return



    # ==================== 4. Collect results by expected filenames (avoid glob false positives) ====================

    print("\n🔄 Collecting search result files...")



    fragger_tsv_dir = os.path.join(output_dir, "fragger_tsv")

    os.makedirs(fragger_tsv_dir, exist_ok=True)

    tsv_moved = 0



    # MSFragger: output is typically <RAW_basename>.tsv (e.g. .raw → .tsv)

    if args.raw:

        fragger_search_roots = _dedupe_abs_dirs(

            [os.getcwd(), output_dir] + list(args.raw)

        )

        for raw_path in args.raw:

            tsv_name = os.path.splitext(os.path.basename(raw_path))[0] + ".tsv"

            dest = os.path.join(fragger_tsv_dir, tsv_name)

            if os.path.isfile(dest):

                tsv_moved += 1

                print(f"  ✓ Already exists (skip move): {dest}")

                continue

            src = _find_named_file(tsv_name, fragger_search_roots)

            if src:

                shutil.move(src, dest)

                tsv_moved += 1

                print(f"  ✓ MSFragger: {src} → {dest}")

            else:

                print(

                    f"  ⚠️ {tsv_name} not found (searched in: {fragger_search_roots})"

                )

    elif args.postprocess_only:

        import glob as _glob

        existing = _glob.glob(os.path.join(fragger_tsv_dir, "*.tsv"))

        tsv_moved = len(existing)

        if tsv_moved:

            print(f"  ✓ postprocess-only: using {tsv_moved} existing MSFragger TSV file(s)")



    print(f"✅ Prepared {tsv_moved} MSFragger .tsv file(s) (in fragger_tsv/)")



    # Comet: output is typically <MGF_basename>.txt

    comet_txt_files = []

    if args.mgf:

        txt_name = os.path.splitext(os.path.basename(args.mgf))[0] + ".txt"

        dest_txt = os.path.join(output_dir, txt_name)

        comet_search_roots = _dedupe_abs_dirs(

            [output_dir, os.getcwd(), args.mgf]

        )

        if os.path.isfile(dest_txt):

            comet_txt_files.append(dest_txt)

            print(f"  ✓ Comet result already in output directory: {dest_txt}")

        else:

            src = _find_named_file(txt_name, comet_search_roots)

            if src:

                if os.path.abspath(src) != os.path.abspath(dest_txt):

                    shutil.move(src, dest_txt)

                comet_txt_files.append(dest_txt)

                print(f"  ✓ Comet: {src} → {dest_txt}")

            else:

                print(

                    f"  ⚠️ {txt_name} not found (searched in: {comet_search_roots})"

                )

    elif args.postprocess_only:

        import glob as _glob

        comet_txt_files = sorted(_glob.glob(os.path.join(output_dir, "*.txt")))

        if comet_txt_files:

            print(f"  ✓ postprocess-only: using {len(comet_txt_files)} existing Comet TXT file(s)")



    print(f"✅ Comet post-processing will use {len(comet_txt_files)} .txt file(s)")



    # 5. MSFragger post-processing

    if tsv_moved > 0:

        msfragger_post = [

            sys.executable,

            str(MSFRAGGER_PEP_SCRIPT),

            "-i",

            fragger_tsv_dir,

            "-o",

            os.path.join(output_dir, "fragger_classI.fasta"),

            "--min_length",

            str(args.msf_min_length),

            "--max_length",

            str(args.msf_max_length),

            "--max_expectation",

            str(args.msf_max_expectation),

            "--min_hyperscore",

            str(args.msf_min_hyperscore),

        ]

        if args.keep_decoy:

            msfragger_post.append("--keep_decoy")

        print("🔄 Running MSFragger Class I peptide extraction (quality and decoy filtering enabled)...")

        subprocess.run(msfragger_post, check=True)



    # 6. Comet post-processing

    if comet_txt_files:

        comet_post = [

            sys.executable,

            str(COMET_PEP_SCRIPT),

            "--input_files",

            *comet_txt_files,

            "--output_fasta",

            os.path.join(output_dir, "comet_all.fasta"),

            "--output_combined_fasta",

            os.path.join(output_dir, "comet_classI.fasta"),

            "--xcorr_threshold",

            str(args.comet_xcorr_threshold),

            "--evalue_threshold",

            str(args.comet_evalue_threshold),

            "--min_deltacn",

            str(args.comet_min_deltacn),

            "--min_length",

            str(args.comet_min_length),

            "--max_length",

            str(args.comet_max_length),

        ]

        if args.keep_decoy:

            comet_post.append("--keep_decoy")

        if args.comet_strict_deltacn:

            comet_post.append("--strict_deltacn")

        print("🔄 Running Comet Class I peptide extraction...")

        try:

            result = subprocess.run(comet_post, check=True, capture_output=True, text=True)

            print(result.stdout)

        except subprocess.CalledProcessError as e:

            print("❌ Comet post-processing failed:")

            print(e.stdout)

            print(e.stderr)

            # Continue execution; do not crash the entire pipeline

    elif not comet_txt_files:

        print("ℹ️ No Comet output files found; skipping Comet post-processing")



    # 7. Merge available Class I results (merged file is generated even if only one source exists)

    fragger_cls = os.path.join(output_dir, "fragger_classI.fasta")

    comet_cls = os.path.join(output_dir, "comet_classI.fasta")

    merged_cls_out = os.path.join(output_dir, "merged_classI_dedup.fasta")

    merged_inputs = []

    if os.path.isfile(comet_cls) and os.path.getsize(comet_cls) > 0:

        merged_inputs.append(comet_cls)

    if os.path.isfile(fragger_cls) and os.path.getsize(fragger_cls) > 0:

        merged_inputs.append(fragger_cls)



    if not args.no_merge_class1 and merged_inputs:

        src_desc = []

        if comet_cls in merged_inputs:

            src_desc.append("comet_classI")

        if fragger_cls in merged_inputs:

            src_desc.append("fragger_classI")

        print(f"\n🔄 Merging Class I ({' + '.join(src_desc)}, deduplicate by peptide sequence)...")

        merge_classi_fastas_dedupe(

            merged_inputs,

            merged_cls_out,

            prefer_first_header=args.class1_merge_headers == "first",

        )

    elif not args.no_merge_class1:

        print(

            "ℹ️ merged_classI_dedup not generated (neither comet_classI nor fragger_classI has valid content)"

        )



    print(f"\n🎉 All steps complete! Results directory: {output_dir}")



if __name__ == "__main__":

    main()

