#!/usr/bin/env python3
"""
RNA-seq Somatic Mutation Calling Pipeline——MAF (Single-End)
============================================================
Goal: Use tumor single-end RNA-seq data only, apply strict filtering, and output a somatic mutation MAF file suitable for neoantigen prediction.
"""

import subprocess
import argparse
import os
import logging
import sys
import glob

# ── Logging setup ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


# ── Utility functions ──────────────────────────────────────────────────────────────────
def run_command(command: str, description: str):
    """Run a shell command; raise on failure and log errors."""
    logger.info(f"▶ Start: {description}")
    logger.debug(f"  CMD: {command}")
    try:
        subprocess.run(command, shell=True, check=True, executable="/bin/bash")
        logger.info(f"✔ Done: {description}")
    except subprocess.CalledProcessError as e:
        logger.error(f"✘ Failed: {description}\n  Error: {e}")
        raise


def file_exists_and_nonempty(path: str) -> bool:
    """Check if file exists and is non-empty (for resume from checkpoint)."""
    return os.path.isfile(path) and os.path.getsize(path) > 0


def validate_fastq_input(path: str, arg_name: str) -> str:
    """
    Accept .fastq/.fq (uncompressed) or .fastq.gz/.fq.gz/.gz (also .fsatq.gz typo).
    Trimmomatic reads these formats directly; downstream uses gzip output.
    """
    if not file_exists_and_nonempty(path):
        raise FileNotFoundError(f"{arg_name} file missing or empty: {path}")

    lowered = path.lower()
    allowed_suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fsatq.gz", ".gz")
    if not lowered.endswith(allowed_suffixes):
        raise ValueError(
            f"{arg_name} supports only .fastq/.fq/.fastq.gz/.fq.gz/.gz (also .fsatq.gz): {path}"
        )
    return path


def resolve_trimmomatic_output(path: str) -> str:
    """Check that Trimmomatic single-end output exists and is non-empty."""
    if file_exists_and_nonempty(path):
        return path
    raise FileNotFoundError(f"Trimmomatic single-end output not found: {path}")


def resolve_existing_path(candidates, description: str) -> str:
    """
    Pick the first existing path from candidates (works from different working directories).
    """
    for p in candidates:
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(
        f"{description} not found. Checked paths: {candidates}"
    )


def skip_if_exists(path: str, description: str) -> bool:
    """Skip step if target file exists; return True if skipped."""
    if file_exists_and_nonempty(path):
        logger.info(f"⏭ Skip (exists): {description} → {path}")
        return True
    return False


def detect_latest_checkpoint(checkpoints: list) -> int:
    """
    Scan stage marker outputs front to back; return furthest stage reached (0 = no intermediate results).
    Intermediate files may be deleted after later steps, so do not rely on previous step filename alone.
    """
    latest = 0
    for stage_id, label, path in checkpoints:
        if file_exists_and_nonempty(path):
            latest = max(latest, stage_id)
            logger.info(f"  ✓ Stage {stage_id} output exists: {label} → {path}")
    return latest


def should_skip_step(
    latest_stage: int,
    min_stage_to_skip: int,
    path: str,
    description: str,
) -> bool:
    """
    Resume: skip if disk shows a later stage (latest_stage >= min_stage_to_skip);
    otherwise skip if path exists (same as skip_if_exists).
    min_stage_to_skip matches this step's stage id (reached after completing this step).
    """
    if latest_stage >= min_stage_to_skip:
        logger.info(
            f"⏭ Skip (resume, stage≥{min_stage_to_skip} detected): {description}"
        )
        return True
    return skip_if_exists(path, description)


def all_files_exist(paths) -> bool:
    """Check that all given files exist and are non-empty."""
    return all(file_exists_and_nonempty(p) for p in paths)


def ensure_required_inputs(step_name: str, required_inputs):
    """
    Verify inputs for a step exist; raise if missing and suggest rerunning earlier steps.
    """
    missing = [p for p in required_inputs if not file_exists_and_nonempty(p)]
    if missing:
        raise FileNotFoundError(
            f"{step_name} required input files missing; cannot rerun this step automatically. Missing: {missing}"
        )


def remove_file_if_exists(path: str):
    """Safely delete a single file (skip if absent)."""
    if os.path.isfile(path):
        os.remove(path)
        logger.info(f"🗑 Deleted intermediate file: {path}")


def cleanup_files(paths, description: str):
    """Batch safe delete of intermediate files (supports glob)."""
    logger.info(f"🧹 Cleaning intermediate files: {description}")
    for p in paths:
        if any(ch in p for ch in ["*", "?", "["]):
            for matched in glob.glob(p):
                remove_file_if_exists(matched)
        else:
            remove_file_if_exists(p)


# ── Main workflow ────────────────────────────────────────────────────────────────────
def main(args):

    # ── Path configuration ──────────────────────────────────────────────────────────────
    ref               = args.ref
    star_index        = args.star_index
    redi_bed          = args.redi_bed          # REDIportal (hg38, chr prefix)
    gnomad            = args.gnomad            # af-only-gnomad.hg38.vcf.gz (indexed)
    dbsnp             = args.dbsnp             # dbsnp_146.hg38.vcf.gz (indexed)
    common_biallelic  = args.common_biallelic  # small_exac_common_3.hg38.vcf.gz
    pon               = args.pon               # panel_of_normals.vcf.gz (optional, strongly recommended)

    out_dir = os.path.abspath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Append file log for this run
    fh = logging.FileHandler(os.path.join(out_dir, "pipeline.log"))
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    sample  = args.sample_name
    threads = args.threads
    fq1     = validate_fastq_input(args.fq1, "--fq1")

    # ── File naming ──────────────────────────────────────────────────────────────
    star_prefix         = os.path.join(out_dir, f"{sample}_")
    star_bam            = f"{star_prefix}Aligned.sortedByCoord.out.bam"
    marked_bam          = f"{out_dir}/{sample}_md.bam"
    split_bam           = f"{out_dir}/{sample}_split.bam"
    recal_table         = f"{out_dir}/{sample}_recal.table"
    final_bam           = f"{out_dir}/{sample}_final.bam"
    f1r2_tar            = f"{out_dir}/{sample}_f1r2.tar.gz"
    read_orientation    = f"{out_dir}/{sample}_read_orientation_model.tar.gz"
    pileup_table        = f"{out_dir}/{sample}_pileups.table"
    contamination_table = f"{out_dir}/{sample}_contamination.table"
    raw_vcf             = f"{out_dir}/{sample}_raw.vcf.gz"
    raw_vcf_stats       = f"{out_dir}/{sample}_raw.vcf.gz.stats"
    filt_vcf            = f"{out_dir}/{sample}_filtered.vcf.gz"
    pass_vcf            = f"{out_dir}/{sample}_pass.vcf.gz"
    pass_vcf_plain      = pass_vcf.replace(".vcf.gz", ".vcf")
    clean_vcf           = f"{out_dir}/{sample}_no_editing.vcf"   # produced directly after REDIportal filter
    out_maf             = f"{out_dir}/{sample}.maf"
    trimmomatic_output  = os.path.join(out_dir, f"{sample}_trimmomatic.fastq.gz")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    trimmomatic_jar = resolve_existing_path(
        [
            args.trimmomatic_jar,
            os.path.join(script_dir, args.trimmomatic_jar),
        ],
        "Trimmomatic jar",
    )
    trimmomatic_adapters = resolve_existing_path(
        [
            args.trimmomatic_adapters,
            os.path.join(script_dir, args.trimmomatic_adapters),
        ],
        "Trimmomatic adapters",
    )

    # ════════════════════════════════════════════════════════════════════════
    # Step 0. Trimmomatic preprocessing (single-end)
    # ════════════════════════════════════════════════════════════════════════
    if not file_exists_and_nonempty(trimmomatic_output):
        trimmomatic_command = (
            f"java -jar {trimmomatic_jar} SE -phred33 -threads {threads} "
            f"{fq1} {trimmomatic_output} "
            f"ILLUMINACLIP:{trimmomatic_adapters}:2:30:10 "
            f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
        )
        run_command(trimmomatic_command, "Trimmomatic SE")
    else:
        logger.info("⏭ Skip (exists): Trimmomatic SE")

    fq1 = resolve_trimmomatic_output(trimmomatic_output)
    logger.info(f"Continuing with Trimmomatic output: {fq1}")

    # ── Resume: scan stage marker outputs (intermediates may be deleted later) ─────────
    checkpoint_specs = [
        (1, "STAR alignment BAM", star_bam),
        (2, "MarkDuplicates BAM", marked_bam),
        (3, "SplitNCigar BAM", split_bam),
        (4, "BaseRecalibrator table", recal_table),
        (5, "ApplyBQSR final BAM", final_bam),
        (6, "Mutect2 raw VCF", raw_vcf),
        (7, "F1R2 orientation model", read_orientation),
        (8, "GetPileupSummaries table", pileup_table),
        (9, "CalculateContamination table", contamination_table),
        (10, "FilterMutectCalls VCF", filt_vcf),
        (11, "PASS VCF", pass_vcf),
        (12, "Decompressed PASS VCF", pass_vcf_plain),
        (13, "REDI-filtered VCF", clean_vcf),
    ]
    logger.info("Scanning output directory for existing intermediate results (resume)…")
    latest = detect_latest_checkpoint(checkpoint_specs)
    logger.info(
        f"Resume check: furthest stage {latest}/13"
        + (" (will continue from incomplete steps)" if latest < 13 else " (all checkpoint products present)")
    )

    # ════════════════════════════════════════════════════════════════════════
    # Step 1. STAR 2-pass alignment (single-end, embedded RG info, required by Mutect2)
    # ════════════════════════════════════════════════════════════════════════
    if not should_skip_step(latest, 1, star_bam, "STAR Alignment"):
        star_cmd = (
            f"STAR --runThreadN {threads} "
            f"--genomeDir {star_index} "
            f"--readFilesIn {fq1} "
            f"--readFilesCommand zcat "
            f"--outBAMsortingThreadN 4 "
            f"--limitBAMsortRAM 10000000000 "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outFileNamePrefix {star_prefix} "
            f"--twopassMode Basic "
            f"--outSAMattributes NH HI AS NM MD "
            f"--outSAMattrRGline ID:{sample} SM:{sample} PL:ILLUMINA LB:lib1 PU:{sample}"
        )
        run_command(star_cmd, "STAR 2-pass Alignment")

    # ════════════════════════════════════════════════════════════════════════
    # Step 2. GATK MarkDuplicates (dedup + index)
    # ════════════════════════════════════════════════════════════════════════
    if not should_skip_step(latest, 2, marked_bam, "MarkDuplicates"):
        run_command(
            f"gatk MarkDuplicates "
            f"-I {star_bam} "
            f"-O {marked_bam} "
            f"-M {out_dir}/{sample}_dup_metrics.txt "
            f"--CREATE_INDEX true "
            f"--VALIDATION_STRINGENCY SILENT",
            "Mark Duplicates"
        )

    if not should_skip_step(latest, 3, split_bam, "SplitNCigarReads"):
        run_command(
            f"gatk SplitNCigarReads "
            f"-R {ref} "
            f"-I {marked_bam} "
            f"-O {split_bam} "
            f"--skip-mq-transform false",
            "SplitNCigarReads"
        )

    if not should_skip_step(latest, 4, recal_table, "BaseRecalibrator"):
        run_command(
            f"gatk BaseRecalibrator "
            f"-R {ref} "
            f"-I {split_bam} "
            f"--known-sites {dbsnp} "
            f"--known-sites {gnomad} "
            f"-O {recal_table}",
            "Base Recalibration"
        )

    if not should_skip_step(latest, 5, final_bam, "ApplyBQSR"):
        run_command(
            f"gatk ApplyBQSR "
            f"-R {ref} "
            f"-I {split_bam} "
            f"--bqsr-recal-file {recal_table} "
            f"-O {final_bam} "
            f"--create-output-bam-index true",
            "Apply BQSR"
        )

    mutect2_expected_outputs = [raw_vcf, raw_vcf_stats, f1r2_tar]
    skip_mutect2 = should_skip_step(latest, 6, raw_vcf, "Mutect2")
    if skip_mutect2 and not all_files_exist(mutect2_expected_outputs):
        logger.warning(
            "⚠ Mutect2 key outputs incomplete (missing raw_vcf/raw_vcf_stats/f1r2); will rerun Step 5 Mutect2."
        )
        ensure_required_inputs("Mutect2", [final_bam])
        skip_mutect2 = False

    if not skip_mutect2:
        pon_arg = f"--panel-of-normals {pon}" if pon and os.path.isfile(pon) else ""
        run_command(
            f"gatk Mutect2 "
            f"-R {ref} "
            f"-I {final_bam} "
            f"--tumor-sample {sample} "
            f"--germline-resource {gnomad} "
            f"--af-of-alleles-not-in-resource 0.0000025 "
            f"--dont-use-soft-clipped-bases true "
            f"--f1r2-tar-gz {f1r2_tar} "
            f"{pon_arg} "
            f"-O {raw_vcf}",
            "Mutect2 Variant Calling (tumor-only)"
        )

    if not should_skip_step(latest, 7, read_orientation, "LearnReadOrientationModel"):
        ensure_required_inputs("LearnReadOrientationModel", [f1r2_tar])
        run_command(
            f"gatk LearnReadOrientationModel "
            f"-I {f1r2_tar} "
            f"-O {read_orientation}",
            "Learn Read Orientation Model (F1R2)"
        )

    if not should_skip_step(latest, 8, pileup_table, "GetPileupSummaries"):
        run_command(
            f"gatk GetPileupSummaries "
            f"-I {final_bam} "
            f"-V {common_biallelic} "
            f"-L {common_biallelic} "
            f"-O {pileup_table}",
            "Get Pileup Summaries"
        )

    if not should_skip_step(latest, 9, contamination_table, "CalculateContamination"):
        run_command(
            f"gatk CalculateContamination "
            f"-I {pileup_table} "
            f"-O {contamination_table}",
            "Calculate Contamination"
        )

    if not should_skip_step(latest, 10, filt_vcf, "FilterMutectCalls"):
        run_command(
            f"gatk FilterMutectCalls "
            f"-R {ref} "
            f"-V {raw_vcf} "
            f"--stats {raw_vcf_stats} "
            f"--contamination-table {contamination_table} "
            f"--orientation-bias-artifact-priors {read_orientation} "
            f"-O {filt_vcf}",
            "Filter Mutect Calls"
        )

    if not should_skip_step(latest, 11, pass_vcf, "SelectVariants PASS"):
        run_command(
            f"gatk SelectVariants "
            f"-R {ref} "
            f"-V {filt_vcf} "
            f"--exclude-filtered true "
            f"-O {pass_vcf}",
            "Select PASS Variants Only"
        )

    if not should_skip_step(latest, 12, pass_vcf_plain, "Decompress PASS VCF"):
        run_command(
            f"bcftools view {pass_vcf} > {pass_vcf_plain}",
            "Decompress PASS VCF for bedtools"
        )

    if not should_skip_step(latest, 13, clean_vcf, "Exclude REDIportal"):
        run_command(
            f"bedtools intersect -v -a {pass_vcf_plain} -b {redi_bed} -header > {clean_vcf}",
            "Exclude Known RNA-Editing Sites (REDIportal)"
        )

    path_file = os.path.join(out_dir, "clean_vcf.path")
    with open(path_file, "w") as f:
        f.write(clean_vcf)
    logger.info(f"📄 clean_vcf path written: {path_file}")

    # Delete intermediates only after final output exists
    if file_exists_and_nonempty(clean_vcf):
        cleanup_files(
            [
                star_bam,
                f"{star_bam}.bai",
                marked_bam,
                f"{marked_bam}.bai",
                split_bam,
                f"{split_bam}.bai",
                recal_table,
                final_bam,
                f"{final_bam}.bai",
                f1r2_tar,
                read_orientation,
                pileup_table,
                contamination_table,
                raw_vcf,
                f"{raw_vcf}.tbi",
                raw_vcf_stats,
                filt_vcf,
                f"{filt_vcf}.tbi",
                pass_vcf,
                f"{pass_vcf}.tbi",
                pass_vcf_plain,
            ],
            "Final intermediate file cleanup"
        )

    logger.info("=" * 60)
    logger.info("✅ Python pipeline complete (Steps 1–10)!")
    logger.info(f"   clean_vcf : {clean_vcf}")
    logger.info("   Next: run_pipeline.sh switches to vep env for vcf2maf annotation")
    logger.info("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="RNA-seq Somatic Mutation_ MAF Pipeline v2.3 (single-end) | Step 1-10 | invoked by run_pipeline.sh",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--fq1",
        required=True,
        help="FASTQ (single-end; .fastq/.fq/.fastq.gz/.fq.gz/.gz, also .fsatq.gz)"
    )
    parser.add_argument("--sample_name", required=True, help="Sample ID (must match BAM header SM)")
    parser.add_argument("--output_dir", default="./output", help="Output directory")
    parser.add_argument("--threads", type=int, default=16, help="Number of parallel threads")
    parser.add_argument(
        "--trimmomatic_jar",
        default="/data/liup/software/Trimmomatic/trimmomatic-0.39.jar",
        help="Trimmomatic jar path (relative to cwd or script directory)",
    )
    parser.add_argument(
        "--trimmomatic_adapters",
        default="/data/liup/software/Trimmomatic/adapters/TruSeq3-SE.fa",
        help="Trimmomatic adapter fasta path (relative to cwd or script directory)",
    )

    parser.add_argument("--ref",
        default="/data/liup/data/reference/hg38/hg38.fa",
        help="Reference genome FASTA (requires .fai and .dict indexes)")
    parser.add_argument("--star_index",
        default="/data/liup/data/reference/STAR_index",
        help="STAR genome index directory")

    parser.add_argument("--dbsnp",
        default="/data/liup/data/reference/dbsnp_146.hg38.vcf.gz",
        help="dbSNP VCF (BQSR known sites, requires .tbi index)")
    parser.add_argument("--gnomad",
        default="/data/liup/data/reference/af-only-gnomad.hg38.vcf.gz",
        help="gnomAD AF-only VCF (Mutect2 germline resource, requires .tbi index)")
    parser.add_argument("--common_biallelic",
        default="/data/liup/data/reference/small_exac_common_3.hg38.vcf.gz",
        help="Common biallelic SNP VCF (GetPileupSummaries, requires .tbi index)")
    parser.add_argument("--pon",
        default="/data/liup/data/reference/1000g_pon.hg38.vcf.gz",
        help="Panel of Normals VCF (strongly recommended; greatly reduces tumor-only false positives)")

    parser.add_argument("--redi_bed",
        default="/data/liup/data/reference/REDIportal_hg38.bed",
        help="REDIportal RNA editing sites BED (hg38, chr prefix)")

    args = parser.parse_args()
    main(args)
