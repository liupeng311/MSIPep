import subprocess
import argparse
import os
import logging
import gzip
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, description):
    """Execute a shell command and log status"""
    logging.info(f"Start: {description}")
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Finished: {description}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {description}: {e}")
        raise

def run_hla_typing(input_fastq, output_dir, prefix):
    """
    Run OptiTypePipeline for single-end RNA-seq data to predict HLA typing.
    """
    OPTITYPE_SCRIPT = "./software/OptiType/OptiTypePipeline.py"
    REFERENCE_HLA = "./software/OptiType/data/hla_reference_rna.fasta"

    os.makedirs(output_dir, exist_ok=True)

    # Unzip fastq.gz to .fastq
    uncompressed_fastq = os.path.join(output_dir, f"{prefix}.fastq")
    with gzip.open(input_fastq, 'rb') as f_in, open(uncompressed_fastq, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Align with razers3
    fished_bam = os.path.join(output_dir, f"{prefix}_RNA.bam")
    run_command(
        f"razers3 -i 95 -m 1 -dr 0 -o {fished_bam} {REFERENCE_HLA} {uncompressed_fastq}",
        "Align to HLA reference using razers3"
    )

    # Sort and index BAM
    sorted_bam = os.path.join(output_dir, f"{prefix}_RNA_sorted.bam")
    run_command(f"samtools sort {fished_bam} -o {sorted_bam}", "Sort BAM")
    run_command(f"samtools index {sorted_bam}", "Index BAM")

    # Convert to FASTQ
    fastq_output = os.path.join(output_dir, f"{prefix}_fished.fastq")
    run_command(f"samtools bam2fq {sorted_bam} > {fastq_output}", "Convert BAM to FASTQ for OptiType")

    # Run OptiType
    run_command(
        f"python {OPTITYPE_SCRIPT} -i {fastq_output} -r -o {output_dir} -p {prefix} -v",
        "Run OptiType HLA typing"
    )

    # Clean up intermediate files
    for f in [uncompressed_fastq, fished_bam, sorted_bam, f"{sorted_bam}.bai", fastq_output]:
        if os.path.exists(f):
            os.remove(f)
            logging.info(f"Deleted intermediate file: {f}")

def main(args):
    # Use directory of input file as output path
    input_dir = os.path.dirname(args.rna_seq_1P)
    args.output_dir = input_dir

    reference_genome = "../../reference/hg38.fa"
    trimmomatic_output_base = os.path.join(args.output_dir, "trimmed_output")

    # Step 1: Trimmomatic for RNA-seq quality control
    trimmomatic_command = (
        f"java -jar ./software/Trimmomatic/trimmomatic-0.39.jar PE -phred33 -threads {args.threads} "
        f"{args.rna_seq_1P} {args.rna_seq_2P} "
        f"-baseout {trimmomatic_output_base}.fastq.gz "
        f"ILLUMINACLIP:./software/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 "
        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
    )
    run_command(trimmomatic_command, "Trimmomatic QC")

    trimmed_1P = f"{trimmomatic_output_base}_1P.fastq.gz"
    trimmed_2P = f"{trimmomatic_output_base}_2P.fastq.gz"

    # Step 2: Alignment and mutation detection pipeline
    bwa_output_sam = os.path.join(args.output_dir, "bwa.sam")
    bwa_output_bam = os.path.join(args.output_dir, "bwa.bam")
    bwa_sorted_bam = os.path.join(args.output_dir, "bwa-sort.bam")
    marked_bam = os.path.join(args.output_dir, "4mark.bam")
    add_rg_bam = os.path.join(args.output_dir, "5Add.bam")
    recal_table = os.path.join(args.output_dir, "6base.table")
    recal_bam = os.path.join(args.output_dir, "6base.bam")
    final_bam = os.path.join(args.output_dir, "7Addor-recal.bam")
    mutation_vcf = os.path.join(args.output_dir, "munation-ref.vcf")
    filtered_vcf = os.path.join(args.output_dir, "munation-ref-filter.vcf")
    final_vcf = os.path.join(args.output_dir, "munation3.vcf")
    annovar_out = os.path.join(args.output_dir, "mutadddb")

    run_command(
        f"bwa mem -t {args.threads} -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' "
        f"{reference_genome} {trimmed_1P} {trimmed_2P} > {bwa_output_sam}",
        "BWA alignment"
    )

    run_command(f"samtools fixmate -@ {args.threads} -O bam {bwa_output_sam} {bwa_output_bam}", "Convert SAM to BAM")
    os.remove(bwa_output_sam)

    run_command(f"samtools sort -@ {args.threads} -O bam -o {bwa_sorted_bam} -T {args.output_dir}/bwa-sort-temp {bwa_output_bam}", "Sort BAM")
    os.remove(bwa_output_bam)

    run_command(f"samtools index {bwa_sorted_bam}", "Index BAM")

    run_command(f"gatk MarkDuplicates -I {bwa_sorted_bam} -O {marked_bam} -M {os.path.join(args.output_dir, '4marked_dup_metrics.txt')}", "Mark duplicates")
    os.remove(bwa_sorted_bam)

    run_command(f"gatk AddOrReplaceReadGroups -I {marked_bam} -O {add_rg_bam} -ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20", "Add read groups")
    os.remove(marked_bam)

    run_command(
        f"gatk BaseRecalibrator -I {add_rg_bam} -R {reference_genome} "
        f"--known-sites ./reference/dbsnp_146.hg38.vcf "
        f"--known-sites ./reference/Mills_and_1000G_gold_standard.indels.hg38.vcf "
        f"-O {recal_table}", "Base recalibration"
    )

    run_command(f"gatk ApplyBQSR -R {reference_genome} -I {add_rg_bam} --bqsr-recal-file {recal_table} -O {recal_bam}", "Apply BQSR")
    os.remove(add_rg_bam)
    os.remove(recal_table)

    run_command(
        f"picard AddOrReplaceReadGroups I={recal_bam} O={final_bam} RGID=sample1 RGLB=library1 "
        f"RGPL=illumina RGPU=unit1 SORT_ORDER=coordinate RGSM=sample1", "Sort and add read groups"
    )
    os.remove(recal_bam)

    run_command(f"samtools index {final_bam}", "Index final BAM")

    run_command(
        f"gatk Mutect2 -R {reference_genome} -I {final_bam} -tumor tumor_sample -O {mutation_vcf} --dont-use-soft-clipped-bases",
        "Call mutations"
    )

    run_command(f"gatk FilterMutectCalls -R {reference_genome} -V {mutation_vcf} -O {filtered_vcf}", "Filter mutations")
    os.remove(mutation_vcf)

    run_command(
        f"cat {filtered_vcf} | perl -alne 'if(/^#/){{print}}else{{next unless $F[6] eq \"PASS\"; next if $F[0] =~ /_/; print}}' > {final_vcf}",
        "Manual filtering VCF"
    )
    os.remove(filtered_vcf)

    run_command(
        f"perl ./software/annovar/table_annovar.pl {final_vcf} ./humandb -buildver hg38 -out {annovar_out} "
        f"-remove -protocol refGene,gnomad_genome,ALL.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08 "
        f"-operation g,f,f,f,f,f,f,f -nastring . -vcfinput",
        "Annotate mutations with ANNOVAR"
    )

    logging.info("RNA-Seq mutation detection pipeline completed.")

    # Clean up intermediate files
    try:
        intermediate_files = [
            trimmed_1P,
            trimmed_2P,
            final_bam
        ]
        intermediate_files.extend([
            os.path.join(args.output_dir, f) for f in os.listdir(args.output_dir) if f.endswith(".bai")
        ])
        for f in intermediate_files:
            if os.path.exists(f):
                os.remove(f)
        logging.info("Intermediate files cleaned up.")
    except Exception as e:
        logging.warning(f"Error cleaning intermediate files: {e}")

    # HLA typing (based on RNA-Seq 1P read)
    run_hla_typing(args.rna_seq_1P, args.output_dir, prefix="HLA_typing")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA-Seq mutation detection and HLA typing pipeline")
    parser.add_argument('-rna_seq_1P', required=True, help="Path to first RNA-Seq read file (FASTQ.gz)")
    parser.add_argument('-rna_seq_2P', required=True, help="Path to second RNA-Seq read file (FASTQ.gz)")
    parser.add_argument('-threads', type=int, default=16, help="Number of threads to use")

    args = parser.parse_args()
    main(args)

#python rna_mutation_pipeline.py   -rna_seq_1P /data/rnaseq/sample_1.fastq.gz   -rna_seq_2P /data/rnaseq/sample_2.fastq.gz \ -threads 16
