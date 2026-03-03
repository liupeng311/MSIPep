#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import argparse
import os
import logging
import shutil

# 日志配置
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, description):
    """执行命令并打印日志"""
    logging.info(f"Start: {description}")
    logging.info(f"Command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Finished: {description}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {description}: {e}")
        raise

def process_sample(sample_label, rna_seq_1P, rna_seq_2P, output_dir, threads, reference_genome):
    """
    对单个样本（肿瘤或正常）从 FASTQ 到“处理完的 BAM”的流程：
    Trimmomatic → BWA → samtools → GATK MarkDuplicates / BQSR → picard AddOrReplaceReadGroups
    返回：最终 BAM 路径
    """
    sample_dir = os.path.join(output_dir, sample_label)
    os.makedirs(sample_dir, exist_ok=True)

    # 1. Trimmomatic 质控
    trimmomatic_output_base = os.path.join(sample_dir, f"{sample_label}_trimmed")
    trimmed_1P = f"{trimmomatic_output_base}_1P.fastq.gz"
    trimmed_2P = f"{trimmomatic_output_base}_2P.fastq.gz"

    trimmomatic_command = (
        f"java -jar ./software/Trimmomatic/trimmomatic-0.39.jar PE -phred33 -threads {threads} "
        f"{rna_seq_1P} {rna_seq_2P} "
        f"-baseout {trimmomatic_output_base}.fastq.gz "
        f"ILLUMINACLIP:./software/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 "
        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
    )
    run_command(trimmomatic_command, f"Trimmomatic QC ({sample_label})")

    # 2. 比对及 BAM 处理
    bwa_output_sam = os.path.join(sample_dir, f"{sample_label}_bwa.sam")
    bwa_output_bam = os.path.join(sample_dir, f"{sample_label}_bwa.bam")
    bwa_sorted_bam = os.path.join(sample_dir, f"{sample_label}_bwa-sort.bam")
    marked_bam = os.path.join(sample_dir, f"{sample_label}_4mark.bam")
    add_rg_bam = os.path.join(sample_dir, f"{sample_label}_5Add.bam")
    recal_table = os.path.join(sample_dir, f"{sample_label}_6base.table")
    recal_bam = os.path.join(sample_dir, f"{sample_label}_6base.bam")
    final_bam = os.path.join(sample_dir, f"{sample_label}_7Addor-recal.bam")

    # BWA 比对
    run_command(
        f"bwa mem -t {threads} -R '@RG\\tID:{sample_label}\\tSM:{sample_label}\\tLB:library1' "
        f"{reference_genome} {trimmed_1P} {trimmed_2P} > {bwa_output_sam}",
        f"BWA alignment ({sample_label})"
    )

    # SAM → BAM
    run_command(
        f"samtools fixmate -@ {threads} -O bam {bwa_output_sam} {bwa_output_bam}",
        f"Convert SAM to BAM ({sample_label})"
    )
    os.remove(bwa_output_sam)

    # 排序
    run_command(
        f"samtools sort -@ {threads} -O bam -o {bwa_sorted_bam} "
        f"-T {sample_dir}/{sample_label}_bwa-sort-temp {bwa_output_bam}",
        f"Sort BAM ({sample_label})"
    )
    os.remove(bwa_output_bam)

    # 建索引
    run_command(f"samtools index {bwa_sorted_bam}", f"Index BAM ({sample_label})")

    # 标记重复
    run_command(
        f"gatk MarkDuplicates -I {bwa_sorted_bam} -O {marked_bam} "
        f"-M {os.path.join(sample_dir, f'{sample_label}_marked_dup_metrics.txt')}",
        f"Mark duplicates ({sample_label})"
    )
    os.remove(bwa_sorted_bam)

    # 添加 read group（GATK）
    run_command(
        f"gatk AddOrReplaceReadGroups -I {marked_bam} -O {add_rg_bam} "
        f"-ID {sample_label} -LB lib1 -PL illumina -PU unit1 -SM {sample_label}",
        f"Add read groups ({sample_label})"
    )
    os.remove(marked_bam)

    # BaseRecalibrator
    run_command(
        f"gatk BaseRecalibrator -I {add_rg_bam} -R {reference_genome} "
        f"--known-sites ./reference/dbsnp_146.hg38.vcf "
        f"--known-sites ./reference/Mills_and_1000G_gold_standard.indels.hg38.vcf "
        f"-O {recal_table}",
        f"Base recalibration ({sample_label})"
    )

    # ApplyBQSR
    run_command(
        f"gatk ApplyBQSR -R {reference_genome} -I {add_rg_bam} "
        f"--bqsr-recal-file {recal_table} -O {recal_bam}",
        f"Apply BQSR ({sample_label})"
    )
    os.remove(add_rg_bam)
    os.remove(recal_table)

    # 再用 picard 统一 read group & 排序（保持与你原脚本风格一致）
    run_command(
        f"picard AddOrReplaceReadGroups I={recal_bam} O={final_bam} "
        f"RGID={sample_label} RGLB=library1 RGPL=illumina RGPU=unit1 RGSM={sample_label} "
        f"SORT_ORDER=coordinate",
        f"Sort and add read groups (picard, {sample_label})"
    )
    os.remove(recal_bam)

    # 最终 BAM 建索引
    run_command(f"samtools index {final_bam}", f"Index final BAM ({sample_label})")

    # 清理部分中间文件（保留最终 BAM 及其索引）
    try:
        intermediate_files = [
            trimmed_1P,
            trimmed_2P,
        ]
        bai_files = [
            os.path.join(sample_dir, f) for f in os.listdir(sample_dir) if f.endswith(".bai")
            and not f.startswith(os.path.basename(final_bam))  # 保留 final_bam.bai
        ]
        intermediate_files.extend(bai_files)
        for f in intermediate_files:
            if os.path.exists(f):
                os.remove(f)
        logging.info(f"Intermediate files cleaned up for {sample_label}.")
    except Exception as e:
        logging.warning(f"Error cleaning intermediate files for {sample_label}: {e}")

    return final_bam


def main(args):
    # 参考基因组路径（根据你自己的目录调整）
    reference_genome = "../../reference/hg38.fa"
    os.makedirs(args.output_dir, exist_ok=True)

    # 处理肿瘤和正常样本
    tumor_sample_label = "tumor_sample"
    normal_sample_label = "normal_sample"

    tumor_bam = process_sample(
        tumor_sample_label,
        args.tumor_rna_1P,
        args.tumor_rna_2P,
        args.output_dir,
        args.threads,
        reference_genome,
    )

    normal_bam = process_sample(
        normal_sample_label,
        args.normal_rna_1P,
        args.normal_rna_2P,
        args.output_dir,
        args.threads,
        reference_genome,
    )

    # 3. Tumor-Normal 突变检测（Mutect2）
    mutation_vcf = os.path.join(args.output_dir, "somatic_mutations_raw.vcf")
    filtered_vcf = os.path.join(args.output_dir, "somatic_mutations_filtered.vcf")
    final_vcf = os.path.join(args.output_dir, "somatic_mutations_pass.vcf")
    annovar_out = os.path.join(args.output_dir, "somatic_mutations_annovar")

    # Mutect2（成对肿瘤-正常）
    run_command(
        f"gatk Mutect2 -R {reference_genome} "
        f"-I {tumor_bam} -I {normal_bam} "
        f"-tumor {tumor_sample_label} -normal {normal_sample_label} "
        f"-O {mutation_vcf} --dont-use-soft-clipped-bases",
        "Call somatic mutations (tumor-normal) with Mutect2"
    )

    # FilterMutectCalls
    run_command(
        f"gatk FilterMutectCalls -R {reference_genome} -V {mutation_vcf} -O {filtered_vcf}",
        "Filter somatic mutations"
    )
    os.remove(mutation_vcf)

    # 额外过滤：只保留 PASS 且去掉带下划线染色体（与你原来脚本一致）
    run_command(
        f"cat {filtered_vcf} | perl -alne 'if(/^#/){{print}}else{{next unless $F[6] eq \"PASS\"; "
        f"next if $F[0] =~ /_/; print}}' > {final_vcf}",
        "Manual filtering VCF"
    )
    os.remove(filtered_vcf)

    # 4. ANNOVAR 注释
    run_command(
        f"perl ./software/annovar/table_annovar.pl {final_vcf} ./humandb -buildver hg38 -out {annovar_out} "
        f"-remove -protocol refGene,gnomad_genome,ALL.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,"
        f"AMR.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08 "
        f"-operation g,f,f,f,f,f,f,f -nastring . -vcfinput",
        "Annotate somatic mutations with ANNOVAR"
    )

    logging.info("Tumor-normal RNA-Seq somatic mutation detection pipeline completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tumor-normal RNA-Seq somatic mutation detection pipeline (Trimmomatic + BWA + GATK + ANNOVAR)"
    )
    parser.add_argument(
        '-tumor_rna_1P',
        required=True,
        help="Tumor RNA-Seq read1 (FASTQ.gz)"
    )
    parser.add_argument(
        '-tumor_rna_2P',
        required=True,
        help="Tumor RNA-Seq read2 (FASTQ.gz)"
    )
    parser.add_argument(
        '-normal_rna_1P',
        required=True,
        help="Normal RNA-Seq read1 (FASTQ.gz)"
    )
    parser.add_argument(
        '-normal_rna_2P',
        required=True,
        help="Normal RNA-Seq read2 (FASTQ.gz)"
    )
    parser.add_argument(
        '-output_dir',
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        '-threads',
        type=int,
        default=16,
        help="Number of threads to use"
    )

    args = parser.parse_args()
    main(args)

    # 示例：
    # python rna_tumor_normal_mutation_pipeline.py \
    #   -tumor_rna_1P /data/tumor/sample_tumor_1.fastq.gz \
    #   -tumor_rna_2P /data/tumor/sample_tumor_2.fastq.gz \
    #   -normal_rna_1P /data/normal/sample_normal_1.fastq.gz \
    #   -normal_rna_2P /data/normal/sample_normal_2.fastq.gz \
    #   -output_dir /data/output/tumor_normal_pipeline \
    #   -threads 16