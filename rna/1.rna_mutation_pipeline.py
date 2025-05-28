
import subprocess
import argparse
import os
import logging

# 配置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, description):
    """执行Shell命令并记录日志和错误处理"""
    logging.info(f"开始：{description}")
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f"完成：{description}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{description} 出错: {e}")
        raise

def main(args):
    # 使用输入文件所在目录作为输出目录
    input_dir = os.path.dirname(args.rna_seq_1P)
    args.output_dir = input_dir

    reference_genome = "../../reference/hg38.fa"

    # 定义Trimmomatic输出基本名（不再使用子目录）
    trimmomatic_output_base = os.path.join(args.output_dir, "trimmed_output")

    # Trimmomatic QC步骤
    trimmomatic_command = (
        f"java -jar ../../software/Trimmomatic/trimmomatic-0.39.jar PE -phred33 -threads {args.threads} "
        f"{args.rna_seq_1P} {args.rna_seq_2P} "
        f"-baseout {trimmomatic_output_base}.fastq.gz "
        f"ILLUMINACLIP:../../software/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 "
        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
    )
    run_command(trimmomatic_command, "Trimmomatic质量控制")

    # 使用Trimmomatic生成的双端文件
    trimmed_1P = f"{trimmomatic_output_base}_1P.fastq.gz"
    trimmed_2P = f"{trimmomatic_output_base}_2P.fastq.gz"

    # 定义中间及最终文件路径（全部在output_dir中）
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

    # BWA比对
    bwa_command = (
        f"bwa mem -t {args.threads} -R '@RG\tID:foo\tSM:bar\tLB:library1' "
        f"{reference_genome} {trimmed_1P} {trimmed_2P} > {bwa_output_sam}"
    )
    run_command(bwa_command, "BWA比对")

    # SAM转BAM
    run_command(f"samtools fixmate -@ {args.threads} -O bam {bwa_output_sam} {bwa_output_bam}", "SAM转BAM")
    os.remove(bwa_output_sam)
    # BAM排序
    run_command(f"samtools sort -@ {args.threads} -O bam -o {bwa_sorted_bam} -T {args.output_dir}/bwa-sort-temp {bwa_output_bam}", "BAM排序")
    os.remove(bwa_output_bam)
    # 生成BAM索引
    run_command(f"samtools index {bwa_sorted_bam}", "生成BAM索引")
    # 标记重复
    run_command(f"gatk MarkDuplicates -I {bwa_sorted_bam} -O {marked_bam} -M {os.path.join(args.output_dir, '4marked_dup_metrics.txt')}", "标记重复")
    os.remove(bwa_sorted_bam)
    # 添加读组信息
    run_command(f"gatk AddOrReplaceReadGroups -I {marked_bam} -O {add_rg_bam} -ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20", "添加读组")
    os.remove(marked_bam)
    # 碱基重校准
    run_command(f"gatk BaseRecalibrator -I {add_rg_bam} -R {reference_genome} --known-sites ../../reference/dbsnp_146.hg38.vcf --known-sites ../../reference/Mills_and_1000G_gold_standard.indels.hg38.vcf -O {recal_table}", "碱基重校准")
    # 应用BQSR
    run_command(f"gatk ApplyBQSR -R {reference_genome} -I {add_rg_bam} --bqsr-recal-file {recal_table} -O {recal_bam}", "应用BQSR")
    os.remove(add_rg_bam)
    os.remove(recal_table)
    # 再次添加读组并排序（Picard）
    run_command(f"picard AddOrReplaceReadGroups I={recal_bam} O={final_bam} RGID=sample1 RGLB=library1 RGPL=illumina RGPU=unit1 SORT_ORDER=coordinate RGSM=sample1", "排序并添加读组")
    os.remove(recal_bam)
    # 索引BAM文件
    run_command(f"samtools index {final_bam}", "索引BAM文件")
    # 突变检测
    run_command(f"gatk Mutect2 -R {reference_genome} -I {final_bam} -tumor tumor_sample -O {mutation_vcf} --dont-use-soft-clipped-bases", "突变检测")
    # 突变过滤
    run_command(f"gatk FilterMutectCalls -R {reference_genome} -V {mutation_vcf} -O {filtered_vcf}", "突变过滤")
    os.remove(mutation_vcf)

    # 手动过滤VCF
    try:
        subprocess.run(
            f"cat {filtered_vcf} | perl -alne 'if(/^#/){{print}}else{{next unless $F[6] eq \"PASS\"; next if $F[0] =~ /_/; print}}' > {final_vcf}",
            shell=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print("Error occurred during VCF filtering:", e)
        raise

    os.remove(filtered_vcf)

    # 使用ANNOVAR注释
    run_command(
        f"perl ../../software/annovar/table_annovar.pl {final_vcf} ../../software/annovar/humandb -buildver hg38 -out {annovar_out} "
        f"-remove -protocol refGene,gnomad_genome,ALL.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08 "
        f"-operation g,f,f,f,f,f,f,f -nastring . -vcfinput",
        "突变注释"
    )

    logging.info("RNA-Seq分析管道已成功完成！")

    try:
        # 仅清理指定的中间文件
        intermediate_files = [
            os.path.join(args.output_dir, "trimmed_output_1P.fastq.gz"),
            os.path.join(args.output_dir, "trimmed_output_2P.fastq.gz"),
            os.path.join(args.output_dir, "7Addor-recal.bam")
        ]
        # 删除所有 .bai 文件
        intermediate_files.extend([
            os.path.join(args.output_dir, f) for f in os.listdir(args.output_dir) if f.endswith(".bai")
        ])

        for f in intermediate_files:
            if os.path.exists(f):
                os.remove(f)
        logging.info("中间文件清理完成。")
    except Exception as e:
        logging.warning(f"中间文件清理时发生错误: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA-Seq分析管道")
 #   parser.add_argument('--reference_genome', required=True, help="参考基因组的路径")
    parser.add_argument('-rna_seq_1P', required=True, help="第一条RNA-Seq读段文件路径")
    parser.add_argument('-rna_seq_2P', required=True, help="第二条RNA-Seq读段文件路径")
    parser.add_argument('-threads', type=int, default=16, help="指定线程数")

    args = parser.parse_args()
    main(args)

#python rna_mutation_pipeline_modified.py -rna_seq_1P /data/rnaseq/sample_1.fastq  -rna_seq_2P /data/rnaseq/sample_2.fastq
