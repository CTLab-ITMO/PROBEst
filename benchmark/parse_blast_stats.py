from Bio import SeqIO
import pandas as pd
import sys

# Путь к файлам
fasta_path = sys.argv[1]
blast_path = sys.argv[2]

# Загрузка длин зондов из FASTA
probe_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
probe_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

# Загрузка BLAST-выхода
cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qseq", "sseq"
]
df = pd.read_csv(blast_path, sep="\t", names=cols)

# Добавляем длину зонда
df["probe_len"] = df["qseqid"].map(probe_lengths)
df["full_qseq"] = df["qseqid"].map(probe_seqs)

# Фильтрация по критерию C75
df["coverage"] = df["length"] / df["probe_len"]
df["C75"] = (df["pident"] >= 75) & (df["coverage"] >= 0.75)

# Отфильтрованные хиты
filtered = df[df["C75"]].copy()
filtered["is_offtarget"] = filtered["sseqid"].str.startswith("off_")

# Подсчёт
target_hits = filtered[~filtered["is_offtarget"]].groupby("qseqid")["sseqid"].nunique().rename("target_hits")
offtarget_hits = filtered[filtered["is_offtarget"]].groupby("qseqid")["sseqid"].nunique().rename("offtarget_hits")

# Мультимаппинги (по >1 хиту на одну sseqid)
multi = (
    filtered.groupby(["qseqid", "sseqid"]).size()
    .reset_index(name="n")
    .query("n > 1")
    .groupby("qseqid").size()
    .rename("multimapping_events")
)

# Последовательности зондов
seqs = pd.Series(probe_seqs, name="probe_seq")

# Объединяем
summary = pd.concat([target_hits, offtarget_hits, multi, seqs], axis=1).fillna(0)
summary[["target_hits", "offtarget_hits", "multimapping_events"]] = summary[["target_hits", "offtarget_hits", "multimapping_events"]].astype(int)

# Сохраняем
summary.to_csv("blast_probe_metrics_c75.csv")
