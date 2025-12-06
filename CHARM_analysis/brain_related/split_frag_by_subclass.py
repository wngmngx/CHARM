import pandas as pd

CHUNK_SIZE = 10_000_000
FRAGMENT_COLUMNS = ["chrom", "start", "end", "CellID", "count"]


def load_cell_to_subclass(metadata_path: str, selected_subclasses):
    metadata = pd.read_csv(metadata_path, sep="\t", usecols=["CellID", "SubClass"])
    return metadata[metadata["SubClass"].isin(selected_subclasses)]


def split_fragments_by_subclass(
    fragment_path: str, metadata_path: str, selected_subclasses
) -> None:
    metadata = load_cell_to_subclass(metadata_path, selected_subclasses)
    output_paths = {}
    for subclass in selected_subclasses:
        out_path = f"split.{fragment_path.replace('.tsv.gz', '')}.{subclass}.tsv"
        with open(out_path, "w"):
            pass
        output_paths[subclass] = out_path

    missing_rows = 0
    reader = pd.read_csv(
        fragment_path,
        sep="\t",
        header=None,
        names=FRAGMENT_COLUMNS,
        chunksize=CHUNK_SIZE,
        compression="gzip",
    )

    for chunk_idx, chunk in enumerate(reader, start=1):
        merged = chunk.merge(metadata, on="CellID", how="left")
        missing_mask = merged["SubClass"].isna()
        if missing_mask.any():
            missing_rows += int(missing_mask.sum())
            merged = merged.loc[~missing_mask]
        if not merged.empty:
            for subclass, group in merged.groupby("SubClass"):
                group_subset = group[FRAGMENT_COLUMNS]
                group_subset.to_csv(
                    output_paths[subclass],
                    mode="a",
                    sep="\t",
                    header=False,
                    index=False,
                )
        if chunk_idx % 5 == 0:
            print(f"Processed {chunk_idx} chunks (~{chunk_idx * CHUNK_SIZE:,} rows).")

    if missing_rows:
        print(f"Skipped {missing_rows:,} fragments without subclass mapping.")
    print("Done.")


if __name__ == "__main__":
    selected_subclasses = ['ITL23GL', 'ASC', 'CTGL', 'ITL6GL', 'PER', 'ITL5GL', 'PVGA',
       'NPGL', 'MGL', 'ITL4GL', 'OGC', 'PTGL', 'SSTGA', 'LAMGA', 'OPC',
       'VLMC', 'VIPGA', 'VEC', 'CLAGL', 'IOL', 'L6bGL', 'VPIA', 'STRGA',
       'MXD', 'D2MSN', 'OBGA1', 'MSGA', 'D1MSN', 'OBNBL', 'OLFGL', 'RGL',
       'LSXGA', 'CNUGA', 'PIRGL', 'OBDOP', 'CA3GL', 'OBGL', 'DGNBL',
       'CRC', 'OBGA2', 'ITHGL', 'CA1GL', 'GRC']
    split_fragments_by_subclass(
        fragment_path="Li2021/combine.mm10.frag.tsv.gz",
        metadata_path="Li2021/Li2021_sup2.tsv",
        selected_subclasses=selected_subclasses,
    )
