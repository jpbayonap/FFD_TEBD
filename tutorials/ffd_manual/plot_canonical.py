from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent


def find_default_csv():
    files = sorted(
        SCRIPT_DIR.glob("ffd_hm_fast_*.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if files:
        return files[0]

    fallback = SCRIPT_DIR / "ffd_hm_N40_order1.csv"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        f"No CSV found in {SCRIPT_DIR}. Expected ffd_hm_fast_*.csv"
    )


def load_data(csv_path):
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    required_cols = {
        "step",
        "t",
        "m",
        "h_real",
        "h_imag",
        "bh_real",
        "bh_imag",
        "max_discarded_fraction",
        "total_discarded_weight",
        "max_bond_dim",
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV is missing columns: {sorted(missing)}")

    if df.empty:
        raise ValueError(f"CSV file is empty: {csv_path}")

    return df


def save_h_density_plot(df, out_path):
    plt.figure(figsize=(7, 4))

    for m, g in df.groupby("m"):
        g = g.sort_values("t")
        plt.plot(g["t"], g["h_real"], label=f"m={m}")

    plt.xlabel("t")
    plt.ylabel(r"$\langle h_m \rangle$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def save_coupled_density_plot(df, out_path):
    plt.figure(figsize=(7, 4))

    for m, g in df.groupby("m"):
        g = g.sort_values("t")
        plt.plot(g["t"], g["bh_real"], label=f"m={m}")

    plt.xlabel("t")
    plt.ylabel(r"$\langle b_m h_m \rangle$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def save_discarded_fraction_plot(diag, out_path):
    plt.figure(figsize=(7, 4))
    plt.semilogy(diag["t"], diag["max_discarded_fraction"])
    plt.xlabel("t")
    plt.ylabel("max discarded fraction")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def save_total_discarded_weight_plot(diag, out_path):
    plt.figure(figsize=(7, 4))
    plt.semilogy(diag["t"], diag["total_discarded_weight"])
    plt.xlabel("t")
    plt.ylabel("total discarded weight per step")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def save_bond_dimension_plot(diag, out_path):
    plt.figure(figsize=(7, 4))
    plt.plot(diag["t"], diag["max_bond_dim"])
    plt.xlabel("t")
    plt.ylabel("max bond dimension")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def save_imaginary_part_plot(df, out_path):
    diag = (
        df.assign(abs_h_imag=df["h_imag"].abs())
        .groupby("step", as_index=False)
        .agg({"t": "first", "abs_h_imag": "max"})
    )

    plt.figure(figsize=(7, 4))
    plt.semilogy(diag["t"], diag["abs_h_imag"])
    plt.xlabel("t")
    plt.ylabel(r"max $|\mathrm{Im}\langle h_m\rangle|$")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="CSV file to plot. Defaults to newest ffd_hm_fast_*.csv.",
    )
    args = parser.parse_args()

    csv_path = args.csv if args.csv is not None else find_default_csv()
    if not csv_path.is_absolute():
        csv_path = (Path.cwd() / csv_path).resolve()

    df = load_data(csv_path)
    diag = df.drop_duplicates("step").sort_values("t")

    stem = csv_path.stem
    out_prefix = SCRIPT_DIR / stem

    outputs = {
        "h": out_prefix.with_name(f"{stem}_h.png"),
        "bh": out_prefix.with_name(f"{stem}_bh.png"),
        "discarded_fraction": out_prefix.with_name(f"{stem}_discarded_fraction.png"),
        "total_discarded_weight": out_prefix.with_name(
            f"{stem}_total_discarded_weight.png"
        ),
        "bond_dim": out_prefix.with_name(f"{stem}_bond_dim.png"),
        "imag": out_prefix.with_name(f"{stem}_imag.png"),
    }

    save_h_density_plot(df, outputs["h"])
    save_coupled_density_plot(df, outputs["bh"])
    save_discarded_fraction_plot(diag, outputs["discarded_fraction"])
    save_total_discarded_weight_plot(diag, outputs["total_discarded_weight"])
    save_bond_dimension_plot(diag, outputs["bond_dim"])
    save_imaginary_part_plot(df, outputs["imag"])

    print(f"loaded {csv_path}")
    for label, path in outputs.items():
        print(f"saved {label}: {path}")


if __name__ == "__main__":
    main()
