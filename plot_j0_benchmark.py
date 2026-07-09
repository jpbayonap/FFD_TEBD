import pandas as pd
import matplotlib.pyplot as plt

csv = "tutorials/ising_manual/j0_benchmark.csv"

df = pd.read_csv(csv)

# Max TEBD error vs exact for each dt and time
err = (
    df.groupby(["dt","t"])["abs_err_exact"].max().reset_index()
)

for dt, sub in err.groupby("dt"):
    plt.plot(sub["t"], sub["abs_err_exact"], label =f"dt={dt}")

plt.yscale("log")
plt.xlabel("t")
plt.ylabel("max_site |Z_TEBD - Z_exact|")
plt.legend()
plt.tight_layout()
plt.savefig("tutorials/ising_manual/j0_error_vs_time.png", dpi=200)
plt.show()