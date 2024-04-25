import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

path = "../data/"
my_dtypes = {"i": np.int32, "t_i": np.float64, "X(t_i)": np.float64}
df_sde = pd.read_csv(
    path + "realization_path_gen_logistic_SDE.csv",
    header=0,
    names=["i", "t_i", "X(t_i)"],
    sep=",",
    # dtype=my_dtypes
)

fig, ax = plt.subplots()
marker_style = dict(
    color="cornflowerblue",
    linestyle=":",
    marker="o",
    markersize=5,
    markerfacecoloralt="gray",
    alpha=0.5,
)
marker_style_00 = dict(
    color="orange",
    linestyle=":",
    marker=".",
    markersize=1,
    markerfacecoloralt="gray",
    alpha=0.7,
)
ax.plot(df_sde["t_i"], df_sde["X(t_i)"], **marker_style_00, label=r"$X(t_i)$")
plt.legend(loc=0)
plt.xlabel(r"$t$")
plt.ylabel(r"$X(t)$")
plt.savefig("gen_log_sde_sample_path.png", dpi=300)
# plt.show()
