import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

path = "../data/"
my_dtypes = {"i": np.int32, "t_i": np.float64, "w(t_i)": np.float64}
df_op_sample_path = pd.read_csv(
    path + "sample_path.csv",
    header=0,
    names=["i", "t_i", "w(t_i)"],
    sep=",",
    dtype=my_dtypes,
)
df_res_sample_path = pd.read_csv(
    path + "res_sample_path.csv",
    sep=",",
    dtype=my_dtypes,
    header=0,
    names=["i", "t_i", "w(t_i)"],
)
#
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
ax.plot(
    df_res_sample_path["t_i"],
    df_res_sample_path["w(t_i)"],
    **marker_style_00,
    label=r"$W(t_i)$"
)
ax.plot(
    df_op_sample_path["t_i"],
    df_op_sample_path["w(t_i)"],
    **marker_style,
    label=r"$W(\tau_i), \quad \tau=r\delta t$"
)
plt.legend(loc=0)
plt.xlabel(r"$t$")
plt.ylabel(r"$W(t)$")
plt.savefig("brownian_sample_likening.png", dpi=300)
# plt.show()
