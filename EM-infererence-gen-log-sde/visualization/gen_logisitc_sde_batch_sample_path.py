import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

path = "../data/"
dir_path = r"../data"
dir = os.listdir(dir_path)

my_dtypes = {"idx": np.int32, "i": np.int32, "t_i": np.float64, "X(t_i)": np.float64}

df_batch_path = pd.read_csv(
    path + dir[-1],
    header=0,
    names=["idx", "i", "t_i", "X(t_i)"],
    sep=",",
    # dtype=my_dtypes
)
fig, ax = plt.subplots()
idx = pd.unique(df_batch_path["idx"])
marker_style_00 = dict(
    color="gray",
    linestyle="-",
    # marker="",
    markersize=1,
    markerfacecoloralt="gray",
    alpha=0.1,
)
x_T = []
for j in idx:
    path_j = df_batch_path[df_batch_path["idx"] == j]
    ax.plot(path_j["t_i"], path_j["X(t_i)"], **marker_style_00, label=r"$X(t_i)$")
    x_T.append(path_j["X(t_i)"].iat[-1])

plt.xlabel(r"$t$")
plt.ylabel(r"$X(t)$")
x_T = np.array(x_T)
ymax = np.max(np.abs(x_T))
binwidth = 0.025
lim = (int(ymax / binwidth) + 1) * binwidth
bins = 100  # np.arange(0, lim + binwidth, binwidth)
divider = make_axes_locatable(ax)
ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)
n, bins, patches = ax_histy.hist(x_T, bins=bins, orientation="horizontal", density=True)
sigma = np.std(x_T)
mu = np.mean(x_T)
y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(
    -0.5 * ((1.0 / sigma) * (bins - mu)) ** 2
)
# ax_histy.plot(y, bins, "--")
# plt.legend(loc=0)
plt.savefig("gen_log_sde_batch_sample_path.png", dpi=300)
g = grid = sns.JointGrid(data=df_batch_path, x="t_i", y="X(t_i)")
g.plot_joint(sns.scatterplot)
g.plot_marginals(sns.kdeplot)
g.plot_marginals(sns.histplot)
g.ax_marg_x.remove()
# g = sns.jointplot(data=df_batch_path, x="t_i", y="X(t_i)", marginal_ticks=True)
# sns.kdeplot(df_batch_path["X(t_i)"], ax=g.ax_marg_x, legend=False)
g.savefig("marginal.png")
# plt.show()
