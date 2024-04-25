import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(
    "../data/realization_path_gen_logistic_SDE.csv", sep=',',
    header=1,
    names=[
        "i",
        "t_i", "x_milstein(t_i)"
    ]
)
path_line_plot = sns.lineplot(
    data=df,
    x="t_i",
    y="x_milstein(t_i)"
)
path_line_plot.axhline(
    y=1.0,
    color='r',
    ls='--',
    label='Carrying capacity'
)
plt.show()
