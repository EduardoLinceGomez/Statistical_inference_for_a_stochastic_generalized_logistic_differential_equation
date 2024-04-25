import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv(
    "output.csv", sep=',',
    header=1,
    names=[
        "i",
        "m_hat", "m_hat_aux",
        "alpha_hat", "alpha_hat_aux",
        "sigma_hat", "sigma_hat_aux"
    ]
)
df.describe()
df["sigma_hat_aux"].hist()
plt.show()
