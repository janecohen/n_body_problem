#%%
import numpy as np
from scipy.integrate._ivp import dop853_coefficients

A = dop853_coefficients.A[:12, :12]
B = dop853_coefficients.B

with open("dop853_coefficients.txt", "w") as f:
    for row in A:
        f.write(" ".join(f"{x:.17e}" for x in row) + "\n")

    f.write(" ".join(f"{x:.17e}" for x in B) + "\n")
# %%