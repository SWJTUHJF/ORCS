import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

ue_tstt, so_tstt = 7480225, 7194922
tstt_01 = [7269329.409151747, 7207530.281915506, 7205796.696519641, 7195742.307461868, 7195742.307461868]
control_ratio_01 = [0.3596783139212424, 0.413477537437604, 0.42252194074818106, 0.42834556803603463, 0.42834556803603463]

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(range(1, 6), tstt_01)
axes[0].set_ylabel("Total travel time", fontsize=12)
axes[0].set_xlabel("Number of iteration", fontsize=12)
axes[0].set_xticks([1, 2, 3, 4, 5])
axes[0].hlines(ue_tstt, xmax=5, xmin=1, color="#FFA500", linestyles="--", label="UE")
axes[0].hlines(so_tstt, xmax=5, xmin=1, color="green", linestyles="--", label="SO")
axes[0].legend(bbox_to_anchor=(0.97, 0.9))
axes[1].plot(range(1, 6), control_ratio_01)
axes[1].set_ylabel("% of SO users", fontsize=12)
axes[1].set_xlabel("Number of iteration", fontsize=12)
axes[1].set_xticks([1, 2, 3, 4, 5])
axes[1].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
plt.show()