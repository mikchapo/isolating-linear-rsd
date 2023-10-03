import matplotlib.pyplot as plt
import numpy as np


fit_stats = np.loadtxt("../output/data_products/eboss_fiducial_fit_stats.dat",
                       dtype={'names': ('label', 'mean', 'std', 'bestfit', 'color'),
                              'formats': ('S12', 'f4', 'f4', 'f4', 'S1')})

fig_width = 5.
plt.figure(figsize=(fig_width, fig_width), dpi=300)
plt.xlim(-2.5, 2.5)
for i in range(-2, 3):
    if i == 0:
        plt.axvline(i, color="k", linestyle="-")
    else:
        plt.axvline(i, color="k", linestyle=":")
for i in range(fit_stats.shape[0]):
    print(fit_stats[i]["label"])
    plt.plot((fit_stats[i]["bestfit"] - fit_stats[i]["mean"]) /
             fit_stats[i]["std"], i, marker="o", color=fit_stats[i]["color"].decode("utf-8"))
    plt.text(-2.7, i, fit_stats[i]["label"].decode("utf-8"), ha="right",
             va="center")
plt.subplots_adjust(left=0.35, top=0.99, bottom=0.1)
plt.yticks([])
plt.savefig("../output/plots/fit_stats.jpg")
