"""  """

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker


def create_pfm(motifs):
    """ constructs a position frequency matrix from a list of binding sites """
    pfm_dict = {"A": [], "C": [], "G": [], "T": []}

    for j in range(len(motifs[0])):
        pfm_dict["A"].append(([i[j] for i in motifs].count("A")))
        pfm_dict["C"].append(([i[j] for i in motifs].count("C")))
        pfm_dict["G"].append(([i[j] for i in motifs].count("G")))
        pfm_dict["T"].append(([i[j] for i in motifs].count("T")))
    return pfm_dict


def get_ic(pfm):
    """ Create information content profiles for Shannon entropy """
    shan_ent = {"A": [], "C": [], "G": [], "T": []}

    df = pd.DataFrame.from_dict(pfm)
    df["sum"] = df.sum(axis=1)
    ppm = df.loc[:, "A":"T"].div(df["sum"], axis=0)

    for index, row in ppm.iterrows():
        # SeqLogo information content profile (Shannon uncertainty measure)
        entropy = 2 + (row["A"] * (0 if (row["A"]) == 0 else np.log2(row["A"]))) + \
                  (row["C"] * (0 if (row["C"]) == 0 else np.log2(row["C"]))) + \
                  (row["G"] * (0 if (row["G"]) == 0 else np.log2(row["G"]))) + \
                  (row["T"] * (0 if (row["T"]) == 0 else np.log2(row["T"])))

        shan_ent["A"].append(row["A"] * entropy)
        shan_ent["C"].append(row["C"] * entropy)
        shan_ent["G"].append(row["G"] * entropy)
        shan_ent["T"].append(row["T"] * entropy)

    return pd.DataFrame.from_dict(shan_ent)


def create_logo(ic, reg_name, out_dir):
    """ Creates sequence logo in pdf format """
    pdf_file = os.path.join(out_dir, f"{reg_name}_logo.pdf")
    if not os.path.exists(pdf_file):
        plt.figure(figsize=(90, 40))

        logo = logomaker.Logo(ic, alpha=.75)
        logo.style_spines(visible=False)
        logo.style_spines(spines=('left', 'bottom'), visible=True)
        logo.style_xticks(rotation=0, fmt='%d', anchor=0)

        logo.ax.set_ylabel("Information content", labelpad=-1)
        logo.ax.xaxis.set_ticks_position('none')
        logo.ax.xaxis.set_tick_params(pad=-1)
        logo.ax.set_xticklabels([x + 1 for x in ic.index.values])
        logo.ax.set_ylim([0, 2])

        plt.savefig(pdf_file, format="pdf")
    return
