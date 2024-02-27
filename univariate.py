import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


class Ttest:
    """The univariate analysis using T-test. 

    Attributes:
        __sequences: The RNA sequences.
        __rna_counts: The RNA counts.
        __annotations: The samples annotations.
        __p_values: The p-values of all genes (ALS vs Control).
        __q_values: The corrected p-values (using fdr_bh).
        __fold_change: The fold change of Control over ALS.
    """

    def __init__(self, seq):
        """Initializes the instance based on the data.

        Args:
            seq: The RNA Sequences.
        """
        self.__sequences = seq
        self.__rna_counts = self.__sequences.get_rna_counts()
        self.__annotations = self.__sequences.get_annotations()
        als_index = self.__annotations["Sample Group"] == "ALS Spectrum MND"
        ctrl_index = self.__annotations["Sample Group"] == "Non-Neurological Control"

        self.__p_values = self.__rna_counts.apply(lambda x: ttest_ind(x[als_index.values], x[ctrl_index.values]).pvalue).fillna(1)
        _, self.__q_values, _, _ = multipletests(self.__p_values.values, method="fdr_bh")

        ctrl_counts = self.__sequences.mean(group="Control")
        als_counts = self.__sequences.mean(group="ALS")
        self.__fold_change = ctrl_counts.div(als_counts).squeeze().apply(lambda x: x if (x == 0) else np.log2(x))

    def plot_metrics(self):
        """Plots the p-values (raw & corrected) and the fold change of each gene."""
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        [ax1, ax2] = axes

        self.__p_values.sort_values().plot(ax=ax1)
        pd.Series(self.__q_values, index=self.__p_values.index).sort_values().plot(ax=ax1)
        ax1.set_ylabel("p-value")
        ax1.legend(["raw", "corrected (fdr_bh)"])

        self.__fold_change.sort_values().plot(ax=ax2)
        ax2.set_ylabel("fold change (Control/ALS)")
        ax2.legend(["fold change"])

        for axe in axes:
            axe.set_xlabel("Genes")
            stop = len(self.__sequences.get_rna_counts().columns)
            x_tick_labels = np.linspace(start=0, stop=stop, num=8, dtype=int)
            axe.set_xticks(x_tick_labels)
            axe.set_xticklabels(x_tick_labels)

        title = "P-values (raw & corrected) and fold change of each gene across all ALS and Control groups"
        fig.suptitle(title)

        plt.show()
        plt.close()

    def volcano_plot(self):
        """Displays a volcano plot of the performed T-test"""
        fig, ax = plt.subplots(1, 1)
        rev_q_values = -np.log10(self.__q_values)

        ax.scatter(self.__fold_change, rev_q_values)
        ax.axhline(y=0.5, color='r', linestyle='--')

        low_pvals = self.__p_values[-np.log10(self.__q_values) > 0.5].index
        for gene in low_pvals:
            plt.annotate(gene, (self.__fold_change[gene], rev_q_values[self.__p_values.index.get_loc(gene)]))

        ax.set_xlabel("log2 fold change ")
        ax.set_ylabel("-log10 p-value")
        fig.suptitle("Volcano Plot")

        plt.show()
        plt.close()