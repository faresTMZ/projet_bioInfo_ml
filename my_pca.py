import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


class MyPCA:
    """The PCA of the samples. 

    Attributes:
        __sequences: The RNA Sequences.
        __data: The input dataframe.
        __pca: The PCA object.
        __pca_data: The data after being scaled and going throught the PCA.
        __per_var: The PCA's explained_variance_ratio pourcentage.
        __pca_df: The new pca data in a dataframe form.
    """

    def __init__(self, seq, data, n_components=0.95):
        """Initializes the instance based on the data.

        Args:
            seq: The RNA Sequences.
            data: The input dataframe.
            n_components: The PCA's argument. 
        """
        self.__sequences = seq
        self.__data = data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(self.__data)
        self.__pca = PCA(n_components=n_components)
        self.__pca.fit(scaled_data)
        self.__pca_data = self.__pca.transform(scaled_data)

    def scree_plot(self):
        """Displays a scree plot of the performed PCA"""
        self.__per_var = np.round(self.__pca.explained_variance_ratio_ * 100, decimals=2)
        labels = ["PC" + str(x) for x in range(1, len(self.__per_var)+1, 2)]

        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        ax.bar(x=range(1,len(self.__per_var)+1), height=self.__per_var)
        ax.set_xticks(np.arange(len(self.__per_var), step=2) + 1)
        ax.set_xticklabels(labels, rotation=65)
        ax.set_ylabel("Percentage of Explained Variance")
        ax.set_xlabel("Principal Components")
        ax.set_title("Scree Plot")

        plt.show()
        plt.close()

    def visualize(self):
        """Plots several different PCs to visualize the PCA space."""
        labels = ["PC" + str(x) for x in range(1, len(self.__per_var)+1)]
        self.__pca_df = pd.DataFrame(self.__pca_data, index=self.__data.index, columns=labels)

        fig, axes = plt.subplots(3, 3, figsize=(18, 10))
        [[ax1, ax2, ax3], [ax4, ax5, ax6], [ax7, ax8, ax9]] = axes
        pc_labels = {n+1: self.__per_var[n] for n in range(7)}

        groups = self.__sequences.get_group_names()
        for group in groups:
            group_index = self.__sequences.get_group_counts(group).index
            pc =  {n+1: self.__pca_df.loc[group_index]["PC"+str(n+1)] for n in range(7)}

            def make_label(axe, n, m):
                axe.scatter(pc[n], pc[m])
                axe.set_xlabel(f"PC{n} - {pc_labels[n]}%")
                axe.set_ylabel(f"PC{m} - {pc_labels[m]}%")

            make_label(ax1, 1, 2)
            make_label(ax2, 1, 3)
            make_label(ax3, 1, 4)
            make_label(ax4, 1, 5)
            make_label(ax5, 1, 6)
            make_label(ax6, 2, 3)
            make_label(ax7, 2, 4)
            make_label(ax8, 2, 5)
            make_label(ax9, 3, 4)

        for axe in axes:
            for ax in axe:
                ax.legend(groups)

        fig.suptitle("PCA Graph of several PC combinations")

        plt.show()
        plt.close()

    def loading_scores(self):
        """Shows the top genes used for the first PC"""
        loading_scores = pd.Series(self.__pca.components_[0], index=self.__data.columns)
        sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
        top_10_genes = sorted_loading_scores[0:10].index.values

        loading_scores[top_10_genes]

