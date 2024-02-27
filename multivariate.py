import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import matthews_corrcoef, make_scorer


class ElasticNet:
    """The multivariate analysis using Elastic-Net. 

    Attributes:
        __rna_counts: The RNA Counts.
        __scaled_counts: The scaled RNA Counts.
        __annotations: The samples annotations.
        __elasticNet: The fine tuned Elastic-Net model.
    """

    def __init__(self, seq):
        """Initializes the instance based on the data.

        Args:
            seq: The RNA Sequences.
        """
        self.__rna_counts = seq.get_rna_counts()
        self.__annotations = seq.get_annotations()
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(self.__rna_counts)
        self.__scaled_counts = pd.DataFrame(scaled_data, columns=self.__rna_counts.columns, index=self.__rna_counts.index)

        X = self.__scaled_counts.to_numpy()
        y = pd.concat([self.__scaled_counts, self.__annotations], axis=1)["Sample Group"].to_numpy()

        self.__elasticNet = LogisticRegressionCV(penalty="elasticnet",
                                                 cv=3,
                                                 solver="saga",
                                                 l1_ratios=[0.25, 0.5, 0.75],
                                                 Cs=[0.1, 0.5],
                                                 scoring=make_scorer(matthews_corrcoef)
        ).fit(X, y)

        model_prediction = self.__elasticNet.predict(X)
        mcc_train_dataset = matthews_corrcoef(y, model_prediction)
        print(f"MCC on train set : {mcc_train_dataset}")

    def save_top_100_genes(self):
        """Saves the top 100 genes to a csv file."""
        coeffs = np.abs(self.__elasticNet.coef_).mean(axis=0)
        top_genes = pd.DataFrame(coeffs, self.__scaled_counts.columns, columns=["Coefficients"])
        top_100_genes = top_genes.sort_values(by="Coefficients")[-100:]
        top_100_genes.to_csv("top_100_genes.csv", sep='\t', float_format='%f', header=False)