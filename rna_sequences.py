import numpy as np
import pandas as pd
import glob
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from IPython.display import display


class RNASequences:
    """The rna counts and annotations of the samples. 

    Attributes:
        __rna_counts: A dataframe containing all the rna counts. Each row
            represents a sample and each column a gene.
        __annotations: A dataframe containing all samples annotations. Each row
            represents a sample and each column an annotation.
        __samples: The concatenation of __rna_counts and __annotations. Contains
            all the informations about each sample.
        __als_counts: A sub dataframe of __rna_counts containing all the samples
            of the "ALS Spectrum MND" group.
        __other_counts: A sub dataframe of __rna_counts containing all the samples
            of the "Other Neurological Disorders" group.
        __control_counts: A sub dataframe of __rna_counts containing all the samples
            of the "Non-Neurological Control group.
        __group_names: A list of string = ["ALS", "Other", "Control"]
    """

    def __init__(self, data_dir="Data/"):
        """Initializes the instance based on the data directory.

        Args:
          data_dir: The directory where the data is located.
        """
        # RNA Counts
        filenames = glob.glob("*.txt", root_dir=data_dir)
        dfs = []

        for filename in filenames:
            df = pd.read_csv(data_dir + filename,
                             sep="\t",
                             names=[filename[:10]],
                             skiprows=1).T
            dfs.append(df)
        
        self.__rna_counts = pd.concat(dfs)
        self.__remove_zeros()

        # Sample Annotations
        self.__annotations = pd.DataFrame(columns=["Subject ID",
                                                   "Sample Group",
                                                   "CNS Subregion"])
        tree = ET.parse(open("Data/GSE124439_family.xml"))
        root = tree.getroot()
        namespace = {"ns": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

        for sample in root.findall("ns:Sample", namespace):
            sample_id = sample.attrib["iid"]
            for channel in sample.iterfind(".//ns:Channel", namespace):
                self.__annotations.loc[sample_id] = [
                    channel[4].text.strip(),
                    channel[5].text.strip(),
                    channel[3].text.strip(),
                ]
        
        self.__check_annotations()

        self.__samples = pd.concat([self.__rna_counts, self.__annotations], axis=1)
        self.__als_counts = self.__get_group_counts("ALS Spectrum MND")
        self.__other_counts = self.__get_group_counts("Other Neurological Disorders")
        self.__control_counts = self.__get_group_counts("Non-Neurological Control")
        self.__group_names = ["ALS", "Other", "Control"]

    def __remove_zeros(self):
        """Removes columns with only zeros."""
        self.__rna_counts = (
            self.__rna_counts
            .loc[:, (self.__rna_counts != 0).any(axis=0)])
    
    def __check_annotations(self):
        """Checks if all the samples have annotations.

        Raises:
            AssertionError: If the indexes of __rna_counts and __annotations
                are different.
        """
        assert (
            self.__rna_counts
            .index
            .difference(self.__annotations.index)
            .empty
        )
    
    def __get_group_counts(self, group):
        """Returns the rna counts of the specified group.

        Args:
            group: A string representing the name of the group.
        """
        samples = self.__annotations[self.__annotations["Sample Group"] == group].index
        return self.__rna_counts.loc[samples]

    def _ipython_display_(self):
        """Displays the samples."""
        display(self.__samples)
    
    def get_samples(self):
        """Returns the rna counts and annotations of all samples."""
        return self.__samples

    def get_rna_counts(self):
        """Returns the rna counts of all samples."""
        return self.__rna_counts
    
    def get_annotations(self):
        """Returns the annotations of all samples."""
        return self.__annotations
    
    def get_group_counts(self, group=None):
        """Returns the counts of a specified test group.

        If there is no specified group returns the counts of all samples.

        Args:
            group: If not None, a string representing the name of the group.
                Groups are : "ALS", "Control" and "Other".
        """
        match group:
            case "ALS":
                return self.__als_counts
            case "Other":
                return self.__other_counts
            case "Control":
                return self.__control_counts
            case _:
                return self.__rna_counts
    
    def get_group_names(self):
        return self.__group_names
    
    def individual_samples_count(self):
        """Returns the amount of individuals per sample count."""
        return self.__annotations["Subject ID"].value_counts().value_counts()
    
    def groups_size(self, cumulative=False):
        """Returns the amount of samples per group."""
        sizes = self.__annotations["Sample Group"].value_counts()
        return np.cumsum(sizes) if cumulative else sizes
    
    def regions_count(self):
        """Returns the amount of samples per region."""
        return self.__annotations["CNS Subregion"].value_counts()
    
    def plot_annotations(self):
        """Plots in three pie plots all the sample's annotations."""
        fig, axes = plt.subplots(1, 3, figsize=(20, 5))

        [individual, group, region] = axes
        counts = [self.individual_samples_count(), self.groups_size(), self.regions_count()]

        for i, count in enumerate(counts):
            count.rename(None).plot(
                ax=axes[i],
                kind="pie",
                autopct=(lambda i: '{:.0f}'.format(i*count.sum()/100)),
            )

        individual.set_title("Amount of individuals per sample count")
        group.set_title("Amount of samples per group")
        region.set_title("Amount of samples per region")
        fig.suptitle("Samples description", fontsize=20)

        plt.show()
        plt.close()

    def mean(self, group=None, by_sample=False, sort=False, rename=True):
        """Returns the mean of rna counts of a group.

        Args:
            group: A string representing the name of the group. 
            by_sample: If False, computes the mean for each gene across all samples.
                If True, computes the mean for each sample across all genes.
            sort: If True, sorts the returned dataframe.
            rename: If True, renames the only column of the returned dataframe to "Mean".
        
        Returns:
            A dataframe with a single row of the computed means. 
        """
        mean = self.get_group_counts(group).mean(axis=by_sample).to_frame()

        if sort:
            mean = mean.sort_values(by=0)
        if rename:
            mean.rename(columns={0: "Mean"}, inplace=True)
        return mean
        
    def median(self, group=None, by_sample=False, sort=False, rename=True):
        """Returns the median of rna counts of a group.

        Args:
            group: A string representing the name of the group. 
            by_sample: If False, computes the median for each gene across all samples.
                If True, computes the median for each sample across all genes.
            sort: If True, sorts the returned dataframe.
            rename: If True, renames the only column of the returned dataframe to "Median".
        
        Returns:
            A dataframe with a single row of the computed means. 
        """
        median = self.get_group_counts(group).median(axis=by_sample).to_frame()

        if sort:
            median = median.sort_values(by=0)
        if rename:
            median.rename(columns={0: "Median"}, inplace=True)
        return median
        
    def std(self, group=None, by_sample=False, sort=False, rename=True):
        """Returns the std of rna counts of a group.

        Args:
            group: A string representing the name of the group. 
            by_sample: If False, computes the std for each gene across all samples.
                If True, computes the std for each sample across all genes.
            sort: If True, sorts the returned dataframe.
            rename: If True, renames the only column of the returned dataframe to "Std".
        
        Returns:
            A dataframe with a single row of the computed means. 
        """
        std = self.get_group_counts(group).std(axis=by_sample).to_frame()

        if sort:
            std = std.sort_values(by=0)
        if rename:
            std.rename(columns={0: "Std"}, inplace=True)
        return std

    def plot_metrics(self, by_sample=False, log=False):
        """Plots the mean, median & std of each group."""
        np.seterr(divide = "ignore")
        fig, axes = plt.subplots(1, 3, figsize=(20, 5))
        [mean, median, std] = axes

        groups = self.get_group_names()
        for group in groups:
            (self.mean(group=group, by_sample=by_sample, sort=True)
             .apply(lambda x: np.log10(x) if log else x)
             .plot(ax=mean))
            (self.median(group=group, by_sample=by_sample, sort=True)
             .apply(lambda x: np.log10(x) if log else x)
             .plot(ax=median))
            (self.std(group=group, by_sample=by_sample, sort=True)
             .apply(lambda x: np.log10(x) if log else x)
             .plot(ax=std))

        for axe in axes:
            if by_sample:
                x_label = "Samples"
                stop = max(self.groups_size())
                x_tick_labels = np.linspace(start=0, stop=stop, num=10, dtype=int)
                title = "Mean, median & std of each sample across all genes"
            else:
                x_label = "Genes"
                stop = len(self.get_rna_counts().columns)
                x_tick_labels = np.linspace(start=0, stop=stop, num=8, dtype=int)
                title = "Mean, median & std of each gene across all samples"
            
            axe.set_xlabel(x_label, fontsize=15)
            axe.set_xticks(x_tick_labels)
            axe.set_xticklabels(x_tick_labels)
            axe.legend(groups)
            fig.suptitle(title + " (log)" if log else title, fontsize=20)

        mean.set_ylabel("Mean", fontsize=15)
        median.set_ylabel("Median", fontsize=15)
        std.set_ylabel("Std", fontsize=15)

        plt.show()
        plt.close()
        np.seterr(divide = "warn")

    def plot_rsd(self):
        """Plots the rsd of each group."""
        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        
        groups = self.get_group_names()
        for group in groups:
            mean = self.mean(group=group, rename=False)
            std = self.std(group=group, rename=False)
            rsd = std.div(mean).sort_values(by=0).rename(columns={0: "RSD"}).apply(lambda x: np.log10(x))
            rsd.plot(ax=ax)

        ax.set_xlabel("Genes", fontsize=15)
        
        stop = len(self.get_rna_counts().columns)
        x_tick_labels = np.linspace(start=0, stop=stop, num=15, dtype=int)
        ax.set_xticks(x_tick_labels)
        ax.set_xticklabels(x_tick_labels)
        ax.legend(groups)
        ax.set_ylabel("RSD", fontsize=15)
        fig.suptitle("Relative Standard Deviation of each gene across all samples (log)", fontsize=20)

        plt.show()
        plt.close()

