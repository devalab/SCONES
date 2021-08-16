import os
import sys
import argparse

import numpy as np
import pandas as pd
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("darkgrid")

NEUTRAL_MAGNITUDE=0.5

def analyze_dataset_distribution(name: str, dataset: pd.DataFrame, fieldname, sign_flip: bool):
    ddG = np.asarray(dataset[fieldname])
    if sign_flip:
        ddG *= -1

    ddG = ddG[~np.isnan(ddG)]

    ddG_negative_data = ddG[ddG < -NEUTRAL_MAGNITUDE]
    ddG_neutral_data = ddG[(-NEUTRAL_MAGNITUDE <= ddG) & (ddG <= NEUTRAL_MAGNITUDE)]
    ddG_positive_data = ddG[ddG > NEUTRAL_MAGNITUDE]

    print("ddG >> number of samples: %d, mean: %.2f, stddev: %.2f" % (len(ddG), np.mean(ddG), np.std(ddG)))
    print("ddG negative >> number of samples: %d, mean: %.2f, stddev: %.2f" %
        (len(ddG_negative_data), np.mean(ddG_negative_data), np.std(ddG_negative_data)))
    print("ddG neutral >> number of samples: %d, mean: %.2f, stddev: %.2f" %
        (len(ddG_neutral_data), np.mean(ddG_neutral_data), np.std(ddG_neutral_data)))
    print("ddG positive >> number of samples: %d, mean: %.2f, stddev: %.2f" %
        (len(ddG_positive_data), np.mean(ddG_positive_data), np.std(ddG_positive_data)))

    kde = stats.gaussian_kde(ddG)
    negative_prob = kde.integrate_box_1d(-np.inf, -NEUTRAL_MAGNITUDE)
    neutral_prob = kde.integrate_box_1d(-NEUTRAL_MAGNITUDE, NEUTRAL_MAGNITUDE)
    positive_prob = kde.integrate_box_1d(NEUTRAL_MAGNITUDE, np.inf)

    print("Percentage signficantly negative samples: %.2f%%" % (negative_prob * 100))
    print("Percentage neutral samples: %.2f%%" % (neutral_prob * 100))
    print("Percentage signficantly positive samples: %.2f%%" % (positive_prob * 100))

    g = sns.displot(ddG)
    ax = g.facet_axis(0, 0)
    ax.set_xlabel(r"$\Delta \Delta G$ (kcal/mol)", fontsize=15)
    ax.set_ylabel("Count", fontsize=15)
    ax.set_xticklabels([int(value) for value in ax.get_xticks()], fontsize=12.5)
    ax.set_yticklabels([int(value) for value in ax.get_yticks()], fontsize=12.5)

    plt.title("distribution of ddG for %s" % name)
    plt.tight_layout()
    plt.show()

def analyze_dataset_mutationwise_composition(name: str, dataset: pd.DataFrame):
    letters = list("ACDEFGHIKLMNPQRSTVWY")
    counts = pd.DataFrame(0, index=letters, columns=letters)
    for _, sample in dataset.iterrows():
        ref_residue, mutated_residue = sample["ref_residue"], sample["mutated_residue"]
        counts[ref_residue][mutated_residue] += 1

    size = 10
    plt.figure(figsize=(size, size))
    ax1 = plt.subplot2grid((size, size), (0, 0), colspan=size-1, rowspan=size-1)
    ax2 = plt.subplot2grid((size, size), (size - 1,0), colspan=size-1, rowspan=1)
    ax3 = plt.subplot2grid((size, size), (0, size - 1), colspan=1, rowspan=size-1)

    sns.heatmap(counts, ax=ax1, annot=True, xticklabels=True, yticklabels=True, fmt='g', cbar=False)
    ax1.xaxis.tick_top()
    ax1.set_ylabel("Mutated residue", fontsize=15)
    ax1.set_title("mutation matrix for %s" % name, fontsize=20)

    colsum = pd.DataFrame(counts.sum(axis=0))
    sns.heatmap(colsum.transpose(), ax=ax2, annot=True, xticklabels=True, yticklabels=False, fmt='g', cbar=False)
    ax2.xaxis.tick_top()
    ax2.set_xlabel("Ref. residue", fontsize=15)
    ax2.set_ylabel("colsum")

    rowsum = pd.DataFrame(counts.sum(axis=1))
    sns.heatmap(rowsum, ax=ax3,  annot=True, xticklabels=False, yticklabels=True, fmt='g', cbar=False)
    ax3.set_xlabel("rowsum")
    ax3.set_ylabel("")

    plt.tight_layout()
    plt.show()
    return counts

def analyze_dataset_frenzcategories_composition(name: str, dataset: pd.DataFrame):
    aa_categories = {
        "small" : ['G', 'A', 'V', 'S', 'T', 'C'],
        "large" : ['F', 'Y', 'W', 'K', 'R', 'H', 'Q', 'E'],
        "negative" : ['D', 'E'],
        "positive" : ['R', 'K'],
        "polar" : ['Y', 'T', 'S', 'H', 'K', 'R', 'E', 'D', 'Q', 'N'],
        "non-charged polar" : ['Y', 'T', 'S', 'N', 'Q', 'H'],
        "hydrophobic" : ['F', 'I', 'L', 'V', 'A', 'G', 'M', 'W'],
        "cysteine" : ['C'],
        "proline" : ['P']
    }

    mutation_categories = {
        "positive to negative" : [("positive", "negative")],
        "negative to positive" : [("negative", "positive")],
        "positive to non-charged polar" : [("positive", "non-charged polar")],
        "negative to non-charged polar" : [("negative", "non-charged polar")],
        "non-charged polar to positive" : [("non-charged polar", "positive")],
        "non-charged polar to negative" : [("non-charged polar", "negative")],
        "negative to hydrophobic" : [("negative", "hydrophobic")],
        "hydrophobic to negative" : [("hydrophobic", "negative")],
        "positive to hydrophobic" : [("positive", "hydrophobic")],
        "hydrophobic to positive" : [("hydrophobic", "positive")],
        "like to like charge" : [("positive", "positive"), ("negative", "negative")],
        "non-charged polar to hydrophobic" : [("non-charged polar", "hydrophobic")],
        "hydrophobic to non-charged polar" : [("hydrophobic", "non-charged polar")],
        "non-charged polar to non-charged polar" : [("non-charged polar", "non-charged polar")],
        "hydrophobic to hydrophobic" : [("hydrophobic", "hydrophobic")],
        "involves proline" : [("proline", ""), ("", "proline")], # "" category represents wildcard matching
        "involves cysteine" : [("cysteine", ""), ("", "cysteine")],
        "small to large" : [("small", "large")],
        "large to small": [("large", "small")],
        "same to same size" : [("small", "small"), ("large", "large")]
    }

    df = pd.DataFrame(0, index=mutation_categories.keys(), columns=["num_samples"])
    for _, sample in dataset.iterrows():
        ref_residue, mutated_residue = sample["ref_residue"], sample["mutated_residue"]
        ref_categories = [key for key, residues in aa_categories.items() if ref_residue in residues]
        mutated_categories = [key for key, residues in aa_categories.items() if mutated_residue in residues]

        sample_categories = set()
        for catname, groupings in mutation_categories.items():
            for ref_group, mutated_group in groupings:
                if (ref_group in ref_categories or ref_group == "") and (mutated_group in mutated_categories or mutated_group == ""):
                    sample_categories.add(catname)

        for category in sample_categories:
            df.at[category, "num_samples"] += 1

    df.reset_index(inplace=True)

    g = sns.catplot(x='index', y='num_samples', kind='bar', data=df, height=5, aspect=1.25, palette=sns.color_palette(["blue"]))
    g.set_xticklabels(rotation=90)
    ax = g.facet_axis(0,0)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, int(ymax * 1.1))
    for p in ax.patches:
        ax.text(p.get_x() + 0.4,
                p.get_height() + 0.5,
                int(p.get_height()),
                color='black', rotation='horizontal', ha='center')
    ax.set_xticklabels(mutation_categories.keys(), fontsize=12.5)
    ax.set_yticklabels([int(value) for value in ax.get_yticks()], fontsize=12.5)
    ax.set_xlabel("")
    ax.set_ylabel("Count", fontsize=15)
    ax.set_title("counts of samples in different categories (dataset: %s)" % name)
    g.tight_layout()
    plt.show()
    return df

def main():
    parser = argparse.ArgumentParser(description='Analyze mutation composition and ddG distribution of a dataset')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to dataset in CSV', required=True, default=None)
    parser.add_argument('--ddG_fieldname', dest='ddG_fieldname', type=str, help='column that contains ddG values', default="ddG")
    parser.add_argument('--sign_flip', dest='sign_flip', type=bool, help='reverse ddG signs', default=False)
    parser.add_argument('--name', dest='name', type=str, help='name of the dataset', default=None)

    args = parser.parse_args()
    if args.dataset_path is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    ddG_fieldname = args.ddG_fieldname
    sign_flip = args.sign_flip
    name = args.name

    if name is None:
        basename = os.path.basename(dataset_path)
        name = os.path.splitext(basename)[0]

    dataset = pd.read_csv(dataset_path)
    analyze_dataset_distribution(name, dataset, ddG_fieldname, sign_flip)

    print("Mutation Matrix")
    print(analyze_dataset_mutationwise_composition(name, dataset))
    print()

    print("Composition")
    print(analyze_dataset_frenzcategories_composition(name, dataset))

if __name__ == "__main__":
    main()
