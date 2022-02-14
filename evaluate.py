"""Evaluate simple classifiers on enzyme class based on DRFP."""

import pickle
from pathlib import Path

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.datasets import load_iris, make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.linear_model import ElasticNet, LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import LinearSVC
from sklearn.tree import DecisionTreeClassifier

from utils import LABELS, OUTPUT, UNKNOWN_LABEL


@click.command()
def main(version: str = "120"):
    """Train and evaluate several simple ML models for classification of top-level EC code."""
    version_dir = OUTPUT.joinpath(version)
    models_dir = version_dir.joinpath("models")
    models_dir.mkdir(exist_ok=True, parents=True)
    metadata_path = version_dir.joinpath("reaction_metadata.tsv")
    fingerprints_path = version_dir.joinpath("reaction_fingerprints.tsv.gz")

    fingerprints = pd.read_csv(fingerprints_path, index_col=0, sep="\t")
    metadata = pd.read_csv(metadata_path, index_col=0, sep="\t")
    idx = metadata["enzyme_class"] != UNKNOWN_LABEL
    X, y = fingerprints[idx], metadata[idx].enzyme_class

    print(metadata[idx].enzyme_class.value_counts())

    rows = []
    classifiers = [
        MLPClassifier(),
        DecisionTreeClassifier(),
        KNeighborsClassifier(n_neighbors=5),
        LogisticRegression(multi_class="multinomial", max_iter=500),
        RandomForestClassifier(),
    ]
    for clf in classifiers:
        name = clf.__class__.__name__.removesuffix("Classifier")
        print("running", name)
        scores = cross_val_score(
            clf, X, y, scoring="roc_auc_ovr", cv=StratifiedKFold(5)
        )
        rows.extend(((name, score) for score in scores))
        models_dir.joinpath(name).with_suffix(".pkl").write_bytes(pickle.dumps(clf))

    df = pd.DataFrame(rows, columns=["classifier", "score"])
    df.to_csv(models_dir.joinpath("evaluation.tsv"), sep="\t", index=False)

    fig, ax = plt.subplots()
    sns.boxplot(data=df, y="classifier", x="score", ax=ax)
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(models_dir.joinpath("clf_results.png"), dpi=300)

    # only predict on LR reactions to reduce duplication
    pred_idx = (~idx) & (metadata["direction"] == "lr")
    results = RandomForestClassifier().fit(X, y).predict(fingerprints[pred_idx])
    prediction_df = metadata.loc[pred_idx, ["smiles"]]
    prediction_df["predicted_class"] = results
    prediction_df["predicted_class_label"] = prediction_df["predicted_class"].map(
        LABELS
    )
    prediction_df = prediction_df[
        ["predicted_class", "predicted_class_label", "smiles"]
    ]
    prediction_df.to_csv(models_dir.joinpath("predictions.tsv"), sep="\t", index=True)


if __name__ == "__main__":
    main()
