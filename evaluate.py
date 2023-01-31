"""Evaluate simple classifiers on enzyme class based on DRFP."""

import pickle

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from more_click import force_option
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from tqdm import tqdm

from utils import LABELS, OUTPUT, UNKNOWN_LABEL, random_state_option, version_option


@click.command()
@version_option
@random_state_option
@force_option
def main(version: str, random_state: int, force: bool):
    """Train and evaluate several simple ML models for classification of top-level EC code."""
    click.secho(f"Using Rhea v{version} and random state of {random_state}", fg="green")
    version_dir = OUTPUT.joinpath(version)
    models_dir = version_dir.joinpath("models")
    models_dir.mkdir(exist_ok=True, parents=True)

    results_path = models_dir.joinpath("predictions.tsv")
    if results_path.is_file() and not force:
        click.echo("::set-output name=RHEA_LEARNED::false")
    else:
        click.echo(f"::set-output name=RHEA_LEARNED::{version}")

        metadata_path = version_dir.joinpath("reaction_metadata.tsv")
        fingerprints_path = version_dir.joinpath("reaction_fingerprints.tsv.gz")

        fingerprints = pd.read_csv(fingerprints_path, index_col=0, sep="\t")
        metadata = pd.read_csv(metadata_path, index_col=0, sep="\t")
        idx = metadata["enzyme_class"] != UNKNOWN_LABEL
        X, y = fingerprints[idx], metadata[idx].enzyme_class

        click.secho("Enzyme Class counts", fg="green")
        print(metadata[idx].enzyme_class.value_counts())

        rows = []
        classifiers = [
            LogisticRegression(
                multi_class="multinomial", max_iter=500, random_state=random_state
            ),
            KNeighborsClassifier(n_neighbors=5),  # no random state necessary
            DecisionTreeClassifier(random_state=random_state),
            RandomForestClassifier(random_state=random_state),
            MLPClassifier(random_state=random_state),
        ]
        it = tqdm(classifiers)
        for clf in it:
            name = clf.__class__.__name__.removesuffix("Classifier")
            it.write(f"running {name}")
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
        prediction_df.to_csv(results_path, sep="\t", index=True)


if __name__ == "__main__":
    main()
