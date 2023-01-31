import warnings

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from drfp import DrfpEncoder
from jinja2 import Environment, FileSystemLoader
from more_click import force_option, verbose_option
from rdkit.Chem import AllChem, rdChemReactions
from scipy.sparse import SparseEfficiencyWarning
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
from tqdm import tqdm

from utils import (
    LABELS,
    MODULE,
    OUTPUT,
    README_PATH,
    TEMPLATES,
    UNKNOWN_LABEL,
    random_state_option,
    version_option,
)

warnings.simplefilter("ignore", SparseEfficiencyWarning)

environment = Environment(
    autoescape=True, loader=FileSystemLoader(TEMPLATES), trim_blocks=False
)
index_template = environment.get_template("README.md.jinja")


@click.command()
@verbose_option
@force_option
@random_state_option
@version_option
@click.option("--progress", is_flag=True)
def main(force: bool, random_state: int, version: str, progress: bool):
    """Generate SMILES and differential reaction fingerprints for each Rhea entry."""
    click.secho(f"Using Rhea v{version} and random state of {random_state}", fg="green")
    output = OUTPUT.joinpath(version)
    output.mkdir(exist_ok=True, parents=True)
    metadata_path = output.joinpath("reaction_metadata.tsv")
    fingerprints_path = output.joinpath("reaction_fingerprints.tsv.gz")
    fingerprints_2d_path = output.joinpath("reaction_fingerprints_2d.tsv")

    if metadata_path.is_file() and fingerprints_path.is_file() and not force:
        metadata_df = pd.read_csv(metadata_path, sep="\t", index_col=0)
        fingerprint_df = pd.read_csv(fingerprints_path, sep="\t", index_col=0)
        click.echo("::set-output name=RHEA_UPDATED::false")
    else:
        click.echo(f"::set-output name=RHEA_UPDATED::{version}")
        # Use a version-specific PyStow directory for reproducible downloading
        # of external resources
        tqdm.write("Downloading/unpacking reaction archive")
        directory = MODULE.ensure_untar(
            version,
            url="ftp://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz",
        )
        tqdm.write("Downloaded/unpacking reaction archive")

        # Columns are: RHEA_ID_MASTER, RHEA_ID_LR, RHEA_ID_RL, RHEA_ID_BI
        directions_df = MODULE.ensure_csv(
            version,
            url="https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv",
            read_csv_kwargs=dict(dtype=str),
        )
        direction_master = set(directions_df["RHEA_ID_MASTER"])
        direction_lr = dict(directions_df[["RHEA_ID_LR", "RHEA_ID_MASTER"]].values)
        direction_rl = dict(directions_df[["RHEA_ID_RL", "RHEA_ID_MASTER"]].values)
        direction_bi = dict(directions_df[["RHEA_ID_BI", "RHEA_ID_MASTER"]].values)

        # Columns are: RHEA_ID	DIRECTION	MASTER_ID	ID
        eccode_df = MODULE.ensure_csv(
            version,
            url="https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv",
            read_csv_kwargs=dict(dtype=str),
        )
        master_to_eccode = dict(eccode_df[["MASTER_ID", "ID"]].values)

        # Search for all rxn files
        paths = sorted(directory.joinpath("rxn").glob("*.rxn"))

        rows = []
        for path in tqdm(
            paths, unit_scale=True, unit="reaction", desc="reading reactions"
        ):
            rhea_id = path.stem
            if rhea_id in direction_master:
                direction = "master"
                master_id = rhea_id
            elif rhea_id in direction_lr:
                direction = "lr"
                master_id = direction_lr[rhea_id]
            elif rhea_id in direction_rl:
                direction = "rl"
                master_id = direction_rl[rhea_id]
            elif rhea_id in direction_bi:
                direction = "bi"
                master_id = direction_bi[rhea_id]
            else:
                tqdm.write(f"failed to find the direction for {path.name}")
                continue

            ec_code = master_to_eccode.get(master_id)
            ec_class = UNKNOWN_LABEL if ec_code is None else ec_code.split(".")[0]

            try:
                rxn = AllChem.ReactionFromRxnFile(path.as_posix())
                rxn.Initialize()
                # rdChemReactions.PreprocessReaction(rxn)
                rows.append(
                    (
                        rhea_id,
                        direction,
                        ec_code,
                        ec_class,
                        rdChemReactions.ReactionToSmiles(rxn),
                    )
                )
            except (ValueError, RuntimeError):
                tqdm.write(f"failed to parse {path.name}")

        metadata_df = pd.DataFrame(
            rows,
            columns=["rhea_id", "direction", "enzyme_code", "enzyme_class", "smiles"],
        ).set_index("rhea_id")
        metadata_df.to_csv(metadata_path, sep="\t")

        click.echo("Calculating differential reaction fingerprints")
        fingerprint_df = pd.DataFrame(
            DrfpEncoder.encode(metadata_df.smiles, show_progress_bar=progress),
            index=metadata_df.index,
        )
        click.echo(f"Saving differential reaction fingerprints to {fingerprints_path}")
        fingerprint_df.to_csv(fingerprints_path, sep="\t")

    click.echo("Rewriting readme")
    README_PATH.write_text(index_template.render(version=version) + "\n")

    fig, (lax, rax) = plt.subplots(1, 2, figsize=(10, 4))

    click.echo("Calculating scree plot")
    pca_full = PCA(random_state=random_state)
    pca_full.fit(fingerprint_df)
    y = np.cumsum(pca_full.explained_variance_ratio_)
    x = np.arange(y.shape[0])
    sns.lineplot(x=x, y=y, ax=lax)
    lax.axhline(0.80, linestyle="--", color="red")
    lax.axhline(0.90, linestyle="--", color="goldenrod")
    lax.axhline(0.95, linestyle="--", color="green")
    lax.set_title(f"PCA Scree Plot (Rhea v{version})")
    lax.set_xlabel("Number Components")
    lax.set_ylabel("Cumulative Explained Variance")

    reducer_2d = Isomap(n_neighbors=2)
    click.echo(f"Transforming to 2D with {reducer_2d.__class__.__name__}")
    transformed_df = pd.DataFrame(
        reducer_2d.fit_transform(fingerprint_df),
        index=fingerprint_df.index.map(str),
        columns=["PC1", "PC2"],
    )

    click.echo("Applying K-means")
    kmeans = KMeans(7, random_state=random_state)
    transformed_df["cluster"] = kmeans.fit_predict(
        PCA(64, random_state=random_state).fit_transform(fingerprint_df)
    ).astype(str)
    transformed_df["class"] = list(metadata_df["enzyme_class"])
    transformed_df["class_label"] = transformed_df["class"].map(LABELS)
    transformed_df.to_csv(fingerprints_2d_path, sep="\t")
    g = sns.scatterplot(
        data=transformed_df[transformed_df["class"] != UNKNOWN_LABEL],
        x="PC1",
        y="PC2",
        hue="class_label",
        alpha=0.2,
        ax=rax,
    )
    g.legend_.set_title(None)
    rax.set_title(f"PCA 2D Reduction (Rhea v{version})")

    plt.tight_layout()
    fig.savefig(output.joinpath("scatter.svg"))
    fig.savefig(output.joinpath("scatter.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
