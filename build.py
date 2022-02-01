import pathlib

import bioversions
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pystow
import seaborn as sns
from drfp import DrfpEncoder
from jinja2 import Environment, FileSystemLoader
from more_click import force_option, verbose_option
from rdkit.Chem import AllChem, rdChemReactions
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from tqdm import tqdm

HERE = pathlib.Path(__file__).parent.resolve()
MODULE = pystow.module("bio", "rhea")

README_PATH = HERE.joinpath("README.md")

environment = Environment(
    autoescape=True, loader=FileSystemLoader(HERE), trim_blocks=False
)
index_template = environment.get_template("README.md.jinja")


@click.command()
@verbose_option
@force_option
@click.option("--random-state", type=int, default=0)
def main(force: bool, random_state: int):
    """Generate SMILES and differential reaction fingerprints for each Rhea entry."""
    # Get the current version of Rhea using :mod:`bioversions`.
    version = bioversions.get_version("rhea")

    # Create a version-specific output directory
    output = HERE.joinpath("output", version)
    output.mkdir(exist_ok=True, parents=True)
    reaction_fingerprints_path = output.joinpath("reaction_fingerprints.tsv.gz")

    if reaction_fingerprints_path.is_file() and not force:
        click.echo("Loading reaction dataframe")
        fingerprint_df = pd.read_csv(reaction_fingerprints_path, sep="\t", index_col=0)
        click.echo("Loaded reaction dataframe")
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

        # Search for all rxn files
        paths = sorted(directory.joinpath("rxn").glob("*.rxn"))

        rows = []
        for path in tqdm(
            paths, unit_scale=True, unit="reaction", desc="reading reactions"
        ):
            try:
                rxn = AllChem.ReactionFromRxnFile(path.as_posix())
                rxn.Initialize()
                # rdChemReactions.PreprocessReaction(rxn)
                rows.append(
                    (
                        path.stem,
                        rdChemReactions.ReactionToSmiles(rxn),
                    )
                )
            except (ValueError, RuntimeError):
                tqdm.write(f"failed: {path.name}")

        smiles_df = pd.DataFrame(rows, columns=["rhea_id", "smiles"])
        smiles_df.to_csv(output.joinpath("reaction_smiles.tsv"), sep="\t", index=False)

        # Add differential reaction fingerprints
        fingerprint_df = pd.DataFrame(
            DrfpEncoder.encode(smiles_df.smiles), index=smiles_df.rhea_id
        )
        fingerprint_df.to_csv(reaction_fingerprints_path, sep="\t")

    README_PATH.write_text(index_template.render(version=version) + "\n")

    fig, (lax, rax) = plt.subplots(1, 2, figsize=(10, 4))

    # Calculating scree plot
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

    pca_2d = PCA(2, random_state=random_state)
    transformed_df = pd.DataFrame(
        pca_2d.fit_transform(fingerprint_df),
        index=fingerprint_df.index.map(str),
        columns=["PC1", "PC2"],
    )

    kmeans = KMeans(7, random_state=random_state)
    transformed_df["cluster"] = kmeans.fit_predict(
        PCA(64, random_state=random_state).fit_transform(fingerprint_df)
    ).astype(str)
    transformed_df.to_csv(output.joinpath("reaction_fingerprints_2d.tsv"), sep="\t")
    sns.scatterplot(
        data=transformed_df,
        x="PC1",
        y="PC2",
        hue="cluster",
        ax=rax,
    )
    rax.set_title(f"PCA 2D Reduction (Rhea v{version})")

    plt.tight_layout()
    fig.savefig(output.joinpath("scatter.svg"))
    fig.savefig(output.joinpath("scatter.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
