import pathlib

import bioversions
import click
import matplotlib.pyplot as plt
import pandas as pd
import pystow
import seaborn as sns
from drfp import DrfpEncoder
from more_click import force_option, verbose_option
from rdkit.Chem import AllChem, rdChemReactions
from sklearn.decomposition import PCA
from tqdm import tqdm

HERE = pathlib.Path(__file__).parent.resolve()
MODULE = pystow.module("bio", "biogrid")


@click.command()
@verbose_option
@force_option
def main(force: bool):
    """Generate SMILES and differential reaction fingerprints for each Rhea entry."""
    # Get the current version of Rhea using :mod:`bioversions`.
    version = bioversions.get_version("rhea")

    # Create a version-specific output directory
    output = HERE.joinpath("output", version)
    output.mkdir(exist_ok=True, parents=True)
    reaction_fingerprints_path = output.joinpath("reaction_fingerprints.tsv.gz")

    if reaction_fingerprints_path.is_file() and not force:
        click.echo("Loading reaction dataframe")
        fingerprint_df = pd.read_csv(reaction_fingerprints_path, sep="\t")
    else:
        # Use a version-specific PyStow directory for reproducible downloading
        # of external resources
        tqdm.write("Downloading/unpacking reaction archive")
        directory = MODULE.ensure_untar(
            version,
            url="ftp://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz",
        )
        tqdm.write("Downloaded/unpacking reaction archive")

        # Search for all rxn files
        paths = list(directory.joinpath("rxn").glob("*.rxn"))

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
        fingerprint_df.to_csv(reaction_fingerprints_path, sep="\t", index=False)

    pca = PCA(2)
    transformed = pd.DataFrame(
        pca.fit_transform(fingerprint_df), index=fingerprint_df.index
    )
    fig, ax = plt.subplots()
    sns.scatterplot(
        data=transformed, x=transformed.columns[0], y=transformed.columns[1], ax=ax
    )
    fig.savefig(output.joinpath("scatter.svg"))
    fig.savefig(output.joinpath("scatter.png"), dpi=300)


if __name__ == "__main__":
    main()
