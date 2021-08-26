import pathlib

import bioversions
import pandas as pd
import pystow
from drfp import DrfpEncoder
from rdkit.Chem import AllChem, rdChemReactions
from tqdm import tqdm

HERE = pathlib.Path(__file__).parent.resolve()


def main():
    """Generate SMILES and differential reaction fingerprints for each Rhea entry."""
    # Get the current version of Rhea using :mod:`bioversions`.
    version = bioversions.get_version("rhea")

    # Create a version-specific output directory
    output = HERE.joinpath("output", version)
    output.mkdir(exist_ok=True, parents=True)

    # Use a version-specific PyStow directory for reproducible downloading
    # of external resources
    module = pystow.module("bio", "rhea", version)

    tqdm.write("Downloading/unpacking reaction archive")
    directory = module.ensure_untar(
        url="ftp://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz",
    )
    tqdm.write("Downloaded/unpacking reaction archive")

    # Search for all rxn files
    paths = list(directory.joinpath("rxn").glob("*.rxn"))

    rows = []
    for path in tqdm(paths, unit_scale=True, unit="reaction"):
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
        except (ValueError, RuntimeError) as e:
            tqdm.write(f"failed: {path.name}")

    rxn_df = pd.DataFrame(rows, columns=["rhea_id", "smiles"])
    rxn_df.to_csv(output.joinpath("reaction_smiles.tsv"), sep="\t", index=False)

    # Add differential reaction fingerprints
    rxn_df["drfp"] = DrfpEncoder.encode(rxn_df.smiles)

    # Save full dataframe as a pickle for reuse later
    rxn_df.to_pickle(output.joinpath("reaction_fingerprints.pkl.gz"))


if __name__ == "__main__":
    main()
