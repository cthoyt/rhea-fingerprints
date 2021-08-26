import pathlib

import bioversions
import pandas as pd
import pystow
from drfp import DrfpEncoder
from rdkit.Chem import AllChem, rdChemReactions
from tqdm import tqdm

HERE = pathlib.Path(__file__).parent.resolve()
OUTPUT = HERE.joinpath("output")
OUTPUT.mkdir(exist_ok=True, parents=True)


def main():
    """Generate SMILES and differential reaction fingerprints for each Rhea entry."""
    version = bioversions.get_version("rhea")
    module = pystow.module("bio", "rhea", version)
    output_directory = OUTPUT.joinpath(version)
    output_directory.mkdir(exist_ok=True, parents=True)
    smiles_path = output_directory.joinpath("reaction_smiles.tsv")
    fingerprints_path = output_directory.joinpath("reaction_fingerprints.pkl.gz")

    directory = module.ensure_untar(
        url="ftp://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz",
    )

    paths = list(directory.joinpath("rxn").glob("*.rxn"))

    rows = []
    for path in tqdm(paths, unit_scale=True, unit="reaction"):
        try:
            rxn = AllChem.ReactionFromRxnFile(path.as_posix())
            rxn.Initialize()
            # rdChemReactions.PreprocessReaction(rxn)
            rows.append((
                path.stem,
                rdChemReactions.ReactionToSmiles(rxn),
            ))
        except (ValueError, RuntimeError) as e:
            tqdm.write(f"failed: {path.name}")

    rxn_df = pd.DataFrame(rows, columns=["rhea_id", "smiles"])
    rxn_df.to_csv(smiles_path, sep="\t", index=False)

    # Add differential reaction fingerprints
    rxn_df["drfp"] = DrfpEncoder.encode(rxn_df.smiles)

    # Save full dataframe as a pickle for reuse later
    rxn_df.to_pickle(fingerprints_path)


if __name__ == '__main__':
    main()
