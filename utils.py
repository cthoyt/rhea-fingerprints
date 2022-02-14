import pathlib

import pystow

__all__ = [
    "HERE",
    "OUTPUT",
    "MODULE",
    "README_PATH",
    "LABELS",
    "UNKNOWN_LABEL",
]

HERE = pathlib.Path(__file__).parent.resolve()
OUTPUT = HERE.joinpath("docs")
MODULE = pystow.module("bio", "rhea")

README_PATH = HERE.joinpath("README.md")

UNKNOWN_LABEL = "Unknown"
LABELS = {
    "1": "Oxidoreductases",
    "2": "Transferases",
    "3": "Hydrolases",
    "4": "Lyases",
    "5": "Isomerases",
    "6": "Ligases",
    "7": "Translocases",
}
