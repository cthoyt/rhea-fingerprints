import pathlib

import bioversions
import click
import pystow

__all__ = [
    "HERE",
    "TEMPLATES",
    "OUTPUT",
    "MODULE",
    "README_PATH",
    "LABELS",
    "UNKNOWN_LABEL",
    "version_option",
    "random_state_option",
]

HERE = pathlib.Path(__file__).parent.resolve()
TEMPLATES = HERE.joinpath("templates")
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

version_option = click.option(
    "--version",
    default=lambda: bioversions.get_version("rhea"),
    help="The version of Rhea. If none given, gets the current version with Bioversions",
)
random_state_option = click.option("--random-state", type=int, default=0)
