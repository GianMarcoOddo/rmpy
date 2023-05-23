from setuptools import setup, find_packages
import pathlib

# Reading the content of the README.md file
here = pathlib.Path(__file__).parent
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="rmpy",
    version="1.2.4",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'yfinance',
        'scipy'],
    author="Gian Marco Oddo, Mattia Aprea, Sofia Compagnoni",
    author_email="gian.marco.oddo@usi.ch, mattia.aprea@usi.ch, sofia.compagnoni@usi.ch",
    description="This Python package provides a comprehensive set of tools for calculating and analyzing Non-Parametric and Parametric Value-at-Risk (NpVaR and pVaR).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/GianMarcoOddo/rmpy",
    license="MIT",
)
