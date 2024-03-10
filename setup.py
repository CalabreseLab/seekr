import os

from setuptools import setup


HERE = os.path.abspath(os.path.dirname(__file__))
about = {}
version = os.path.join(HERE, "seekr", "__version__.py")
with open(version, "r", encoding="utf-8") as f:
    exec(f.read(), about)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    "cython",
    "tqdm",
    "numpy",
    "pandas",
    "requests",
    "networkx>=3.1",
    "python-igraph==0.10.4",
    "leidenalg==0.9.1",
    "matplotlib>=3.7",
    "seaborn",
    "scipy",
    "statsmodels",
]

# ensure matplotlib to be greater than 3.5.3 which is compatible with python 3.7 and above
# and could accomodate the fontmanager function

test_requirements = ["pytest"]

setup(
    name=about["__title__"],
    python_requires=">=3.9",
    version=about["__version__"],
    install_requires=requirements,
    tests_require=test_requirements,
    description=about["__description__"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=about["__url__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    license=about["__license__"],
    packages=["seekr"],
    package_data={
        'seekr': ['data/*'],  # Include all files in the 'data' sub-directory
    },
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.9",
    ],
    entry_points={
        "console_scripts": [
            "seekr_download_gencode = seekr.console_scripts:console_download_gencode",
            "seekr_filter_gencode = seekr.console_scripts:console_filter_gencode",
            "seekr_kmer_counts = seekr.console_scripts:console_kmer_counts",
            "seekr_pearson = seekr.console_scripts:console_pearson",
            "seekr_norm_vectors = seekr.console_scripts:console_norm_vectors",
            "seekr = seekr.console_scripts:console_seekr_help",
            "seekr_find_dist = seekr.console_scripts:console_find_dist",
            "seekr_find_pval = seekr.console_scripts:console_find_pval",
            "seekr_adj_pval = seekr.console_scripts:console_adj_pval",
            "seekr_kmer_heatmap = seekr.console_scripts:console_kmer_heatmap",
            "seekr_kmer_dendrogram = seekr.console_scripts:console_kmer_dendrogram",
            "seekr_kmer_leiden = seekr.console_scripts:console_kmer_leiden",
            "seekr_kmer_count_barplot = seekr.console_scripts:console_kmer_count_barplot",
            "seekr_kmer_msd_barplot = seekr.console_scripts:console_kmer_msd_barplot",
            "seekr_kmer_comp_textplot = seekr.console_scripts:console_kmer_comp_textplot",
            "seekr_kmer_indi_textplot = seekr.console_scripts:console_kmer_indi_textplot",
        ]
    },
)
