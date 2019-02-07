import os
import sys

from setuptools import setup


# This is a temporary hack for installing igraph on Mac with Anaconda
# If this stays in long term, the rest of the setup should be wrapped in a context manager.
# That way, if setup fails, we can return state on exit.
IS_MAC = sys.platform == 'darwin'
MAYBE_CONDA = any('CONDA' in key for key in os.environ)
SET_TARGET = IS_MAC and MAYBE_CONDA
if SET_TARGET:
    current_target = os.environ.get('MACOSX_DEPLOYMENT_TARGET', '')
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.14'

HERE = os.path.abspath(os.path.dirname(__file__))
about = {}
version = os.path.join(HERE, 'seekr', '__version__.py')
with open(version, 'r', encoding='utf-8') as f:
    exec(f.read(), about)

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = [
    'cython',
    'tqdm',
    'numpy',
    'pandas',
    'requests',
    'networkx',
    'python-igraph',
    'louvain',
    'leidenalg',
    'ushuffle'
]

test_requirements = [
    'pytest'
]

setup(name=about['__title__'],
      version=about['__version__'],

      install_requires=requirements,
      tests_require=test_requirements,
      description=about['__description__'],
      long_description=long_description,
      long_description_content_type="text/markdown",
      url=about['__url__'],
      author=about['__author__'],
      author_email=about['__author_email__'],
      license=about['__license__'],
      packages=['seekr'],
      zip_safe=False,
      classifiers=['Intended Audience :: Science/Research',
                   'License :: OSI Approved :: MIT License',
                   'Natural Language :: English',
                   'Programming Language :: Python :: 3.6'],
      entry_points = {'console_scripts':
          ['seekr_download_gencode = seekr.console_scripts:console_download_gencode',
           'seekr_canonical_gencode = seekr.console_scripts:console_canonical_gencode',
           'seekr_kmer_counts = seekr.console_scripts:console_kmer_counts',
           'seekr_pearson = seekr.console_scripts:console_pearson',
           'seekr_norm_vectors = seekr.console_scripts:console_norm_vectors',
           'seekr_rand_rnas = seekr.console_scripts:console_gen_rand_rnas',
           'seekr_graph = seekr.console_scripts:console_graph',
           'seekr = seekr.console_scripts:console_seekr_help']})

if SET_TARGET:
    if current_target:
        os.environ['MACOSX_DEPLOYMENT_TARGET'] = current_target
    else:
        del os.environ['MACOSX_DEPLOYMENT_TARGET']
