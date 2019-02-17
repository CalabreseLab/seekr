import os
import sys

from setuptools import setup


# This is a temporary hack for installing igraph on Mac with Anaconda
# If this stays in long term, the rest of the setup should be wrapped in a context manager.
# That way, if setup fails, we can return state on exit.
IS_MAC = sys.platform == 'darwin'
MAYBE_CONDA = any('CONDA' in key for key in os.environ)
NO_TARGET = 'MACOSX_DEPLOYMENT_TARGET' not in os.environ
SET_TARGET = IS_MAC and MAYBE_CONDA and NO_TARGET
if SET_TARGET:
    msg = ('It looks like you are using OSX and Anaconda. '
           'Please run the command:\n\n'
           '$ export MACOSX_DEPLOYMENT_TARGET=10.14\n\n'
           'Then, rerun:\n\n'
           '$ pip install seekr\n\n'
           'For more details see: \n\n'
           'https://github.com/igraph/python-igraph/issues/208\n\n'
           'This issue will hopefully be fixed soon.')
    raise OSError(msg)

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
    'ushuffle',
    'matplotlib',
    'seaborn'
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
           'seekr_visualize_distro = seekr.console_scripts:console_visualize_distro',
           'seekr_norm_vectors = seekr.console_scripts:console_norm_vectors',
           'seekr_rand_rnas = seekr.console_scripts:console_gen_rand_rnas',
           'seekr_graph = seekr.console_scripts:console_graph',
           'seekr_domain_pearson = seekr.console_scripts:console_domain_pearson',
           'seekr = seekr.console_scripts:console_seekr_help']})
