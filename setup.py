from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='seekr',
      version='1.1.3',
      install_requires=['tqdm', 'numpy', 'pandas'],
      description='A library for counting small kmer frequencies in nucleotide sequences.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/CalabreseLab/seekr',
      author='Jessime Kirk',
      author_email='jessime.kirk@gmail.com',
      license='MIT',
      packages=['seekr'],
      zip_safe=False,
      entry_points = {'console_scripts':
          ['kmer_counts = seekr.console_scripts:console_kmer_counts',
           'pearson = seekr.console_scripts:console_pearson',
           'norm_vectors = seekr.console_scripts:console_norm_vectors']})
