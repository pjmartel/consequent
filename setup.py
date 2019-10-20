from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='consequent',
      version='0.1',
      description='Biological sequence tools',
      long_description=readme(),
      url='https://github.com/pjmartel/consequent',
      author='Paulo Martel',
      author_email='pmartel@ualg.pt',
      license='MIT',
      packages=['consequent'],
      scripts=[
          'bin/align_check.py',
          'bin/align_check_mp.py',
          'bin/dot_plot.py'
      ],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'biopython',
      ],
      zip_safe=False)
