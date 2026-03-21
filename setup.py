#!/usr/bin/env python

from setuptools import find_packages, setup
from mutmap.__init__ import __version__

setup(name='mutmap',
      version='{}'.format(__version__),
      description='MutMap: pipeline to identify causative mutations responsible for a phenotype',
      author='Yu Sugihara',
      author_email='sugihara.yu.85s@kyoto-u.jp',
      url='https://github.com/YuSugihara/MutMap',
      license='GPL',
      packages=find_packages(exclude=['test', 'docs']),
      python_requires='>=3.5',
      install_requires=[
            'matplotlib',
            'numpy',
            'pandas',
      ],
      extras_require={
            'plot-style': ['seaborn'],
      },
      entry_points={'console_scripts': [
            'mutmap = mutmap.mutmap:main',
            'mutplot = mutmap.mutplot:main',
            ]
        }
    )
