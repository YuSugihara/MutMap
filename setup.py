#!/usr/bin/env python

from distutils.core import setup

setup(name='mutmap',
      version='2.0.4',
      description='MutMap: pipeline to identify causative mutations responsible for a phenotype',
      author='Yu Sugihara',
      author_email='yu57th@gmail.com',
      url='https://github.com/YuSugihara/MutMap',
      license='GPL',
      packages=['mutmap'],
      entry_points={'console_scripts': [
            'mutmap = mutmap.mutmap:main',
            'mutplot = mutmap.mutplot:main',
            ]
        }
    )
