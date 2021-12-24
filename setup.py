#!/usr/bin/env python

from distutils.core import setup
from mutmap.__init__ import __version__

setup(name='mutmap',
      version='{}'.format(__version__),
      description='MutMap: pipeline to identify causative mutations responsible for a phenotype',
      author='Yu Sugihara',
      author_email='sugihara.yu.85s@kyoto-u.jp',
      url='https://github.com/YuSugihara/MutMap',
      license='GPL',
      packages=['mutmap'],
      entry_points={'console_scripts': [
            'mutmap = mutmap.mutmap:main',
            'mutplot = mutmap.mutplot:main',
            ]
        }
    )
