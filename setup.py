"""
Peptide sequence foldability prediciton tool
"""
from setuptools import setup

dependencies = ['numpy', 'pandas', 'matplotlib', 'colorama']

setup(
    name='prefold',
    version='0.1.0',
    url='https://github.com/aretas2/preFold',
    license='BSD',
    author='Aretas Gaspariunas',
    author_email='aretasgasp@gmail.com',
    description='Peptide sequence foldability prediciton tool',
    long_description=open('README.md').read(),
    packages=['prefold_py'],
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    data_files = [('', ['prefold_py/datasets/HB_DATA_NORM.DAT', 'prefold_py/datasets/PKA_DATA_VOET.DAT'])],
    entry_points={
        'console_scripts': ['prefold = prefold_py.prefold_cli:main'],
    },
    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Developers :: Scientists',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Software Development :: Libraries :: Python Modules'
        ]
)
