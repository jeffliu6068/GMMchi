from setuptools import setup

setup(
    name='GMMchi',
    version='0.1',
    author='Jeff Liu',
    author_email='jeffliu6068@gmail.com',
    packages=['GMMchi'],
    url='http://pypi.python.org/pypi/GMMchi/',
    license='LICENSE.txt',
    description='GMM with chi-square protocol',
    long_description=open('README.txt').read(),
    install_requires=[
        'pandas',
        'scipy',
        'numpy',
        'matplotlib',
        'seaborn',
        'sklearn',
        'tqdm',
    ],
)
