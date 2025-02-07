from setuptools import setup, find_packages

from riboss import __version__


setup(
    name='riboss',
    version=__version__,
    description='Comparing translational potential of open reading frames within individual transcripts',

    url='https://github.com/lcscs12345/riboss',
    author='CS Lim',
    author_email='lcscs12345@gmail.com',

    packages=find_packages(exclude=['tests', 'tests.*']),

    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
)
