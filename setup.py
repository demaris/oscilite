from setuptools import setup

setup(name='cml_basic'
    version='0.1',
    description='basic diffusive coupled map lattice and analysis utilities',
    url='http://github.com/demaris/cml_basic',
    author='demaris',
    author_email='davidl.demaris@gmail.com',
    license='MIT',
    packages=['cml_basic'],
    data_files=[("",["LICENSE.txt"])],
    install_requires=[
        numpy,
        scipy.stats,
        scipy.signal,
        matplotlib,
    ],
    zip_safe=False)