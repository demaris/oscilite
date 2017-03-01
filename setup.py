from setuptools import setup

setup(name='cml_basic',
    version='0.1.2',
    description='basic diffusive coupled map lattice and analysis utilities',
    url='http://github.com/demaris/cml_basic',
    author='demaris',
    author_email='davidl.demaris@gmail.com',
    license='MIT',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 2.7',
    ],	
    packages=['cml_basic'],
    data_files=[("",["LICENSE.txt"])],
    zip_safe=False)
