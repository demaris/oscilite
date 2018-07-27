from setuptools import setup

setup(name='cml_mini',
    version='0.0.3',
    description='diffusive and competitive coupled map lattice, visualization, and analysis utilities',
    url='http://github.com/demaris/cml_mini',
    author='demaris',
    author_email='davidl.demaris@gmail.com',
    license='MIT',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 2.7',
    ],	
    packages=['cml_mini'],
    data_files=[("",["LICENSE"])],
    zip_safe=False)
