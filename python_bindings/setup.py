from setuptools import setup, Extension
from Cython.Build import cythonize

import os
os.environ["CC"] = "g++"

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mass_dawg",
    version="1.1.1", 
    url="https://github.com/zmcgrath96/mass_DAWG",
    author="Zachary McGrath", 
    author_email="zmcgrath96@gmail.com", 
    description="A datastructure for efficient storage of mass spectrometry data", 
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="bioinformatics mass_spectrometry ms dawg graph prefix",
    ext_modules=cythonize(Extension(
        "mass_dawg", 
        ["mass_dawg.pyx"], 
        language="c++", 
        extra_compile_args=["-std=c++11"]
    ))
)