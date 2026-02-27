#!/usr/bin/env python
"""Setup configuration for template_cleaner package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="template_cleaner",
    version="1.0.0",
    author="Shamik Ghosh",
    author_email="thequarkexpress@gmail.com",
    description="Template-based E-to-B leakage cleaner for CMB B-mode analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/1cosmologist/template_cleaner",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "healpy",
    ],
)
