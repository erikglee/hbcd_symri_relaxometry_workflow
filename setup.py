from setuptools import setup, find_packages

setup(
    name="HBCD_SYMRI_REALXOMETRY_WORKFLOW",
    version="0.1.0",
    description="A brief description of the package",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/erikglee/hbcd_symri_relaxometry_workflow",
    packages=find_packages(),
    install_requires=[
        # List your package dependencies here
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)