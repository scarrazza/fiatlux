import sys, os

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

# load long description from README
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="fiatlux",
    version="0.1.1",
    description="API for LUXqed methodology in global PDF fits",
    author="Stefano Carrazza",
    license="AGPLv3",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/fiatlux",
    cmake_args=['-DPYTHON_ONLY:BOOL=ON'],
    extras_require={"test": ["pytest"]},
    python_requires=">=3.6",
    long_description=long_description,
    long_description_content_type="text/markdown",
)
