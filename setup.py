import sys

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

setup(
    name="fiatlux",
    version="0.1.0",
    description="API for LUXqed methodology in global PDF fits",
    author="Stefano Carrazza",
    license="AGPLv3",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/fiatlux",
    cmake_args=['-DPYTHON_ONLY:BOOL=ON'],
    extras_require={"test": ["pytest"]},
    python_requires=">=3.6",
)
