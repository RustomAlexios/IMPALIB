# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from skbuild import setup

setup(
    name="pyimpa",
    version="1.0.0",
    description="A package for MPA-based solvers to linear programs",
    author='BP-OPT Team',
    license="MIT",
    package_dir={"": "src", "": "test"},
    packages=['impa', 'python_tsp'],
    python_requires=">=3.9",
    install_requires=["bitstring", "tqdm", "numpy", "elkai", "python_tsp"],
    scripts=["src/impa/main_kc_mwm.py"]
)