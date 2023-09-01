from skbuild import setup


setup(
    name="pyimpa",
    version="1.0.0",
    description="A package for MPA-based solvers to linear programs",
    author='BP-OPT Team',
    license="MIT",
    package_dir={"": "src"},
    packages=['impa'],
    python_requires=">=3.7",
    install_requires=["bitstring", "tqdm"],
    scripts=["src/impa/main_wrapper_optimized.py"]
)
