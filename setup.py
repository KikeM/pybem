import versioneer
from setuptools import setup, find_packages

setup(
    name="pybem",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Robust solution to propeller equations.",
    author_email="enrique.millanvalbuena@gmail.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=["numpy", "scipy", "fluids"],
)
