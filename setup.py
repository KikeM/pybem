from setuptools import setup, find_packages

setup(
    name             = 'pybem',
    author_email     = 'enrique.millanvalbuena@gmail.com',
    packages         = find_packages(),
    install_requires = [
        'numpy',
        'scipy',
        'pandas',
        'fluids'
    ]

)