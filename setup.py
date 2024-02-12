import setuptools

setuptools.setup(
    name="xyz2atoms",
    version="1.0.0",
    author="Mykola Melnychenko",
    url="https://github.com/MelnychenkoM/xyz2atoms",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    install_requires = ['numpy', 'plotly', 'pandas'],
    python_requires='>=3.5',
)