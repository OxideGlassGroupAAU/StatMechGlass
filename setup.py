import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='StatMechGlass',
    version='0.1.5',
    author='Mikkel Bødker',
    author_email='mikkelboedker@gmail.com',
    description='Calculates structure of oxide glasses using statistical mechanics',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/OxideGlassGroupAAU/StatMechGlass",
    install_requires=[
        'numpy>=1.16',
        'scipy>=1.2',
        'matplotlib>=3.4.2',
        'sklearn',
    ],
    keywords='glass, statistical mechanics',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires='>=3.6',
    include_package_data=True,
)
