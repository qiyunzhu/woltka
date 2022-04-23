# Installation

## Operating system

There is no restriction as far as we are aware of. Tested and working on Linux, macOS and Windows systems.

## Software environment

Woltka is written in **Python 3**. One needs at least Python 3.6 to run the program.

We recommend [Conda](https://docs.conda.io/en/latest/) for managing Python version and packages. One can create a conda environment and install necessary dependencies (Cython and BIOM format) by:

```bash
conda create -n woltka python=3
conda activate woltka
```

If you already have a [QIIME 2](https://qiime2.org/) environment, these steps can be omitted as the dependencies are already included. See [details](../woltka/q2).

## Installation

Option 1: Install the latest release using [pip](https://pypi.org/project/pip/):

```bash
pip install woltka
```

Option 2: Install the current development from GitHub:

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
```

Option 3: Install from a local copy:

Download this [repository](https://github.com/qiyunzhu/woltka/archive/master.zip) or any of the previous [releases](https://github.com/qiyunzhu/woltka/releases). Unzip and navigate to the package directory. Then execute:

```bash
python setup.py install
```

Type `woltka` to check if installation is successful, in which case command-line help information will be displayed on the screen.

## Acceleration

Woltka has a [Numba](https://numba.pydata.org/)-accelerated version (the [numba](https://github.com/qiyunzhu/woltka/tree/numba) branch) in parallel to the master branch. In this version, the coordinates-based read-gene matching (see [details](ordinal,md)) is significantly faster and consumes less memory. To use this feature, install Woltka using the following way instead:

```bash
conda install -c conda-forge numba biom-format
pip install git+https://github.com/qiyunzhu/woltka.git@numba
```

[**Note**] This feature only accelerates read-gene matching. It does not help with the analysis of microbiome structure.

## Upgrade

Just add `--upgrade` or `-U` to the pip command:

```bash
pip install -U woltka
```

## Uninstallation

```bash
pip uninstall woltka
```

If you no longer need the conda environment:

```bash
conda env remove -n woltka --all
```

## Compatibility

If in the future some dependencies have changes that are not compatible with the current release of Woltka, the following "safe" commands can be used to install the current versions of dependencies.

```bash
conda create -n woltka python=3.8.2
conda activate woltka
conda install -c conda-forge cython=0.29.6 biom-format=2.1.8
```

## Test

You may test whether Woltka functions correctly by running it on some small test datasets. You will need to download the [repository](https://github.com/qiyunzhu/woltka/archive/master.zip) to get them. See [example usage](../README.md#example-usage) for how to run these small tests manually. Alternatively, you may execute the following command in the repository directory:

```bash
python -m unittest
```

This will run all tests automatically. If it prints "OK" in the end, everything is okay. Otherwise it will report specific error(s), and in such case, please contact the development team for solutions.
