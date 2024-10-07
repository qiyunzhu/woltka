# Installation

## Operating system

There is no restriction as far as we are aware of. Tested and working on **Linux**, **macOS** and **Windows** systems.

## Software environment

Woltka is written in **Python 3**. One needs at least Python 3.6 to run the program.

We recommend [Conda](https://docs.conda.io/en/latest/) for managing Python version and packages. One can create a conda environment for Woltka by:

```bash
conda create -n woltka python=3
conda activate woltka
```

If you already have a [QIIME 2](https://qiime2.org/) environment, these steps can be omitted as the dependencies are already included. See [details](../woltka/q2).

You may also use [Mamba](https://mamba.readthedocs.io/en/latest/) instead of Conda. The procedures discussed in this page that apply to Conda also apply to Mamba.

## Installation

Option 1: Install the latest release from the [Bioconda](https://bioconda.github.io/) channel:

```bash
conda install -c conda-forge -c bioconda woltka
```

- **Note**: It is necessary to also specify the [conda-forge](https://conda-forge.org/) channel, in order to install the latest version of [biom-format](https://biom-format.org/), a dependency of Woltka.

Option 2: Install the latest release using [pip](https://pypi.org/project/pip/):

```bash
pip install woltka
```

Option 3: Install the current development from GitHub:

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
```

Option 4: Install from a local copy:

Download this [repository](https://github.com/qiyunzhu/woltka/archive/main.zip) or any of the previous [releases](https://github.com/qiyunzhu/woltka/releases). Unzip and navigate to the package directory. Then execute:

```bash
python setup.py install
```

Type `woltka` to check if installation is successful, in which case command-line help information will be displayed on the screen.

## Upgrade

If you installed Woltka using Conda, do:

```bash
conda update -c conda-forge -c bioconda woltka
```

If you installed Woltka using pip, do:

```bash
pip install -U woltka
```

## Uninstallation

If you installed Woltka using Conda, do:

```bash
conda remove woltka
```

If you installed Woltka using pip, do:

```bash
pip uninstall woltka
```

If you no longer need the Conda environment:

```bash
conda env remove -n woltka
```

## Compatibility

If in the future some dependencies have changes that are not compatible with the current release of Woltka, the following "safe" commands can be used to install the current versions of dependencies.

```bash
conda create -n woltka python=3.12.7
conda activate woltka
conda install -c conda-forge numba=0.60.0 biom-format=2.1.16
conda install -c bioconda woltka=0.1.7
```

## Test

You may test whether Woltka functions correctly by running it on some small test datasets. You will need to download the [repository](https://github.com/qiyunzhu/woltka/archive/master.zip) to get them. See [example usage](../README.md#example-usage) for how to run these small tests manually. Alternatively, you may execute the following command in the repository directory:

```bash
python -m unittest
```

This will run all tests automatically. If it prints "OK" in the end, everything is okay. Otherwise it will report specific error(s), and in such case, please contact the development team for solutions.
