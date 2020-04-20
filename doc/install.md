# Installation

## Operating system

There is no restriction as far as we are aware of. Tested and working on Linux, macOS and Windows systems.


## Software environment

Woltka is written in Python 3. One needs at least Python 3.6 to run the program.

We recommend [Conda](https://docs.conda.io/en/latest/) for managing Python version and packages. One can create a conda environment and install necessary dependencies by:

```bash
conda create -n woltka python=3
conda activate woltka
conda install -c conda-forge cython biom-format
```

If you already have a [QIIME 2](https://qiime2.org/) environment, these steps can be omitted as the dependencies are already included. See [details](../woltka/q2).

## Installation

Option 1: Install from GitHub:

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
```

Option 2: Install from a local copy:

Download this [repository](https://github.com/qiyunzhu/woltka/archive/master.zip). Unzip. Then execute:

```bash
python setup.py install
```

Type `woltka` to check if installation is successful, in which case command-line help information will be displayed on the screen.

## Upgrade

Just add `--upgrade` or `-U` to the pip command:

```bash
pip install -U git+https://github.com/qiyunzhu/woltka.git
```

## Compatibility

If in the future some dependencies have changes that are not compatible with the current release of Woltka, the following "safe" commands can be used to install the current versions of dependencies.

```bash
conda create -n woltka python=3.8.2
conda activate woltka
conda install -c conda-forge cython=0.29.6 biom-format=2.1.8
```
