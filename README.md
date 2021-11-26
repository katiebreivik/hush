<p align="center">
<a href="https://github.com/rodluger/showyourwork">
</a>
<br>
<br>
<a href="mailto:sarahgthiele@gmail.com?cc=kbreivik@flatironinstitute.org">
      <img src="https://img.shields.io/badge/contact-authors-blueviolet.svg?style=flat" alt="Email the authors"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/katiebreivik/hush/actions/workflows/showyourwork.yml">
<img src="https://github.com/katiebreivik/hush/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/dag.pdf">
<img src="https://img.shields.io/badge/article-dag-blue.svg?style=flat" alt="Article graph"/>
</a>
</p>

## Welcome to _hush_!

This repository contains all of the code necessary to create the results and figures in [Thiele+2021](https://arxiv.org). In this study, Milky Way-like galaxies of double white dwarf binaries are simulated, which are orbiting in the frequency band observable by the future space-based gravitational wave detector [LISA](https://www.elisascience.org). This study is the first binary population simulation which incorporates the [metallicity-dependent binary fraction](https://iopscience.iop.org/article/10.3847/1538-4357/ab0d88). 

### Getting started:

To begin, clone the _hush_ repository to a local directory. Dependencies are laid out in [environment.yml](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/environment.yml). 

There are two levels of complexity at which you can run this pipeline for yourself. 

The first is to run the entire project pipeline start to finish, which can be done by simply moving into your `hush` directory and running the command `make pdf`. 

__!!WARNING!!: This process will generate 60+ GB of data on your local disk, and requires sufficient computational power to run the simulations.__ If you have questions about computational requirements, feel free to reach out to the authors by selecting the contact badge above.

The full pipeline will...
- download all the COSMIC + Ananke simulations which lie in [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo, and which form the base data for the simulations. These can also be downloaded on your own. Once these files are in your hush directory within `/src/data/`, the `make pdf` command will not re-download them each time.
- Simulate the double white dwarf populations and produce their corresponding data files.
- run the [scripts](https://github.com/katiebreivik/hush/tree/main/src/figures) which create the paper figure data which will be contained in a single HDF5 file as `results.hdf`
- generate figures and compile the paper whose text is contained in [ms.tex](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/src/ms.tex)

The second, and less computationally intense option is to work with just the data used to generate the paper figures. `results.hdf` can be downloaded from [this](https://zenodo.org/record/5722715#.YaA2Sy0ZPOQ) Zenodo. The other file needed, `FIRE.h5`, holds the positions and metallicities to create a metallicity-dependent star formation history and can be downloaded from [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo.

### Accessing the paper from this GitHub repo:

The "pdf" and "tarball" badges at the top of this README will compile and download the article PDF, or download a tarball containing all of the manuscript files, respectively.

This repository was created using the <img width = "150" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/> framework. if you run into any trouble, please check out the [showyourwork documentation](https://showyourwork.readthedocs.io) or reach out to us. If you think you've encountered a bug, head to the [issues page](https://github.com/rodluger/showyourwork/issues) and raise a new one if needed.
