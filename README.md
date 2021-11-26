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

This repository contains all of the code necessary to create the results and figures in [Thiele+2021](https://arxiv.org). In this study, we simulate populations of double white dwarf binaries which are orbiting in the frequency band observable by the future space-based gravitational wave detector [LISA](https://www.elisascience.org). Specificially, we investigate the effects of incorporating a [metallicity-dependent binary fraction](https://iopscience.iop.org/article/10.3847/1538-4357/ab0d88) for the first time. 

### Getting started and usage:

To begin, clone the _hush_ repository to a local repository. Dependencies are detailed in [environment.yml](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/environment.yml). 

There are two levels of complexity at which you can run the paper pipeline for yourself. 

The first is to run the entire project pipeline start to finish, which can be done by simply moving into your `hush` directory and running the command `make pdf`. 

__!!WARNING!! This process will generate 60+ GB of data on your local disk, and requires ~20 GB of ram to run the simulations.__ If you have questions about computational requirements, feel free to reach out to the authors using the contact badge above.

Running the `make pdf` command will...
- download all the COSMIC + Ananke simulations which lie in [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo, and which form the base data for the simulations. These files will be downloaded in your hush directory within `/src/data/`. These can also be downloaded on your own. Once downloaded, the `make pdf` command will not re-download them each time.
- simulate the double white dwarf populations and produce their corresponding data files which again will be in `/src/data/`.
- run the [scripts](https://github.com/katiebreivik/hush/tree/main/src/figures) to create the paper figure data, which will be contained in a single HDF5 file, `results.hdf`.
- generate figures and compile the paper, whose text is contained in [ms.tex](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/src/ms.tex)

The second, and less computationally intense option is to skip the galactic simulations and work with just the already-simulated data which we used to generate the paper figures. `results.hdf` can be downloaded from [this](https://zenodo.org/record/5722715#.YaA2Sy0ZPOQ) Zenodo. The other figure file needed, `FIRE.h5`, holds the data to implement a metallicity-dependent star formation history, and can be downloaded from [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo. The paper figures can be generated using `CEsep.py`, `LISA_SNR.py`, `Mc_vs_dist.py`, `PSD.py`, `form_eff.py`, `lisa_nums.py`, `model_comp.py`, and `sfh_vs_fb.py`. 

### Accessing the paper from this GitHub repo:

Instead of running the pipeline, you can also just compile the paper directly from this GitHub. It will reference our existing data - which was created using this pipeline - in their respective Zenodo sites. When clicked, the "pdf" and "tarball" badges at the top of this README will compile and download the article PDF, or download a tarball containing all of the manuscript files, respectively.

This repository was created using the [showyourwork](https://github.com/rodluger/showyourwork) framework. If you'd like to _show your own work_, check out the [showyourwork docs](https://showyourwork.readthedocs.io) here!
