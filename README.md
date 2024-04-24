<p align="center">
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
<a href="https://github.com/katiebreivik/hush/actions/workflows/build.yml">
<img src="https://github.com/katiebreivik/hush/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
</p>

## Welcome to _hush_!

This repository contains all of the code necessary to create the results and figures in [Thiele+2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...945..162T/abstract). In this study, we simulate populations of double white dwarf binaries which are orbiting in the frequency band observable by the future space-based gravitational wave detector [LISA](https://www.elisascience.org). Specificially, we investigate the effects of incorporating a [metallicity-dependent binary fraction](https://iopscience.iop.org/article/10.3847/1538-4357/ab0d88) for the first time.

### Getting started and usage:

To begin, clone the _hush_ repository to a local repository. Dependencies are detailed in [environment.yml](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/environment.yml). To build the article locally, you should have the latest version of **showyourwork!** installed (`pip install showyourwork`).

Running the command `showyourwork` in the `hush` directory will dowload the `results.hdf` data set from [this](https://zenodo.org/record/5722715#.YaA2Sy0ZPOQ) Zenodo and the `FIRE.h5` data set from [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo. Once that happens, the scripts to produce each figure will be run: `CEsep.py`, `LISA_SNR.py`, `Mc_vs_dist.py`, `PSD.py`, `form_eff.py`, `lisa_nums.py`, `model_comp.py`, and `sfh_vs_fb.py`.

<!---
You can also run the entire project pipeline start to finish, which can be done running the command `showyourwork` from the `hush` directory.

**!!WARNING!! This process will generate 60+ GB of data on your local disk, and requires ~20 GB of ram to run the simulations.** If you have questions about computational requirements, feel free to reach out to the authors using the contact badge above.

Running the `showyourwork` command will...

- download all the COSMIC + Ananke simulations which lie in [this](https://zenodo.org/record/5722451#.YZ152fHMLyg) Zenodo, and which form the base data for the simulations. These files will be downloaded in your hush directory within `/src/data/`. These can also be downloaded on your own. Once downloaded, the `showyourwork` command will not re-download them each time.
- simulate the double white dwarf populations and produce their corresponding data files which again will be in `/src/data/`.
- run the [scripts](https://github.com/katiebreivik/hush/tree/main/src/figures) to create the paper figure data, which will be contained in a single HDF5 file, `results.hdf`.
- generate figures and compile the paper, whose text is contained in [ms.tex](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/src/ms.tex)
--->

### Accessing the paper from this GitHub repo:

Instead of running the pipeline, you can also just compile the paper directly from this GitHub. It will reference our existing data - which was created using this pipeline - in their respective Zenodo sites. When clicked, the "pdf" and "tarball" badges at the top of this README will compile and download the article PDF, or download a tarball containing all of the manuscript files, respectively.

This repository was created using the [showyourwork](https://github.com/rodluger/showyourwork) framework. If you'd like to _show your own work_, check out the [showyourwork docs](https://showyourwork.readthedocs.io) here!
