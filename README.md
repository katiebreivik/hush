<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/katiebreivik/hush/actions/workflows/showyourwork.yml">
<img src="https://github.com/katiebreivik/hush/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/dag.pdf">
<img src="https://img.shields.io/badge/article-dag-blue.svg?style=flat" alt="Article graph"/>
</a>
<a href="https://github.com/katiebreivik/hush/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
<a href="mailto:sarahgthiele@gmail.com?cc=kbreivik@flatironinstitute.org">
      <img src="https://img.shields.io/badge/contact-authors-blueviolet.svg?style=flat" alt="Email the authors"/>
</a>
</p>

## Welcome to _hush_ !

This repository contains all of the code necessary to create the results and figures in [Thiele+2021](https://arxiv.org). In this study, a Milky Way-like galaxies of double white dwarf binaries are simulated, which are orbiting in the frequency band observable by the future space-based gravitational wave detector [LISA](https://www.elisascience.org), and is the first binary population simulation which incorporates the [metallicity-dependent binary fraction](https://iopscience.iop.org/article/10.3847/1538-4357/ab0d88). 

### Getting started:

To begin, clone the _hush_ repository to a local directory. Dependencies are laid out in [environment.yml](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/environment.yml). 

There are two levels at which to run this pipeline for yourself. The first is from scratch: to run the entire pipeline start to finish, simply run the command `make pdf`. __!! WARNING: This process will generate over 60 GB of data on your local disk, and requires sufficient computational power to run the simulations !!__ If you have questions about computational requirements, feel free to reach out to the authors by selecting

- download all the COSMIC + Ananke simulations which lie in [this Zenodo](https://zenodo.org/record/5722451#.YZ152fHMLyg) which form the base data for the simulations. 
- run the pipeline to simulate the double white dwarf populations and produce their corresponding data files
- run the [scripts](https://github.com/katiebreivik/hush/tree/main/src/figures) which create the paper figures
- generate figures and compile the paper whose text is contained in [ms.tex](https://github.com/katiebreivik/hush/blob/1eaf321cc5bc97dbc260139181cf2618bc16f833/src/ms.tex)

Your new repository is set up and ready to go. Click on the badges at the top to take you to the compiled article PDF or to a tarball containing all of the manuscript files. Both the PDF and the tarball are automatically updated every time you push changes to this repo; note that builds usually take a few minutes (or more, depending on what you're doing).

The first thing you might want to do is customize the `src/ms.tex` file, which is currently just filled with placeholder text. You should also delete the current placeholder figure scripts in the `src/figures` directory, and add the scripts needed to build your own figures. If your workflow has external dependencies (which it most likely will), you must add them to the `environment.yml` file so `showyourwork` can build the paper from scratch. See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments) for details. Finally, change or edit the `LICENSE` as needed and replace the text in this `README.md` with some useful information for the reader!

If you run into any trouble, please check out the [showyourwork documentation](https://showyourwork.readthedocs.io). If you think you've encountered a bug, please check out the [issues page](https://github.com/rodluger/showyourwork/issues) and raise a new one if needed.
