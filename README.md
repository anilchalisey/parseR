README
================
Anil Chalisey

The parseR package provides tools for quality control of raw reads, read alignment, expression estimation, differential gene expression identification, and biological pathway and transcription motif factor enrichment analysis. It aims to provide an environment and pipeline that facilitates analysis of RNA-seq data in a straightforward and reproducible manner.

To our knowledge, this is the first package of this type that runs on both Linux/Unix-based systems and Windows 10. To run on the latter operating system, the package leverages the Windows-subsystem for linux (WSL) tool released in the Windows 10 Creator's Edition. In addition, this is one of the very few tools that runs a complete pipeline from read QC to motif enrichment analysis in a desktop or local server environment without the need to upload data to a web-based server.

A detailed guide to setting up WSL or Linux to run bioinformatic analyses may be found in the GitHub wiki for parseR [here](https://github.com/anilchalisey/parseR/wiki/Setting-up-WSL-Bash-on-Windows-10). A guide to using parseR may be found [here](https://github.com/anilchalisey/parseR/wiki/parseR-user-guide), and is also available as a vignette once the package is installed.
