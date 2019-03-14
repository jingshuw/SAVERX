# SAVERX

R package for transfer learning of scRNA-seq denoising. Take a look at our free [SAVER-X web-server](https://singlecell.wharton.upenn.edu/saver-x/) for the transfer learning online computation! We also encourage you to read our [pre-print manucript](https://www.biorxiv.org/content/10.1101/457879v2) for more information.

## Installation

First, install the supporting Python package sctransfer. See the source code of the package [here](http://github.com/jingshuw/sctransfer)

```
pip install sctransfer
```

Next, open R and install the R package SAVERX
```
library(devtools)
install_github("jingshuw/SAVERX")
```

## Basic Usage

Our current pre-trained models can be downloaded [here](https://www.dropbox.com/sh/4u22cfuswcfcwvu/AAC6CgsO7dvQSNInTF0wWMDva?dl=0)
