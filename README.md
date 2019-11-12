# SAVERX

R package for transfer learning of scRNA-seq denoising. Take a look at our free [SAVER-X web-server](https://singlecell.wharton.upenn.edu/saver-x/) for the transfer learning online computation! We also encourage you to read our [pre-print manucript](https://www.biorxiv.org/content/10.1101/457879v2) for more information. You can also refer to our earlier denoising method [SAVER](http://github.com/mohuangx/SAVER).

## Updates
11/12/2019: fix the bug in SAVER-X by requiring the reticulate package to be updated to version 1.13 

8/14/2019: SAVER-X now accepts either a data file or a matrix object in R as input, and creates a new folder under current working directory to store intemediate and final results. It returns the file name storing the denoised results.

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

Our input can be either a data file or a dense/sparse matrix. The data file can be '.txt', '.csv' or '.rds' file. The '.rds' file can store either a matrix or a sparse matrix of class 'dgCMatrix'. 
As a toy example, you may download one of our demo datasets on the web server, shekhar_downsampled.csv, the down-sampled mouse retina data from [here](https://www.dropbox.com/sh/kctbw41kdh6jmnb/AAAO5Icu97Ep6uoWFdHRKIcMa?dl=0). SAVER-X will generate a new folder based on the current time to store intemediate files and final results and return the file name storing the denoised matrix. 
### SAVER-X without pretraining

```
library(SAVERX)
file <- saverx("./testdata/shekhar_downsampled.csv")
denoised.data <- readRDS(file)
```

### SAVER-X with a pretrained model
For the demo dataset, we have a pre-trained model for the mouse retina, please download the file, mouse_Retina.hdf5, and you may save it in './mouse_retina.hdf5'
```
library(SAVERX)
file <- saverx("./testdata/shekhar_downsampled.csv", data.species = "Mouse", 
use.pretrain = T, pretrained.weights.file = "./mouse_retina.hdf5", model.species = "Mouse")
denoised.data <- readRDS(file)
```

For both cases, object 'denoised.data' is your final denoised matrix. When dealing with large datasets, you can set 'is.large.data = T' to reduce RAM. Set 'clearup.python.session = F' if you plan to run SAVER-X more than once in your R session. 

