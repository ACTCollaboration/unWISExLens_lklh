#!/usr/bin/env bash

# script from Joe Zuntz
if [ -d data ]
then
    echo Data already downloaded
elif ! command -v wget &> /dev/null
then
    echo wget not installed. Please obtain it to download the data.
else
    wget https://portal.nersc.gov/project/act/act_x_unWISE_xcorr+3x2pt/data_unWISExLens.tar.gz
    tar -zxvf data_unWISExLens.tar.gz
    rm data_unWISExLens.tar.gz
fi