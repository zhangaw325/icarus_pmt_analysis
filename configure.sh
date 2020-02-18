#!/bin/bash

currentdir=$(pwd)

export PATH="$currentdir/build/:$PATH"
export LD_LIBRARY_PATH="$currentdir/build/:$LD_LIBRARY_PATH"
