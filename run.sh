#!/bin/sh
rm -f Plots/*
julia test.jl && avconv -r 20 -i "Plots/plot%d.png"  Plots/output.avi
