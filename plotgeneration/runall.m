%RUNALL run all scripts and generate plots
% Script for generating all plots for paper
clear all; close all;
addpath('../')

girpccond
pccondnum

runhsdtests
runsdqtests
runhdqtests


generateplots
SuiteSparseTest

rmpath('../')
