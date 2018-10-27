%% demo for calculate the area and volume in all frames of Patient 2 data 
clear all; clc; close all;
PatientFileName = 'patient2.dcm';
thickness = 1; % cm 
outStepSize = 0.1; % cm 
sliceSpace = 0.1; % cm
[VolumeLV, ejectionFraction] = FindPatientLVEjectionFraction(...
    PatientFileName, thickness, outStepSize,sliceSpace);
