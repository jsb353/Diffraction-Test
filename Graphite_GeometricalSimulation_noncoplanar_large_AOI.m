% EXAMPLE SCRIPT TO GENERATE DIFFRACTION PATTERN FROM NONCOPLANAR DIFFRACTION GEOMETRY

% CLEAN WORKSPACE
home
clear
close all

% PATH TO DIFFWIZ LIBRARY
addpath(genpath('C:\Users\maher\Google drive\DIFFWIZ\'))
addpath(genpath('C:\Users\Cosmin\Desktop\Cr2AlC\'))
addpath(genpath('C:\Users\Cosmin\Desktop\Grand-Diffraction-Master\'))

% DEFINE LATTICE OR LOAD LATTICE FROM EXISTING STRUCTURE LIBRARY
load graphite.mat

% DEFINE DIRECTION OF CRYSTAL NORMAL
Lattice.Normal = [0 1 0 ]; 

% DEFINE X-RAYS
Probe.Type = 'electrons';
Probe.Energy = 60000; % [eV]
Probe.DiffractionGeometry = 'noncoplanar';
Probe.psi = 12; % This's the grazing angle

% DEFINE DETECTOR
Detector.Shape = 's';
Detector.Size = 1000; % diameter or length in mm
Detector.SpotFWHMx = 10; % Diffraction spot size on detector
Detector.SpotFWHMy = 10;
Detector.DistanceToSample = 200; % Sample-detector distance in mm
Detector.Offset = [0 0]; % offset from center of beam in mm

% MILLER INDICES TO LOOP OVER 
hkl = -6:6;

% MAIN FUNCTION
I = GeometricalSimulation2(Lattice, Probe, Detector, hkl,1);






