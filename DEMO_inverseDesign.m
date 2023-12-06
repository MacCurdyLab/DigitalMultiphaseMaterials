%% DEMO_inverseDesign
%Lawrence Smith | lasm4254@colorado.edu

%this file is a demonstration of inverse design of digital multiphase
%materials. Raw data from empirical tests that characterize modulus and
%toughness for various inkjet composites are imported. Surface fits to the
%data are extracted. The user defines a target trajectory through material
%property 2-space. Then constrained nonlinear optimization is used to
%determine a set of points in material design space which most closely
%produce the target tuples of material properties.

clear; clc; close all

load rawDataset.mat

%Generate a surface fit for modulus data
sf1 = fit(Fs_Fr_1,log10(Modulus),'poly22');

%Generate a surface fit for toughness data
sf2 = fit(Fs_Fr_2,log10(Toughness),'poly23');

%Define a pair of vectors containing target modulus and toughness 
nseg = 5;
E_target = linspace(100,10,nseg)';
T_target = linspace(2.5,1.5,nseg)';

options = optimset('Display','none');

%set initial guess for volume fractions
xR = [0.4 0.4];

%initialize vectors for storing the recipe
F_soft = [];
F_rigid = [];

%loop over the entries in the target material property vectors
for i = 1:nseg
    
    %this matrix algebra problem constrains optimization solutions so that
    %the sum of rigid and soft fractions is between 0 and 0.6
    A = [1 1; -1 -1];
    b = [1,-0.6;];

    %these bounds constrain the fraction of rigid and soft components
    %each to the interval [0 1]
    lb = [0 0]';
    ub = [1 1]';

    %perform constrained optimization
    [xR,eR] = fmincon(@(x) errorFunc(x,E_target(i),T_target(i),sf1,sf2),xR,A,b,[],[],lb,ub,[],options);

    F_soft = [F_soft; xR(1)];
    F_rigid = [F_rigid; xR(2)];
    fprintf('Solved Tuple # %i of %i with residual %1.1e\n',i,nseg,eR)
    
end

F_liquid = 1-F_soft-F_rigid;
E_achieved = sf1(F_soft,F_rigid);
T_achieved = sf2(F_soft,F_rigid);

T = table(E_target,E_achieved,T_target,T_achieved,F_soft,F_rigid,F_liquid)

%plot reachable space


function e = errorFunc(x,targetE,targetT,sf1,sf2)

Fs = x(1);
Fr = x(2);

E = sf1(Fs,Fr);
T = sf2(Fs,Fr);

e = (E-targetE)^2 + (T-targetT)^2;

end
