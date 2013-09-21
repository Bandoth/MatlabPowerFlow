% Joshua Hartwig, Christine Martin, Nan XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Computational Methods of Power Systems
% Power Flow Analysis

% This program will perform a power flow analysis of a single conductor 
% three phase transmission line.
% This program assumes the three conductors are at equal height and 
% equally spaced.

%%%%%%%%%%%%% INPUT SECTION %%%%%%%%%%%%%
% Set characteristics of the line in terms of ohms and meters
LineRadius = 0.015189;      % radius of conductors in meters
Rl = 0.00006723253;         % resistance per length in ohms/meter
D_betCond = 9;              % distances between conductors in meters
H_toGround = 30;            % conductor height above ground in meters

% Set characteristics of the circuit
NumBuses = 5;               % Number of buses in the circuit
% Bus distances in meters in 2nd row.  1st row for names unused.  Must be
% in same order as incidence matrix.
BusDistances = [12,     23,     34,     41,     45,     51;
                350000, 305000, 250000, 290000, 300000, 320000];
            
Sbase = 100;    % Sbase in MVA
Vbase = 400;    % Vbase in kV

% Set power generation, voltage, and load information of buses
% Use 99999 for unknown values
Pg    = [99999,    550,    0,      0,      0];     % Pg in MW
Qg    = [99999,    99999,  0,      0,      0];     % Qg in MVAR
Pd    = [0,        350,    290,    225,    320];   % Pd in MW
Qd    = [0,        220,    195,    120,    205];   % Qd in MVAR
V_mag = [400,      395,    99999,  99999,  99999]; % V magnitude in kV
delta = [0,        99999,  99999,  99999,  99999]; % Voltage angle in degrees
%%%%%%%%%%%%% END INPUT SECTION %%%%%%%%%%%%%

% This program will read in a csv file containing an incidence matrix.
% The first row must be the names of the branches in order of the incidence
% matrix.  You must list them as "yp12, ys12, yp21, yp23, ys23, yp32" etc.
% The following rows are the incidence matrix itself.
fileID = fopen('Incidence Matrix.csv','r');    % Open incidence matrix as var "fileID"
% Import as strings in cell data for names
IncidenceMatNames = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','delimiter',',');
fclose(fileID); % Close file
%celldisp(IncidenceMatNames) % Debug for checking strings
IncidenceMatrix = csvread('Incidence Matrix.csv', 1, 0, [1,0,5,17]);    % Import incidence matrix

% Calculate equivalent distance and radius (GMD and GMR?)
Deq = nthroot((D_betCond.^2 * (D_betCond * 2)), 3);
Radeq = LineRadius * exp(-1 / 4);

% Begin calculating Zl
L = (2 * 10.^-7) * log(Deq / Radeq);    % Calculate inductance of line per meter
Omega = 2 * pi * 60;    % Calculate omega at 60Hz
Xl = Omega * L;         % Calculate Xl
Zl = Rl + 1i * Xl;      % Calculate Zl for the transmission line

% Begin calculating Yl
Epsilon0 = 8.85 * 10.^-12;  % Set constant epsilon0
HToReflection = H_toGround * 2; % Set height to conductor reflection
% Calculate heights to neighboring reflections
H12 = nthroot((HToReflection.^2 + D_betCond.^2), 2);
H23 = nthroot((HToReflection.^2 + D_betCond.^2), 2);
H31 = nthroot((HToReflection.^2 + (D_betCond * 2).^2), 2);
% Calculate capacitance to neutral
Cn = (2 * pi * Epsilon0) / (log(Deq / LineRadius) - log(nthroot((H12 * H23 * H31), 3) / HToReflection));
Yl = 1i * Omega * Cn;   % Calculate Yl for the transmission line

% Calculate gamma
Gamma = sqrt(Zl * Yl);

% Calculate Yp and Ys for each bus connection and set up for Ybranch
YpYsMat(1,:) = BusDistances(1,:);
numBranches = size(IncidenceMatNames);
YbranchVector = zeros(1, numBranches(1,2));
for index = 1 : NumBuses+1
    len = BusDistances(2, index);
    
    CurrentYp = Yl * len * (tanh(Gamma * len / 2) / (Gamma * len));
    CurrentYp_pu = CurrentYp * 1600;
    
    CurrentZs = Zl * len * (sinh(Gamma * len) / (Gamma * len));
    CurrentZs_pu = CurrentZs / 1600;
    
    CurrentYs = 1 / CurrentZs;
    CurrentYs_pu = CurrentYs * 1600;
    
    YpYsMat(2,index) = CurrentYp_pu;
    YpYsMat(3,index) = CurrentYs_pu;
    
    YbranchVector(1, (index-1)*3+1) = CurrentYp_pu;
    YbranchVector(1, (index-1)*3+2) = CurrentYs_pu;
    YbranchVector(1, (index-1)*3+3) = CurrentYp_pu;
end

% Create Ybranch
Ybranch = diag(YbranchVector);

% Create bus admittance matrix Ybus
Ybus = IncidenceMatrix * Ybranch * IncidenceMatrix';

% Convert inputs into per unit and radians
% Simultaneously find unknowns
Pg_pu = zeros(size(Pg));
Qg_pu = zeros(size(Qg));
Pd_pu = zeros(size(Pd));
Qd_pu = zeros(size(Qd));
V_pu  = zeros(size(V_mag));
delta_rad = zeros(size(delta));

Pg_unknown = zeros(size(Pg));
Qg_unknown = zeros(size(Qg));
Pd_unknown = zeros(size(Pd));
Qd_unknown = zeros(size(Qd));
V_unknown  = zeros(size(V_mag));
delta_unknown = zeros(size(delta));

Unknowns = zeros(6, NumBuses);

for index = 1 : NumBuses
    if Pg(index) ~= 99999
        Pg_pu(index) = Pg(index) / Sbase;
        Pg_unknown(index) = 0;
    else
        Pg_pu(index) = 99999;
        Pg_unknown(index) = 1;
    end
    
    if Qg(index) ~= 99999
        Qg_pu(index) = Qg(index) / Sbase;
        Qg_unknown(index) = 0;
    else
        Qg_pu(index) = 99999;
        Qg_unknown(index) = 1;
    end
    
    if Pd(index) ~= 99999
        Pd_pu(index) = Pd(index) / Sbase;
        Pd_unknown(index) = 0;
    else
        Pd_pu(index) = 99999;
        Pd_unknown(index) = 1;
    end
    
    if Qd(index) ~= 99999
        Qd_pu(index) = Qd(index) / Sbase;
        Qd_unknown(index) = 0;
    else
        Qd_pu(index) = 99999;
        Qd_unknown(index) = 1;
    end
    
    if V_mag(index) ~= 99999
        V_pu(index) = V_mag(index) / Vbase;
        V_unknown(index) = 0;
    else
        V_pu(index) = 99999;
        V_unknown(index) = 1;
    end
    
    if delta(index) ~= 99999
        delta_rad(index) = delta(index) * (2*pi) / 360;
        delta_unknown(index) = 0;
    else
        delta_rad(index) = 99999;
        delta_unknown(index) = 1;
    end
end

Unknowns(1,:) = Pg_unknown;
Unknowns(2,:) = Qg_unknown;
Unknowns(3,:) = Pd_unknown;
Unknowns(4,:) = Qd_unknown;
Unknowns(5,:) = V_unknown;
Unknowns(6,:) = delta_unknown;



% Output onto screen
disp('Line radius');
disp(num2str(LineRadius));
disp(' ');
disp('Distance between conductors');
disp(num2str(D_betCond));
disp(' ');
disp('Height to Ground');
disp(num2str(H_toGround));
disp(' ');
disp('Rl = ');
disp(num2str(Rl));
disp(' ');
disp('Number of buses = ');
disp(num2str(NumBuses));
disp(' ');
disp('Bus Distances');
disp(num2str(BusDistances));
disp(' ');
disp('Deq');
disp(num2str(Deq));
disp(' ');
disp('Radius eq');
disp(num2str(Radeq));
disp(' ');
disp('L = ');
disp(num2str(L));
disp(' ');
disp('Omega = ');
disp(num2str(Omega));
disp(' ');
disp('Xl = ');
disp(num2str(Xl));
disp(' ');
disp('Zl = ');
disp(num2str(Zl));
disp(' ');
disp('Epsilon0 = ');
disp(num2str(Epsilon0));
disp(' ');
disp('Height to conductor reflections = ');
disp(num2str(HToReflection));
disp(' ');
disp('H12 = ');
disp(num2str(H12));
disp(' ');
disp('H23 = ');
disp(num2str(H23));
disp(' ');
disp('H31 = ');
disp(num2str(H31));
disp(' ');
disp('Cn = ');
disp(num2str(Cn));
disp(' ');
disp('Yl = ');
disp(num2str(Yl));
disp(' ');
disp('Gamma = ');
disp(num2str(Gamma));
disp(' ');
disp('Incidence Matrix = ');
disp(num2str(IncidenceMatrix));
disp(' ');
disp('Yp and Ys = ');
disp(num2str(YpYsMat));
disp(' ');
disp('YbranchVector = ');
disp(num2str(YbranchVector));
disp(' ');
disp('Ybranch = ');
disp(num2str(Ybranch));
disp(' ');
disp('Ybus = ');
disp(num2str(Ybus));
disp(' ');
disp('Pg Original = ');
disp(num2str(Pg));
disp(' ');
disp('Qg Original = ');
disp(num2str(Qg));
disp(' ');
disp('Pd Original = ');
disp(num2str(Pd));
disp(' ');
disp('Qd Original = ');
disp(num2str(Qd));
disp(' ');
disp('V original = ');
disp(num2str(V_mag));
disp(' ');
disp('Delta in degrees Original = ');
disp(num2str(delta));
disp(' ');
disp('Pg_pu = ');
disp(num2str(Pg_pu));
disp(' ');
disp('Qg_pu = ');
disp(num2str(Qg_pu));
disp(' ');
disp('Pd_pu = ');
disp(num2str(Pd_pu));
disp(' ');
disp('Qd_pu = ');
disp(num2str(Qd_pu));
disp(' ');
disp('V_pu = ');
disp(num2str(V_pu));
disp(' ');
disp('Delta in radians = ');
disp(num2str(delta_rad));
disp(' ');
disp('Unknowns = ');
disp(num2str(Unknowns));
disp(' ');
disp('XXXXXXXXXXX');
%disp(num2str());
disp(' ');
disp('Program End');

% polar and rectangular conversions
% [Real, Imag] = pol2cart(pi/5, 1);
% [Angle, Mag] = cart2pol(Real, Imag);
