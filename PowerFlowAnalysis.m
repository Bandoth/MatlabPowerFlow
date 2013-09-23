%% Header Section
% Joshua Hartwig, Christine Martin, Nan Ding
% Computational Methods of Power Systems
% Power Flow Analysis

% This program will perform a power flow analysis of a single conductor 
% three phase transmission line.
% This program assumes the three conductors are at equal height and 
% equally spaced.

function PowerFlowAnalysis()
%% Input Section
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
BusDistances = [12,     23,     34,     14,     45,     15;
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

%% Ybus Calculations

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
    %CurrentZs_pu = CurrentZs / 1600;
    
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

% Output Ybus Calculations onto screen
displayYbusCalcs();



%% Initializations for Gauss Seidel

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

% Fill unknowns matrix
Unknowns(1,:) = Pg_unknown;
Unknowns(2,:) = Qg_unknown;
Unknowns(3,:) = Pd_unknown;
Unknowns(4,:) = Qd_unknown;
Unknowns(5,:) = V_unknown;
Unknowns(6,:) = delta_unknown;

% Constant Definitions for easy access
ROW_Pg    = 1;
ROW_Qg    = 2;
ROW_Pd    = 3;
ROW_Qd    = 4;
ROW_V     = 5;
ROW_delta = 6;
BUSTYPE_Ref = 0;
BUSTYPE_Ctrl = 1;
BUSTYPE_Load = 2;

% Identify bus types
BusTypes = zeros(1, NumBuses);

for index = 1 : NumBuses
    if ((Unknowns(ROW_Pg, index) == 1) && (Unknowns(ROW_Qg, index) == 1))
        BusTypes(index) = BUSTYPE_Ref;
    elseif ((Unknowns(ROW_Qg, index) == 1) && (Unknowns(ROW_delta, index) == 1))
        BusTypes(index) = BUSTYPE_Ctrl;
    else
        BusTypes(index) = BUSTYPE_Load;
    end
end

% Output initial conditions onto screen
displayInitialConditions();

% Make initial guesses for unknown Vi for all buses
for index = 2 : NumBuses
    if Unknowns(ROW_V, index) == 1
        V_pu(index) = V_pu(1) * 0.9;
    end
end

% Make initial guesses for unknonw delta for all buses
for index = 2 : NumBuses
    if Unknowns(ROW_delta, index) == 1
        delta_rad(index) = -10 * 2 * pi / 360; % set to -10 degrees
    end
end

% Calculate inital Q for Control Buses
Q_estimate = zeros(1, NumBuses);
for index = 2 : NumBuses
    if BusTypes(index) == BUSTYPE_Ctrl
        SummationTerm = 0;
        for k = 1 : NumBuses
            [AngleRad, Mag] = cart2pol(real(Ybus(index, k)), imag(Ybus(index, k)));
            gamma_current = AngleRad;
            y_mag = Mag;
            SummationTerm = SummationTerm + y_mag * V_pu(k) * sin(delta_rad(k) - delta_rad(index) + gamma_current);
        end

        Q_estimate(index) = -V_pu(index) * SummationTerm;
    end
end

% Output Initial Guesses onto screen
displayInitialGuesses();



%% Gauss Seidel

BusConverged = 1;
BusNotConverged = 0;

ConvergenceStatus = zeros(1, NumBuses);
ConvergenceStatus(1) = BusConverged;

iterToConverge = zeros(1, NumBuses);

ConvergenceCondition = 0.000001;

iterCount = 1;
ConvergenceComplete = false;

while (~ConvergenceComplete)
    % Calculating values for next iteration
    for busIndex = 2 : NumBuses
        if (ConvergenceStatus(busIndex) == BusNotConverged)
            Pbus = Pg_pu(busIndex) - Pd_pu(busIndex);
            Qbus = Qg_pu(busIndex) - Qd_pu(busIndex);
            y_term = 1 / Ybus(busIndex, busIndex);
            
            %%% Control Bus Calculations
            if (BusTypes(busIndex) == BUSTYPE_Ctrl)
                % Update Qbus with control bus estimate from previous
                % iteration
                Qbus = Q_estimate(busIndex);
                
                % Calculate delta_rad(u+1)
                [Real_V, Imag_V] = pol2cart(delta_rad(busIndex), V_pu(busIndex));
                Current_V_phasor = Real_V + 1i * Imag_V;
                
                % Calculate Summation Term
                SummationTerm = 0;
                for k = 1 : NumBuses
                    if (k ~= busIndex)
                        [Real_V, Imag_V] = pol2cart(delta_rad(k), V_pu(k));
                        OtherBus_V_phasor = Real_V + 1i * Imag_V;
                        SummationTerm = SummationTerm + Ybus(busIndex, k) * OtherBus_V_phasor;
                    end
                end
                
                new_V_phasor = y_term * ((Pbus - 1i * Qbus) / conj(Current_V_phasor) - SummationTerm);
                [AngleRad, ~] = cart2pol(real(new_V_phasor), imag(new_V_phasor));
                last_delta_rad = delta_rad(busIndex);
                delta_rad(busIndex) = AngleRad;
                
                DELTA_delta = last_delta_rad - delta_rad(busIndex);
                
                % Calculate Q(u+1)
                % Calculate Summation Term
                SummationTerm = 0;
                for k = 1 : NumBuses
                    [AngleRad, Mag] = cart2pol(real(Ybus(busIndex, k)), imag(Ybus(busIndex, k)));
                    gamma_current = AngleRad;
                    y_mag = Mag;
                    SummationTerm = SummationTerm + y_mag * V_pu(k) * sin(delta_rad(k) - delta_rad(busIndex) + gamma_current);
                end
                
                last_Qestimate = Q_estimate(busIndex);
                Q_estimate(busIndex) = -V_pu(busIndex) * SummationTerm;
                
                DELTA_Qg = last_Qestimate - Q_estimate(busIndex);
                
                % Check Convergence Complete on current bus
                if ((abs(DELTA_delta) <= ConvergenceCondition) && (abs(DELTA_Qg) <= ConvergenceCondition))
                    Qg_pu(busIndex) = Q_estimate(busIndex) + Qd_pu(busIndex);
                    ConvergenceStatus(busIndex) = BusConverged;
                    iterToConverge(busIndex) = iterCount;
                end
                
                % End Control Bus calculations
                
            %%% Load Bus Calculations
            elseif (BusTypes(busIndex) == BUSTYPE_Load)
                % Calculate V(u+1) and delta_rad(u+1) simultaneously
                [Real_V, Imag_V] = pol2cart(delta_rad(busIndex), V_pu(busIndex));
                Current_V_phasor = Real_V + 1i * Imag_V;
                
                % Calculate Summation Term
                SummationTerm = 0;
                for k = 1 : NumBuses
                    if (k ~= busIndex)
                        [Real_V, Imag_V] = pol2cart(delta_rad(k), V_pu(k));
                        OtherBus_V_phasor = Real_V + 1i * Imag_V;
                        SummationTerm = SummationTerm + Ybus(busIndex, k) * OtherBus_V_phasor;
                    end
                end
                
                new_V_phasor = y_term * ((Pbus - 1i * Qbus) / conj(Current_V_phasor) - SummationTerm);
                [AngleRad, Mag] = cart2pol(real(new_V_phasor), imag(new_V_phasor));
                last_V = V_pu(busIndex);
                V_pu(busIndex) = Mag;
                last_delta_rad = delta_rad(busIndex);
                delta_rad(busIndex) = AngleRad;
                
                DELTA_V = last_V - V_pu(busIndex);
                DELTA_delta = last_delta_rad - delta_rad(busIndex);
                
                % Check Convergence Complete on current bus
                if ((abs(DELTA_V) <= ConvergenceCondition) && (abs(DELTA_delta) <= ConvergenceCondition))
                    ConvergenceStatus(busIndex) = BusConverged;
                    iterToConverge(busIndex) = iterCount;
                end

                % End Load Bus calculations
                
            end
        end
    end
    
    %%% Check if final iteration complete
    ConvergenceComplete = true;
    for convergeCounter = 1 : NumBuses
        if (ConvergenceStatus(convergeCounter) == BusNotConverged)
            ConvergenceComplete = false;
        end
    end
    
    if ~ConvergenceComplete
        iterCount = iterCount + 1;
    end
end

displayConvergenceResults();

%% Power Flows

[~, NumColumns] = size(YpYsMat);
YpYsMatSize = NumColumns;
Sbuses = zeros(NumBuses);

for i = 1 : NumBuses
    [Real, Imag] = pol2cart(delta_rad(i), V_pu(i));
    Vi_phasor = Real + 1i * Imag;

    for k = 1 : NumBuses
        [Real, Imag] = pol2cart(delta_rad(k), V_pu(k));
        Vk_phasor = Real + 1i * Imag;
        
        % Find correct Yp and Ys
        busCombination = i * 10 + k + 0i;
        busCombination2 = k * 10 + i + 0i;
        for j = 1 : YpYsMatSize
            if ((busCombination == YpYsMat(1, j)) || (busCombination2 == YpYsMat(1, j)))
                Yp_ik = YpYsMat(2, j);
                Ys_ik = YpYsMat(3, j);
                
                Sbuses(i, k) = Vi_phasor * (conj(Vi_phasor) - conj(Vk_phasor)) * conj(Ys_ik) + conj(Yp_ik) * V_pu(i)^2;
            end
        end
    end
end

Pg1_sum = 0;
Qg1_sum = 0;
for i = 1 : NumBuses
    Pg1_sum = Pg1_sum + real(Sbuses(1, i));
    Qg1_sum = Qg1_sum + imag(Sbuses(1, i));
end

Pg_pu(1) = Pg1_sum;
Qg_pu(1) = Qg1_sum;

disp(' ');
disp('Sbuses = ');
disp(num2str(Sbuses));
disp(' ');
disp('Pg_pu = ');
disp(num2str(Pg_pu));
disp(' ');
disp('Qg_pu = ');
disp(num2str(Qg_pu));


%% Power Losses

Plosses = zeros(NumBuses);
Qlosses = zeros(NumBuses);

for i = 1 : NumBuses
    for k = 1 : NumBuses
        Plosses(i, k) = real(Sbuses(i, k)) + real(Sbuses(k, i));
        Qlosses(i, k) = imag(Sbuses(i, k)) + imag(Sbuses(k, i));
    end
end

disp(' ');
disp('Plosses = ');
disp(num2str(Plosses));
disp(' ');
disp('Qlosses = ');
disp(num2str(Qlosses));


%% Net Power Checks


%[AngleRad, Mag] = cart2pol(real(Ybus(1,1)), imag(Ybus(1,1)));
%AngleDegrees = AngleRad * 360 / (2 * pi);
%[Real, Imag] = pol2cart(AngleRad, Mag);
%phasor = Real + 1i * Imag;


disp(' ');
disp('XXXXXXXXXXX');
%disp(num2str());

disp('Program End');

%% Display Functions

    function displayYbusCalcs()
        disp('**********************************************************');
        disp('******************Ybus Calculations***********************');
        disp('**********************************************************');
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
        disp('**********************************************************');
        disp('****************End Ybus Calculations*********************');
        disp('**********************************************************');
        disp(' ');
end

    function displayInitialConditions()
        disp(' ');
        disp('**********************************************************');
        disp('*************Initial Conditions of Buses******************');
        disp('**********************************************************');
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
        disp('Bus Types = ');
        disp(num2str(BusTypes));
        disp(' ');
        disp('**********************************************************');
        disp('************End Initial Conditions of Buses***************');
        disp('**********************************************************');
        disp(' ');
    end

    function displayInitialGuesses()
        disp(' ');
        disp('**********************************************************');
        disp('*******************Initial Guesses************************');
        disp('**********************************************************');
        disp(' ');
        disp('V_pu with Initial Guesses');
        disp(num2str(V_pu));
        disp(' ');
        disp('delta_rad with Initial Guesses');
        disp(num2str(delta_rad));
        disp(' ');
        disp('Qg_pu with Initial Guesses');
        disp(num2str(Qg_pu));
        disp(' ');
        disp('**********************************************************');
        disp('*****************End Initial Guesses**********************');
        disp('**********************************************************');
        disp(' ');
    end

    function displayConvergenceResults()
        disp(' ');
        disp('**********************************************************');
        disp('*****************Convergence Results**********************');
        disp('**********************************************************');
        disp(' ');
        disp('Convergence Condition');
        disp(num2str(ConvergenceCondition));
        disp(' ');
        disp('Number of Iterations for buses to convergence');
        disp(num2str(iterToConverge));
        disp(' ');
        disp('Number of Iterations for whole system to convergence');
        disp(num2str(iterCount));
        disp(' ');
        disp('Qg_pu approximations');
        disp(num2str(Qg_pu));
        disp(' ');
        disp('V_pu approximations');
        disp(num2str(V_pu));
        disp(' ');
        disp('delta_rad approximations');
        disp(num2str(delta_rad));
        disp(' ');
        
        delta_degrees_local = delta_rad ./ (2 * pi) .* 360;
        disp('delta approximations in degrees');
        disp(num2str(delta_degrees_local));
        disp(' ');
        
        disp('**********************************************************');
        disp('***************End Convergence Results********************');
        disp('**********************************************************');
        disp(' ');
    end

end
