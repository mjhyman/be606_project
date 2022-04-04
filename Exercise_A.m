%% Mapping from filenames to exercise number
clear; clc; close all;
% 'Exercise A,B' = exercise 2 on/off. Velocity = 13.4 Km/hr
% 'Exercise C,D' = exercise 3 on/off. velocity = 14.4 Km/hr
% 'Exercise E,F' = exercise 4 on/off. velocity = 15.7 Km/hr

% Exercise 2
HR2_on = importdata('Exercise A.txt');
HR2_off = importdata('Exercise B.txt');
% Exercise 3
HR3_on = importdata('Exercise C.txt');
HR3_off = importdata('Exercise D.txt');
% Exercise 4
HR4_on = importdata('Exercise E.txt');
HR4_off = importdata('Exercise F.txt');

%%% Retrieve initial heart rates for solving ODE for each exercise
HR2_on_initcond = HR2_on(1,2);
HR3_on_initcond = HR3_on(1,2);
HR4_on_initcond = HR4_on(1,2);
HR2_off_initcond = HR2_off(1,2);
HR3_off_initcond = HR3_off(1,2);
HR4_off_initcond = HR4_off(1,2);

%%% Define oxygen demand (these are taken from table 1)
D2_on = 156;
D2_off = 72;

D3_on = 166;
D3_off = 72;

D4_on = 175;
D4_off = 70; % The paper says this is non-constant. Maybe try adjusting this?

%%% Define time span for each exercise (taken from S1_dataset)
HR2on_tend = find_tend(HR2_on);
HR2off_tend = find_tend(HR2_off);

HR3on_tend = find_tend(HR3_on);
HR3off_tend = find_tend(HR3_off);

HR4on_tend = find_tend(HR4_on);
HR4off_tend = find_tend(HR4_off);

%% HR for each exercise
% Exercise 2
a = HR2_on; 
b = HR2_off;
% Exercise 3
c = HR3_on;
d = HR3_off;
% Exercise 4
e = HR4_on;
f = HR4_off;

%resize matrix to remove zeros
at = a(1:962,1);
ahr = a(1:962,2);

bt = b(1:705,1);
bhr = b(1:705,2);

ct = c(1:922,1);
chr = c(1:922,2);

dt = d(1:734,1);
dhr = d(1:734,2);

et = e(1:933,1);
ehr = e(1:933,2);

ft = f(1:786,1);
fhr = f(1:786,2);

%add time from exercise time vector to recovery vector
btnew = bt + at(end,1);
dtnew = dt+ct(end,1);
ftnew = ft + et(end,1);

%append time and heart rate vectors
ex2t = [at; btnew];
ex2hr = [ahr; bhr];
ex3t = [ct; dtnew];
ex3hr = [chr; dhr];
ex4t = [et; ftnew];
ex4hr = [ehr; fhr];

%% Solve ODE for exercise 2
%%% Solve for on state
tspan = [0, HR2on_tend];
x0 = HR2_on_initcond;  % CHANGE THIS VARIABLE FOR EXERCISES
D = D2_on;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_on,HR2on_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Solve for off state
tspan = [HR2on_tend, HR2off_tend];
x0 = HR2on_fit(end);    % set this equal to the final value of ON state
D = D2_off;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_off,HR2off_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Plot ODE solution for exercise 2
figure;
% Plot ON state
%plot(t_on, HR2on_fit, 'LineWidth', 3); hold on;
% Overlay OFF state
%plot(t_off, HR2off_fit, 'LineWidth', 3);

%%% Overlay the raw data for exercise 2
%scatter(ex2t, ex2hr);

%%% plot features
xlabel('Time (s)','FontSize',14);
ylabel('Heart Rate (bpm)','FontSize',14);
title('Exercise 2 - Fitted ODE');
grid on;
xlim([t_on(1), t_off(end)]);

%% Solve ODE for exercise 3
%%% Solve for on state
tspan = [0, HR3on_tend];
x0 = HR3_on_initcond;  % CHANGE THIS VARIABLE FOR EXERCISES
D = D3_on;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_on,HR3on_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Solve for off state
tspan = [HR3on_tend, HR3off_tend];
x0 = HR3on_fit(end);    % set this equal to the final value of ON state
D = D3_off;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_off,HR3off_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Plot ODE solution for exercise 3
figure;
% Plot ON state
%plot(t_on, HR3on_fit, 'LineWidth', 3); hold on;
% Overlay OFF state
%plot(t_off, HR3off_fit, 'LineWidth', 3);

%%% Overlay the raw data for exercise 3
%scatter(ex3t, ex3hr);

%%% plot features
xlabel('Time (s)','FontSize',14);
ylabel('Heart Rate (bpm)','FontSize',14);
title('Exercise 3 - Fitted ODE');
grid on;
xlim([t_on(1), t_off(end)]);

%% Solve ODE for exercise 4
%%% Solve for on state
tspan = [0, HR4on_tend];
x0 = HR4_on_initcond;  % CHANGE THIS VARIABLE FOR EXERCISES
D = D4_on;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_on,HR4on_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Solve for off state
tspan = [HR4on_tend, HR4off_tend];
x0 = HR4on_fit(end);    % set this equal to the final value of ON state
D = D4_off;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_off,HR4off_fit] = ode23(@(t,x) odeFun(t,x,D), tspan, x0);

%%% Plot ODE solution for exercise 4
figure;
% Plot ON state
plot(t_on, HR4on_fit, 'LineWidth', 3); hold on;
% Overlay OFF state
plot(t_off, HR4off_fit, 'LineWidth', 3);

%%% Overlay the raw data for exercise 4
%scatter(ex4t, ex4hr);

%%% plot features
legend
xlabel('Time (s)','FontSize',14);
ylabel('Heart Rate (bpm)','FontSize',14);
title('Exercise 4 - Fitted ODE');
grid on;
xlim([t_on(1), t_off(end)]);

%% Eigenvalue 2
%{
% Define constants (non-normalized)
A = 3.217e-8;   % ( (beats/min)^(-3.38) ) / minute
B = 1.63;       % slope for leaving/approaching HR_min (dimensionless)
C = 1.75;       % slope for approaching/leaving HR_max (dimensionless)
E = 1.0;        % gives plateu shape (dimensionless)
% Define HRmin and HRmax from data
HRmin = 40;
HRmax = 185;
%%% Define D(v) between min/max HR
D = linspace(HRmin, HRmax, 100);
% Define eigenvalue for HRmin < D < HRmax
Lbound = -A.*((D-HRmin).^B).*((HRmax-D)).^C;
% Plot eigenvalue
figure;
plot(D, Lbound);
xlabel('Oxygen Demand (D)');
ylabel('Eigenvalue');
title({'Eigenvalue','HRmin < D < HRmax'});
%%% Define D(v) > max HR
D = linspace(HRmax, (HRmax+10), 100);
% Define eigenvalue for D > HRmax
Lunbound = -A.*((D-HRmin).^B).*((HRmax-D)).^C;
% Define real and imaginary components
Lreal = real(Lunbound);
Limag = imag(Lunbound);
% Plot real and imaginary components of eigenvalue
figure;
scatter(Lreal, Limag, 'g*');
xlabel('Real');
ylabel('Imaginary');
title({'Eigenvalue','D > HRmax'});
% Plot eigenvalue
% figure;
% plot(D, Lunbound);
% xlabel('Oxygen Demand (D)');
% ylabel('Eigenvalue');
% title({'Eigenvalue','D > HRmax'});
%}
%% Define ODE
function [dhrdt] = odeFun(~, x, D)
% Input ~: this is a dummy variable that Matlab requires as a placeholder
% Input x: this is A array
%       x(1) = heart rate (HR)  (beats/minute)
%       x(2) = oxygen demand (D(v,t)) (beats/minute)

% Define absolute HR min/max
HRmin = 50;
HRmax = 190;

% Define constants (non-normalized)
A = 3.217e-8;   % ( (beats/min)^(-3.38) ) / minute
B = 1.63;       % slope for leaving/approaching HR_min (dimensionless)
C = 1.75;       % slope for approaching/leaving HR_max (dimensionless)
E = 1.0;        % gives plateu shape (dimensionless)

%%% Define non-normalized ODE
% d/dt(hr) = A * [hr - HRmin]^B * [HRmax - hr]^C * [D(v,t) - hr]^E
HR = x(1);
dhrdt = A.*(( HR-HRmin ).^(B)) .* ((HRmax-x(1)).^C) .* ((D-HR).^E);
        
end
%% Find the last recording for each dataset
function [tend] = find_tend(recording)
% data: the raw dataset

% this will be the first zero element
[~,ind] = min(recording(:,2));

% the previous index is the final data point
ind = ind-1;

% take respective timestamp
tend = recording(ind,1); 

end
