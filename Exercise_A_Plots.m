%% Solve ODE for exercise 4
%%% Solve for on state
tspan = [0, HR4on_tend];
x0 = HR4_on_initcond;  % CHANGE THIS VARIABLE FOR EXERCISES
D = D4_on;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_on1,HR4on_fit1] = ode23(@(t,x) odeFun2(t,x,D), tspan, x0);

%%% Solve for off state
tspan = [HR4on_tend, HR4off_tend];
x0 = HR4on_fit(end);    % set this equal tothe final value of ON state
D = D4_off;             % CHANGE THIS VARIABLE FOR EXERCISES
[t_off1,HR4off_fit1] = ode23(@(t,x) odeFun2(t,x,D), tspan, x0);

%%% Plot ODE solution for exercise 4
figure;
% Plot A normal
plot(t_on, HR4on_fit, 'blue', 'LineWidth', 3); hold on;
plot(t_off, HR4off_fit, '--', 'LineWidth', 3);
% Plot A reduced
plot(t_on1, HR4on_fit1, 'red', 'LineWidth', 3);
plot(t_off1, HR4off_fit1, '--', 'LineWidth', 3);

%%% Overlay the raw data for exercise 4
%scatter(ex4t, ex4hr);

%%% plot features
legend;
xlabel('Time (s)','FontSize',14);
ylabel('Heart Rate (bpm)','FontSize',14);
title('Exercise 4 - Fitted ODE');
grid on;
xlim([t_on(1), t_off(end)]);

%% Define ODE2
function [dhrdt] = odeFun2(~, x, D)
% Input ~: this is a dummy variable that Matlab requires as a placeholder
% Input x: this is A array
%       x(1) = heart rate (HR)  (beats/minute)
%       x(2) = oxygen demand (D(v,t)) (beats/minute)

% Define absolute HR min/max
HRmin = 50;
HRmax = 190;

% Define constants (non-normalized)
A = 2e-8;   % ( (beats/min)^(-3.38) ) / minute
B = 1.63;       % slope for leaving/approaching HR_min (dimensionless)
C = 1.75;       % slope for approaching/leaving HR_max (dimensionless)
E = 1.0;        % gives plateu shape (dimensionless)

%%% Define non-normalized ODE
% d/dt(hr) = A * [hr - HRmin]^B * [HRmax - hr]^C * [D(v,t) - hr]^E
HR = x(1);
dhrdt = A.*(( HR-HRmin ).^(B)) .* ((HRmax-x(1)).^C) .* ((D-HR).^E);
        
end