function [dhrdt] = odeFun3(~, x, D)
% Input ~: this is a dummy variable that Matlab requires as a placeholder
% Input x: this is A array
%       x(1) = heart rate (HR)  (beats/minute)
%       x(2) = oxygen demand (D(v,t)) (beats/minute)

% Define absolute HR min/max
%Change Hrmin from 40 -> 30 or HrMax from 185 -> 175
HRmin = 40;
HRmax = 181; %181
 
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
