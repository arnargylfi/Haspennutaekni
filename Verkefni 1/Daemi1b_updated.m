clear; close all;
% Constants
Ub = 14000;         % Breakdown voltage
p = 1;              % Pressure in bar
A = 112.5;          % Constant for air in Pa^(-1) m^(-1)
B = 2737;           % Constant for air in V Pa^(-1) m^(-1)
gamma = 0.01;       % Secondary electron emission coefficient
d=0.003;
d1=0.0012;
d2=0.0006;            % Thickness of dielectric material (initial)
d3=d2;
eps1=4;
eps2=1;
eps3=4;
Err=1*10^-7;        % Minimum error acceptance
% Call the function to solve for d_min
%d_min = solve_min_diameter(Ub,p,A,B,gamma,d);
count=0;
% Convert pressure from bar to Pascals
p = p * 1e5;  % 1 bar = 100,000 Pascals

E2=Ub/((eps2/eps1)*d1+d2+(eps2/eps3)*d3)
V2=E2*d2
Ub2=B*p*d2/(log(A*p*d2)-log(log(1+1/gamma)));   

while abs(V2-Ub2)>Err & d2>0
    d2=d2-2*10^-8
    d1=(d-d2)/2;
    d3=d1;
    E2=Ub/((eps2/eps1)*d1+d2+(eps2/eps3*d3))
    V2=E2*d2
    Ub2=B*p*d2/(log(A*p*d2)-log(log(1+1/gamma)))   
    count=count+1
end

d2
count

% % Display the solution
%     fprintf('The minimum thickness d_min is: %.5f m\n', d2);
