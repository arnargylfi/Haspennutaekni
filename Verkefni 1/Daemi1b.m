clear; close all;
% Constants
U = 14000;         % Breakdown voltage
p = 1e5;              % Pressure in bar
A = 112.5;          % Constant for air in Pa^(-1) m^(-1)
B = 2737;           % Constant for air in V Pa^(-1) m^(-1)
gamma = 1/p;       % Secondary electron emission coefficient
d=0.003;
d2=0.00002;            % Thickness of dielectric material (initial)
d1=(d-d2)/2
d3=d-d1-d2
eps1=4;
eps2=1;
eps3=4;
Err=1*10^-6;        % Minimum error acceptance
% Call the function to solve for d_min
%d_min = solve_min_diameter(Ub,p,A,B,gamma,d);
count=0;
% Convert pressure from bar to Pascals
%p = 1e2;  % 1 bar = 100,000 Pascals
E1=U/(d1+(eps1/eps2)*d2+(eps1/eps3)*d3)
E3=U/((eps3/eps1)*d1+(eps3/eps2)*d2+d3)
E2=U/((eps2/eps1)*d1+d2+(eps2/eps3)*d3) 
V2=E2*d2                                % The main result 
%Ub2=B*(p/1000)*d2*100/(log(A*(p/1000)*d2*100)-log(log(1+1/gamma)))   

% Since gamma is not defined properly, we will not iterate the Pachen
% equation, but instead use the values above and use the curve.

%Ub2=B*p*d2/(log(A*p*d2)-log(log(1+1/gamma)));   

%while abs(V2-Ub2)>Err & V2>0 & Ub2>0 & d2>0
% while V2>0 & d2>0 & count<1000
%     d2=d2-1*10^-5
%     d1=(d-d2)/2;
%     d3=d1;
%     E2=U/((eps2/eps1)*d1+d2+(eps2/eps3*d3))
%     V2=E2*d2
%     %Ub2=B*(p/1000)*d2*100/(log(A*(p/1000)*d2*100)-log(log(1+1/gamma)))   
%     count=count+1
% end

 
count
d2
E2
V2
%Ub2
% % Display the solution
%     fprintf('The minimum thickness d_min is: %.5f m\n', d2);
