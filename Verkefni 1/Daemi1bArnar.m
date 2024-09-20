clear; close all;
% Constants
U = 14e3;         % Voltage across thing
p = 100e3;              % Pressure in Pa
A = 11.25;          % Constant for air in Pa^(-1) m^(-1)
B = 273.75;           % Constant for air in V Pa^(-1) m^(-1)
gamma = 1/p;       % Secondary electron emission coefficient
d_heild=3e-3;
width=0.015e-3;            % Thickness of dielectric material (initial)
d1 = (d_heild-width)/2;
d3=d1;
eps = [4,1,4];
d = [d1 width d3];


V = 0;
Ub = Inf;
k = 0;
widths = [];
breakdowns = [];
voltages = [];
while width<d_heild/20
    width = width + 1e-8;
    d1 = (d_heild-width)/2;
    d3=d1;
    Ub = breakdown(width,B,p,A,gamma);
    [~,E2,~] = EField(d,eps,U);
    V = E2*width;
    k = k+1;
    widths(k) = width;
    breakdowns(k) = Ub;
    voltages(k) = V;
    if abs(Ub-V)<1
        BDwidth = width;
        BDVoltage = V;
    end

end

plot(1000*widths,breakdowns,'DisplayName','Breakdown Voltage')
hold on
plot(1000*widths, voltages,'DisplayName','Actual Voltage Across')
plot(BDwidth*1000, BDVoltage, 'ro', 'MarkerSize', 10, 'DisplayName', 'Intersection');
text(BDwidth*1000, BDVoltage, [' X = ' num2str(BDwidth*1000),'mm'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
legend('show')
xlabel('Bubble width (mm)')
ylabel('Voltage [V]')
hold off

        
function Ub = breakdown(width,B,p,A,gamma)
    Ub=B*p*width/(log(A*p*width)-log(log(1+1/gamma)));   
end

function [E1,E2,E3] = EField(d,eps,U)
    Elist = zeros(3,1);
    for i = 1:3
        Elist(i)  = U/(d(1)*eps(i)/eps(1)+d(2)*eps(i)/eps(2)+d(3)*eps(i)/eps(3));
    end
    E1 = Elist(1);
    E2 = Elist(2);
    E3 = Elist(3);
end

    
    