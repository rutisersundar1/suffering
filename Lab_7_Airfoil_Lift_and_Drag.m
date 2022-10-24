clc
clear all

%Varibles
data = xlsread('Airfoil_Data.xlsx');
historic = xlsread('Historical Data for NACA0015.csv');

x_historic = historic(:,1);
velocity_ratio = historic(:,2);

chord = 0.2; %m
x_coord = data(1:39,1) * chord;
airfoil_x = x_coord;
airfoil_x(end+1) = 1; %Add point at (1, 0) to close the airfoil for display
y_coord = data(1:39,2) * chord;
airfoil_y = y_coord;
airfoil_y(end+1) = 0; %Added to close the airfoil for display

zero_deg_pressure = data(1:39,3);
five_deg_pressure = data(1:39,4);
ten_deg_pressure = data(1:39,5);
fifteen_deg_pressure = data(1:39,6);

zero_deg_data = readtable('NACA_0015_airfoil_data\0015 0 degrees.txt');
five_deg_data = readtable('NACA_0015_airfoil_data\0015 5 degrees.txt');
ten_deg_data = readtable('NACA_0015_airfoil_data\0015 10 degrees.txt');
fifteen_deg_data = readtable('NACA_0015_airfoil_data\0015 15 degrees.txt');
Polar_data_comp = readtable('NACA_0015_airfoil_data\Polar data 0015.txt');
Historical_data_comp = readtable('NACA_0015_airfoil_data\Polar data 0015 R3.2.txt');

%Airfoil Plot
% figure(1)
% plot(airfoil_x,airfoil_y,'r')
% axis equal

%Calculating Pressure Coefficients
% Historic

for i = 1:length(x_historic)
    C_p_Historic(i) = ((1-velocity_ratio(i)));
end

%Calculate pressure coefficients at varying AoA
C_p_0deg = zero_deg_pressure'     ./ 0.621607;
C_p_5deg = five_deg_pressure'     ./ 0.631065;
C_p_10deg = ten_deg_pressure'     ./ 0.625189;
C_p_15deg = fifteen_deg_pressure' ./ 0.629732;

figure(2)
plot(x_coord, C_p_0deg)
hold on
%plot(x_coord(21:end, C_p_0deg(21:end)))
plot(zero_deg_data.x(1:50,1),zero_deg_data.Cpv(1:50,1),'--r')
hold off
set(gca, 'YDir', 'reverse')
legend('Experimental','Computational')
title('Pressure Coefficients at 0 degrees of attack')
xlabel('x/c')
ylabel('Pressure Coefficient')

figure(3)
plot(x_coord,C_p_5deg)
hold on
plot(five_deg_data.x,five_deg_data.Cpv,'r')
hold off
set(gca, 'YDir', 'reverse')
legend('Experimental','Computational')
title('Pressure Coefficients at 5 degrees of attack')
xlabel('x/c')
ylabel('Pressure Coefficient')


figure(4)
plot(x_coord,C_p_10deg)
set(gca, 'YDir','reverse')
legend('Experimental')
title('Pressure Coefficients at 10 degrees of attack')
xlabel('x/c')
ylabel('Pressure Coefficient')


figure(5)
plot(x_coord,C_p_15deg)
set(gca, 'YDir','reverse')
legend('Experimental')
title('Pressure Coefficients at 15 degrees of attack')
xlabel('x/c')
ylabel('Pressure Coefficient')
%%%%
% Lift and Drag Coefficients for angles of attack
% Historic
angles = linspace(1,26,1);
C_N_Historic = -trapz(C_p_Historic,x_historic);
C_A_Historic = trapz(C_p_Historic,x_historic);
%trapz(x_historic, C_p_Historic)
C_L_Historic = -C_A_Historic*sind(angles)+C_N_Historic*cosd(angles);
C_D_Historic = C_A_Historic*cosd(angles)+C_N_Historic*sind(angles);

%Prepare sections of airfoil for segmenting data
%Nightmare structure so that we can give wing segments to functions...
%Yes, all of this is awful and could be made way more efficient. 
segments = struct(...
    "chord", chord, ...
    "upper_x", x_coord(1:20), ...
    "upper_x_flip", flip(x_coord(1:20), 1),...
    "lower_x", x_coord(21:end), ...
    "lower_y", y_coord(21:end), ...
    "upper_y", y_coord(1:20),...
    "upper_y_flip", flip(y_coord(1:20), 1));
%Segment C_p data into front/back sections

Cpseg0deg = segmentCp(C_p_0deg, segments);
Cpseg5deg = segmentCp(C_p_5deg, segments);
%trapz(segments.back_l_y, Cpseg5deg.back_l)
Cpseg10deg = segmentCp(C_p_10deg, segments);
Cpseg15deg = segmentCp(C_p_15deg, segments);

%Calculate C_N, C_L, C_A, C_D from C_p data
Coeffs0deg = calculateCs(Cpseg0deg, 0, segments);
Coeffs5deg = calculateCs(Cpseg5deg, 5, segments);
Coeffs10deg = calculateCs(Cpseg10deg, 10, segments);
Coeffs15deg = calculateCs(Cpseg15deg, 15, segments);
%Coeffs: .N, .L, .D, .A

%Check proper front/back segmenting
%hold on
%plot(front_x, front_y, 'Color', [1 0 0]);
%plot(back_x, back_y, 'Color', [0 1 0]);
%plot(back_x_lower, back_y_lower, 'Color', [0 0 1]);
%hold off

% C_l vs Alpha Plot
figure(6)

plot(Polar_data_comp.Var1(7:37), Polar_data_comp.Var2(7:37))
hold on

plot(0, Coeffs0deg.L,'r*')
plot(5, Coeffs5deg.L,'b*')
plot(10, Coeffs10deg.L,'g*')
plot(15, Coeffs15deg.L,'c*')
hold off
legend('Computational','Experimental 0 degrees','Experimental 5 degrees','Experimental 10 degrees','Experimental 15 degrees')
title('Coefficient of Lift versus degree of attack')
xlabel('Angle of attack (degrees)')
ylabel('Lift Coefficient')

% C_l vs C_d
figure(7)

plot(Polar_data_comp.Var2(7:37),Polar_data_comp.Var3(7:37))
hold on 
plot(Coeffs0deg.L, Coeffs0deg.D,'rx')
plot(Coeffs5deg.L, Coeffs5deg.D,'bx')
plot(Coeffs10deg.L, Coeffs10deg.D,'gx')
plot(Coeffs15deg.L, Coeffs15deg.D,'cx')
hold off

legend('Computational','Experimental 0 degrees','Experimental 5 degrees','Experimental 10 degrees','Experimental 15 degrees')
title('Coefficient of Lift versus Coefficient of Drag')
xlabel('Lift Coefficient')
ylabel('Drag Coefficient')

% Historical Plots
% C_l vs Alpha Plot
% figure(8)
% plot(Historical_data_comp.Var1(7:37),Historical_data_comp.Var2(7:37))
% hold on 
% plot(Historical_data_comp.Var1(7:37),C_L_Historic)
% hold off
% % C_l vs C_d
% figure(9)
% plot(Historical_data_comp.Var3(7:37),Historical_data_comp.Var2(7:37))
% hold on
% plot(C_D_Historic,C_L_Historic);

%Segment Cp data into upper/lower sections of wing
function Cp = segmentCp(CpData, seg)
    Cp = struct( ...
        "upper_flipped", flip(CpData(:, 1:20), 2),...
        "upper", CpData(:, 1:20), ...
        "lower", CpData(:, 21:end));
end

%Calculate coefficients N A L D from Cp data
function Coeffs = calculateCs(Cp, alpha, seg) 
    C_N = (trapz(seg.lower_x, Cp.lower)...
        - trapz(seg.upper_x_flip, Cp.upper_flipped)) / seg.chord;
    
    C_A = (1/seg.chord) * (trapz(seg.lower_y, Cp.lower)...
        - trapz(seg.upper_y_flip, Cp.upper_flipped)) / seg.chord;
    
    C_L = -C_A * sind(alpha) + C_N * cosd(alpha);

    C_D = C_A * cosd(alpha) + C_N * sind(alpha);

    Coeffs = struct( ...
        "N", C_N, ...
        "A", C_A, ...
        "L", C_L, ...
        "D", C_D);
end