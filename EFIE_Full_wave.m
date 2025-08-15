%% MATLAB Implementation for EFIE-Based Scattering Field Computation (Full-Wave Method)
clear all; close all; clc;

PI = 3.141592653589793;
Epsilon_0 = 8.854e-12;
Mu_0 = 12.56637061e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0);
f = 970e6;
Lambda = c / f;
DeltaX = Lambda / 4.0;
Omega = 2.0 * PI * f;
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0);
GrossStep = 10.0;
GrossNoSteps = 384;
NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);

Xsource = 0.0;
Ysource = 442.0;
H = 2.4;

fileID = fopen('X.04', 'r');
data = textscan(fileID, '%f %f', 385);
fclose(fileID);
X_terrain = data{1};
Y_terrain = data{2};

J = complex(zeros(1, NoLinesubs));
Et = complex(zeros(1, NoLinesubs));
Z = complex(zeros(NoLinesubs, NoLinesubs));
E_inc = complex(zeros(NoLinesubs, 1));

fprintf('Constructing impedance matrix and incident field...\n');
for p = 1:NoLinesubs
    for q = 1:NoLinesubs
        if p == q
            R_self = R_p_q(p, p+1, X_terrain, Y_terrain, DeltaX, GrossStep);
            Z(p,p) = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
                     (R_self - 1i*(2*R_self/pi)*log(1.781*Beta_0*R_self/(4*exp(1))));
        else
            Rpq = R_p_q(q, p, X_terrain, Y_terrain, DeltaX, GrossStep);
            Z(p,q) = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
                     (besselj(0, Beta_0*Rpq) - 1i*bessely(0, Beta_0*Rpq)) * ...
                     R_p_q(q, q+1, X_terrain, Y_terrain, DeltaX, GrossStep);
        end
    end
    R_source = R_source_p(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep);
    E_inc(p) = -(Beta_0^2/(4*Omega*Epsilon_0)) * ...
               (besselj(0, Beta_0*R_source) - 1i*bessely(0, Beta_0*R_source));
end

fprintf('Solving linear system for surface current...\n');
J = Z \ E_inc;

fprintf('Performing electric field calculation...\n');
for idx = 1:NoLinesubs
    Et(idx) = complex(0, 0);
    for n = 1:idx
        Rn = R_surf_obs(n, idx, X_terrain, Y_terrain, DeltaX, GrossStep, H);
        Z_val = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
               (besselj(0, Beta_0*Rn) - 1i*bessely(0, Beta_0*Rn));
        Et(idx) = Et(idx) + J(n) * R_p_q(n, n+1, X_terrain, Y_terrain, DeltaX, GrossStep) * Z_val;
    end
    R_obs = R_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H);
    Ei = -(Beta_0^2/(4*Omega*Epsilon_0)) * ...
        (besselj(0, Beta_0*R_obs) - 1i*bessely(0, Beta_0*R_obs));
    Et(idx) = Ei - Et(idx);
end

fprintf('Outputting surface current data...\n');
fileID = fopen('SurfaceCurrent(FullWave).dat', 'w');
for idx = 1:NoLinesubs
    fprintf(fileID, '%.6f\t%.6e\n', ...
        x_coord(idx, DeltaX), abs(J(idx)));
end
fclose(fileID);

fprintf('Outputting electric field data...\n');
fileID = fopen('ElectricField(FullWave).dat', 'w');
for idx = 1:NoLinesubs
    R_obs = R_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H);
    dB_field = 20*log10(abs(Et(idx))/sqrt(R_obs));
    fprintf(fileID, '%.6f\t%.6f\n', x_coord(idx, DeltaX), dB_field);
end
fclose(fileID);

figure('Units','normalized','Position',[0.1 0.3 0.4 0.4])
plot(x_coord(1:NoLinesubs, DeltaX), abs(J), 'LineWidth',1.5)
xlabel('Distance Along Surface (m)')
ylabel('|J| (A/m)')
title('Surface Current Distribution (Full Wave Solution)')
grid on
saveas(gcf, 'SurfaceCurrent(FullWave).fig')

figure('Units','normalized','Position',[0.5 0.3 0.4 0.4])
R_obs = arrayfun(@(idx) R_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H), 1:NoLinesubs);
dB_fields = 20*log10(abs(Et)./sqrt(R_obs));
plot(x_coord(1:NoLinesubs, DeltaX), dB_fields, 'LineWidth',1.5)
xlabel('Distance Along Surface (m)')
ylabel('Electric Field (dB)')
title('Scattered Field Distribution (Full Wave Solution)')
grid on
saveas(gcf, 'ElectricField(FullWave).fig')

function x = x_coord(a, DeltaX)
    x = a * DeltaX;
end

function y = y_coord(a, DeltaX, GrossStep, Y_terrain)
    pos = a * DeltaX;
    section = floor(pos / GrossStep) + 1;
    if section >= length(Y_terrain)
        error('Terrain data index out of range: position %.2f exceeds data range', pos);
    end
    dx = pos - (section-1)*GrossStep;
    prop = dx / GrossStep;
    y = Y_terrain(section) + prop*(Y_terrain(section+1)-Y_terrain(section));
end

function R = R_source_p(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep)
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    R = sqrt((Xsource - xp)^2 + (Ysource - yp)^2);
end

function R = R_source_obs(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H)
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    R = sqrt((Xsource - xp)^2 + (Ysource - yp - H)^2);
end

function R = R_p_q(p, q, X_terrain, Y_terrain, DeltaX, GrossStep)
    if p == q
        error('Error calculating distance for the same node: p=%d, q=%d', p, q);
    end
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    xq = x_coord(q, DeltaX);
    yq = y_coord(q, DeltaX, GrossStep, Y_terrain);
    R = sqrt((xq - xp)^2 + (yq - yp)^2);
end

function R = R_surf_obs(n, idx, X_terrain, Y_terrain, DeltaX, GrossStep, H)
    xn = x_coord(n, DeltaX);
    yn = y_coord(n, DeltaX, GrossStep, Y_terrain);
    xobs = x_coord(idx, DeltaX);
    yobs = y_coord(idx, DeltaX, GrossStep, Y_terrain) + H;
    R = sqrt((xobs - xn)^2 + (yobs - yn)^2);
end
