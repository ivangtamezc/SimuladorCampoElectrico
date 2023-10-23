clc;

% declaración de constantes
Q = input("Escriba la carga (placas): ");
Qg = input("Escriba la carga (globulo): ");
Xi = -4;
Yi = -2;
sep = 0.1;
dis = 10;

longP = input("Escriba la longitud positiva:");
anchoP = input("Escriba el ancho positivo: ");

longN = input("Escriba la longitud negativa: ");
anchoN = input("Escriba el ancho negativo: ");

[xe, ye] = meshgrid(-8:0.4:8);
hold on
U = zeros(size(xe));
V = zeros(size(ye));

% declaración de cargas
qP = Q;
qN = -Q;

% cálculo de diferencial de carga
sigma = 4;

dAp = (longP * anchoP);
dAn = (longN * anchoN);

dQp = dAp * sigma;
dQn = dAn * sigma;

% definición de las coordenadas para las láminas

nY = Yi:sep:Yi+longP;
nY2 = Yi: sep: Yi+longN;

mX = Xi*ones(size(nY));
mX2 = (Xi+dis)*ones(size(nY2));

if longP > longN
    mtx = mX2;
else
    mtx = mX;
end

% cálculo y suma de los campos eléctricos

for i = 1:size(mtx, 2)
    
    K = 8.99 * 10^9;

    Rx = xe - mX(i);
    Ry = ye - nY(i);

    R = sqrt(Rx.^2 + Ry.^2).^3;
    
    Ex = K .* qN .* Rx ./ R;
    Ey = K .* qN .* Ry ./ R;

    Rx = xe - mX2(i);
    Ry = ye - nY2(i);

    R = sqrt(Rx.^2 + Ry.^2).^3;

    Ex = Ex + K.*qP.*Rx./R;
    Ey = Ey + K.*qP.*Ry./R;

    E = sqrt(Ex.^2 + Ey.^2);

    U = U + Ex./E;
    V = V + Ey./E;

end

% gráficos

g = quiver(xe, ye, U, V, 'AutoScaleFactor', 0.6);
set(g, 'color', [1 0 0], 'LineWidth', 1.2);
plot(mX, nY, 'ob', 'LineWidth',5);
plot(mX2, nY2, 'or', 'LineWidth',5);
grid on

% globulos

Ti = 0;
deltat = 0.01;
x0 = 0;
y0 = 4;
v0 = 0;
tiempo = 100;
g = 9.8; 

vx = [0]; 
vy = [-10]; 

tiempo1 = [Ti]; 
tiempo2 = [Ti]; 
 
px = [x0];  
py = [y0]; 

RxInicial = x0 - Xi;
RyInicial = y0 - Yi;

rInicial = sqrt((RxInicial)^2 + (RyInicial)^2)^2;
iteraciones = tiempo/deltat;

masa = 10000000;

ax = [-(E(x0 + 1, y0 + 1) * Qg)/masa];
ay = [(- (masa * g)-(E(x0 + 1, y0 + 1) * Qg)) / masa]; 

for i = 2:(iteraciones+1) 
    vy = [vy, vy(i-1)+deltat*ay(i-1)];
    vx = [vx, vx(i-1)+deltat*ax(i-1)];

    py = [py, py(i-1)+vy(i-1)*deltat+(1/2)*ay(i-1)*(deltat^2)];
    px = [px, px(i-1)+vx(i-1)*deltat+(1/2)*ax(i-1)*(deltat^2)];

    if py(i) > 10 || py(i) < -10
        break
    elseif px(i) > 4 || px(i) < -4
        break
    end

    Xidx = round(abs(px(i-1)));
    Yidx = round(abs(py(i-1)));

    ax = [ax, ax+ -(E(Xidx + 1, Yidx + 1) * Qg)/masa];
    ay = [ay, ay + ((masa * g)-(E(Xidx + 1, Yidx + 1) * Qg)) / masa];

end 

plot(px, py, '-o','MarkerIndices',1:5:length(py))