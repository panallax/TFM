%%

%%%%%%%% DATOS CASQUETE %%%%%%%
d = c/40;  %Distancia mínima entre nodos
r = 1;     %Radio casquete
h = 0.3;   %Altura casquete

%%%%%%% DATOS NUBE PUNTOS %%%%%
x = 1;
y = 1;
z = 1;
zmin = 2;
n = 1000;

%%%%%%% DATOS REFLECTANCIA %%%%%
rof = 998e-9; %*e-9
c = 1.5;
rop = 1064e-9;
cpl = 1.585;
cps = 0;
a = 5e-3;
%%
Nodes = generate_mesh(r,h,d);
points = generate_random_points(x, y, z, zmin, n);

pos = 0:0.16:3;
frecs = linspace(0,150,4096);

SL2 = mod_reflectores_Alex(rof,c,rop,cpl,cps,a,frecs);
E = fft(sen);
sum
for j=1:length(pos)
    for i=1:length(frecs)
        E_ = E(i)*ones(length(Nodes),1);
        k = 2*pi/(frecs(i)*c);
        T = mat_T(Nodes, points, k);
        a = SL2(i)*eye(n);
        R = a*T;
        F = T'*R*T*E_;
        sum_F(i,j) = sum(F);
    end
end

scatter3(points(:,1), points(:,2), points(:,3), "filled")
hold on
scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")