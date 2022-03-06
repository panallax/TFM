%%

%%%%%%%% DATOS CASQUETE %%%%%%%
d = 0.01;  %Distancia mínima entre nodos
r = 1;     %Radio casquete
h = 0.3;   %Altura casquete

%%%%%%% DATOS NUBE PUNTOS %%%%%
x = 1;
y = 1;
z = 1;
zmin = 2;
n = 1000;

%%%%%%% DATOS REFLECTANCIA %%%%%
k = 1;
Ew = 1;
%%

Nodes = generate_mesh(r,h,d);
points = generate_random_points(x, y, z, zmin, n);
T = mat_T(Nodes, points, k);


scatter3(points(:,1), points(:,2), points(:,3), "filled")
hold on
scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")