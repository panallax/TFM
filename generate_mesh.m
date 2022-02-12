d = 0.01;
r = 1;
h = 0.3;

h0 = r;
dth = d/r;
theta = 0;
Nodes = zeros(1,3);
M = zeros(1,3);
z=0;
while z <= h
    z = r*(1 - cos(theta))
    Ra = r*sin(theta);
    phi = 0;
    for i=1:round(2*pi/(d/Ra))
        M(i,1) = Ra*cos(phi);
        M(i,2) = Ra*sin(phi);
        M(i,3) = z;
        phi = phi + d/Ra;
    end
    Nodes = [Nodes; M];
    theta = theta + dth;
end

plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3))