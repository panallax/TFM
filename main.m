%%
%%%%%%% DATOS NUBE PUNTOS %%%%%
x = 1;
y = 1;
z = 1;
zmin = 2;
n = 1;

%%%%%%% DATOS REFLECTANCIA %%%%%
rof = 998e-9; %*e-9
c = 1.5;
rop = 1064e-9;
cpl = 1.585;
cps = 0;
a = 5e-3;

%%%%%%%% DATOS CASQUETE %%%%%%%
d = c/40;  %Distancia mínima entre nodos
r = 1;     %Radio casquete
h = 0.3;   %Altura casquete
%%
tic
Nodes = generate_mesh(r,h,d);
points = generate_random_points(x, y, z, zmin, n);
points = [1.5 0 1];
sen = load('sen_Alex.mat');
pos = 0:0.16:3;
frecs = linspace(0.01,150,4096);

SL2 = mod_reflectores_Alex(rof,c,rop,cpl,cps,a,frecs);
E = fft(sen.sen);
sum_F = zeros(length(frecs),length(pos));

for j=1:length(pos)
    Nodes(:,1) = Nodes(:,1) + 0.16;
    for i=1:length(frecs)
        E_ = E(i)*ones(length(Nodes),1);
        k = 2*pi/(frecs(i)*c);
        T = mat_T(Nodes, points, k);
        R = SL2(i)*eye(n);
        F = T.'*R*T*E_;
        sum_F(i,j) = sum(F);
    end
end

sum_F_t = zeros(length(frecs), length(pos));
for i=1:length(pos)
    sum_F_t(:,i) = ifft(sum_F(:,i));
end

imagesc(abs(sum_F_t))
toc
scatter3(points(:,1), points(:,2), points(:,3), "filled")
hold on
scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")



% nImages = length(pos);
% 
% fig = figure;
% for idx = 1:nImages
%     Nodes(:,1) = Nodes(:,1) + 0.16;
%     scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")
%     drawnow
%     frame = getframe(fig);
%     im{idx} = frame2im(frame);
% end
