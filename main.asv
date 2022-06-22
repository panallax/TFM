%%
%%%%%%% DATOS NUBE PUNTOS %%%%%
x = 5;
y = 2;
z = 6;
zmin = 11;
ro = 10;
n_r = round(x*y*z*ro);

%%%%%%% DATOS REFLECTANCIA %%%%%
rof = 998e-9; %*e-9
c = 1.5;
rop = 1064e-9;
cpl = 1.585;
cps = 0;

%%%%%% ATENUACIÓN %%%%%%%%%
at = 2.15e-5;
a_r = 5e-3; %
a_t = 15e-3;
%%%%%%%% DATOS CASQUETE %%%%%%%
d = c/20/2;                     %Distancia mínima entre nodos
r = 3;       %                   %Radio casquete
Rc = 14;     %                   %Radio de curvatura
h = Rc*(1-sqrt(1-(r/Rc)^2));    %Altura casquete
%%%%%%%%% TEJIDOS E INTERFASES %%%%%%
x_t = -0.8:0.15:3.8;
y_t = -0.8:0.15:0.8;
z_1 = 10.8;
z_2 = 12.2;
ro_t = 500;
%%

tic
Nodes = generate_mesh(Rc,h,d);
% reflector_points = generate_random_points(x, y, z, zmin, n_r); 
[tissue_points,n_t] = tissue_generator(x_t,y_t,z_1,z_2, ro_t);
interfase_points_1 = interfase_generator(x_t,y_t,z_1);
interfase_points_2 = interfase_generator(x_t,y_t,z_2);
points = [tissue_points; interfase_points_1; interfase_points_2];
sub = squareform(pdist(points));
n = length(points);
% points = [tissue_points; reflector_points];
% n = n_r + n_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRINGIR FRECUENCIAS 1
%%%% Espacio de análisis, 7mm, Tiempo de análisis=7*2/1.5=9.3us
%%%% Podemos retrasar 2*9/1.5= 12us
%%%% [~,indt]=min(abs(sen.t-9.3))=1396, reducimos a 2048 puntos
%%%%%%%%%%%%%%%%%%%%%% (nos quedamos con el intervaslo de tiempo
%%%%%%%%%%%%%%%%%%%%%% estrictamente necesario: 2048 puntos)
sen = load('sen_Alex.mat');
frecs = linspace(0.01,150,2048);
senc = sen.sen(1:2048);
pos = 0:0.016:3;

%%%%%%% ATENUACIÓN EN K %%%%%%%%%%%%
cc = c*(1+1i*c/2/pi*at*frecs);
k = 2*pi*frecs./cc;

% SL2_r = mod_reflectores_Alex(rof,c,rop,cpl,cps,a_r,frecs);
SL2 = mod_reflectores_Alex(rof,c,rop,cpl,cps,a_t,frecs);

E = fft(senc);
E(round(length(E)/2):end) = 0;
%%%%%%%%%%%%%%%%%%% RESTRINGIR FRECUENCIAS 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v imax]=max(abs(E(1:length(E)/2)));
freqind0=find(abs(E(1:length(E)/2))>v/5);
Ev=zeros(1,2048);
Ev((freqind0):max(freqind0)) = E(min(freqind0):max(freqind0)).*tukeywin(max(freqind0)-min(freqind0)+1,0.4)';
freqind = find(abs(Ev)>0);

sum_F = zeros(length(frecs),length(pos));
sum_FR = zeros(length(frecs),length(pos));

% delete(gcp('nocreate'))
% numCores = feature('numcores');
% parpool(numCores)

for j = 1:length(pos)
    disp(j)
    Nodes(:,1) = Nodes(:,1) + 0.016;
%     [num2str(j) '/' num2str(length(pos))];
    for i = 1:length(frecs)
        if find(freqind==i)
%             E_ = E(i)*ones(length(Nodes),1);
            T = mat_T(Nodes, points, k(i));
            intern_dist = exp(-1i*k(i).*sub)./sub;
            intern_dist(1:1+size(intern_dist,1):end) = 0;
            R = SL2(i)^2.*intern_dist + SL2(i).*eye(n);
%             R = sparse(R);
            F = T.'*(R*sum(T,2));
            sum_F(i,j) = sum(F)*E(i);
        end
    end
    sum_FR(:,j)=sum_F(:,j).*exp(1i*2*pi*(frecs')*12);
end

sum_F_t = zeros(length(frecs), length(pos)); % 
for i=1:length(pos)
    sum_F_t(:,i) = ifft(sum_FR(:,i));
end
toc
imagesc(pos,linspace(zmin,z+zmin, 2048), abs(sum_F_t))
daspect([1 1 1])
colormap('jet')


% scatter3(points(:,1), points(:,2), points(:,3), "filled")
% hold on
% scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")




%%%%%% CAPTURADOR DE GIFS

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
% 

%%%%% CODIGO TENSORIAL %%%%%

% tic
% nodes_steps = repmat(points, 1, 1, length(pos));
% nodes_steps = permute(nodes_steps,[1 3 2]);
% vec = 1:length(pos);
% M = repmat(vec, length(points(:, 1)), 1);
% nodes_steps(:, :, 1) = points(:, 1)  + M*(pos(2)-pos(1));
% nodes_steps_ =  reshape(nodes_steps,[size(nodes_steps, 1)* size(nodes_steps, 2), size(nodes_steps, 3)]);
% dist = pdist2(nodes_steps_, Nodes, "euclidean");
% dist = reshape(dist, [size(nodes_steps, 1), size(nodes_steps, 2), size(Nodes, 1)]);
% dist = permute(dist,[1 3 2]);
% 
% T_total = exp(permute(k(:) .* permute(dist, [4 1 2 3]) * -1i, [2,3,4,1])) ./ dist;
% T_total = permute(T_total, [2,1,3,4]);
% I = eye(length(Nodes));
% SL2_total =  SL2(:) .* permute(I, [3 2 1]);
% SL2_total = permute(SL2_total, [2,3,1]);
% 
% for i = 1:length(pos)
%     sub_T = squeeze(T_total(:, :, i, :));
%     prod_st = pagemtimes(permute(sub_T, [2,1,3]), SL2_total);
%     prod_st = sum(squeeze(pagemtimes(prod_st, sum(sub_T,2))), 1);
%     sum_FR(:, i) = prod_st(1, :).*E.*exp(1i*2*pi*(frecs)*12);
% end
% 
% sum_F_t = zeros(length(frecs), length(pos));
% for i=1:length(pos)
%     sum_F_t(:,i) = ifft(sum_FR(:,i));
% end
% 
% imagesc(abs(sum_F_t))
% toc