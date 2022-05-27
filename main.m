%%
%%%%%%% DATOS NUBE PUNTOS %%%%%
x = 5;
y = 2;
z = 1;
zmin = 2;
n = 1;

%%%%%%% DATOS REFLECTANCIA %%%%%
rof = 998e-9; %*e-9
c = 1.5;
rop = 1064e-9;
cpl = 1.585;
cps = 0;

%%%%%% ATENUACIÓN %%%%%%%%%
at = 2.15e-5;
a = 5e-3;

%%%%%%%% DATOS CASQUETE %%%%%%%
d = c/20/2;                     %Distancia mínima entre nodos
r = 3;                          %Radio casquete
Rc = 14;                        %Radio de curvatura
h = Rc*(1-sqrt(1-(r/Rc)^2));    %Altura casquete
%%

tic
Nodes = generate_mesh(Rc,h,d);
% points = generate_random_points(x, y, z, zmin, n);
points = [1.5 0 14];
sen = load('sen_Alex.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRINGIR FRECUENCIAS 1
%%%% Espacio de análisis, 7mm, Tiempo de análisis=7*2/1.5=9.3us
%%%% Podemos retrasar 2*9/1.5= 12us
%%%% [~,indt]=min(abs(sen.t-9.3))=1396, reducimos a 2048 puntos
%%%%%%%%%%%%%%%%%%%%%% (nos quedamos con el intervaslo de tiempo
%%%%%%%%%%%%%%%%%%%%%% estrictamente necesario: 2048 puntos)

frecs = linspace(0.01,75,2048);
senc = sen.sen(1:2048);
pos = 0:0.016:3;

%%%%%%% ATENUACIÓN EN K %%%%%%%%%%%%
cc = c*(1+1i*c/2/pi*at*frecs);
k = 2*pi*frecs./cc;

SL2 = mod_reflectores_Alex(rof,c,rop,cpl,cps,a,frecs);
E = fft(sen.sen);
E(round(length(E)/2):end) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRINGIR FRECUENCIAS 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v imax]=max(abs(E(1:length(E)/2)));
freqind0=find(abs(E(1:length(E)/2))>v/5);
Ev=zeros(1,2048);
Ev((freqind0):max(freqind0)) = E(min(freqind0):max(freqind0)).*tukeywin(max(freqind0)-min(freqind0)+1,0.4)';
freqind = find(abs(Ev)>0);

sum_F = zeros(length(frecs),length(pos));

delete(gcp("nocreate"))
parpool("local", 4)

for j = 1:length(pos)
    Nodes(:,1) = Nodes(:,1) + 0.016;
    [num2str(j) '/' num2str(length(pos))];
    parfor i = 1:length(frecs)
        if find(freqind==i)
            E_ = E(i)*ones(length(Nodes),1);
            T = mat_T(Nodes, points, k(i));
            R = SL2(i)*eye(n);
            F = T.'*R*T*E_;
            sum_F(i,j) = sum(F);
        end
    end
    sum_FR(:,j)=sum_F(:,j).*exp(1i*2*pi*(frecs')*12);
end

sum_F_t = zeros(length(frecs), length(pos));
for i=1:length(pos)
    sum_F_t(:,i) = ifft(sum_FR(:,i));
end

imagesc(abs(sum_F_t))
toc
% scatter3(points(:,1), points(:,2), points(:,3), "filled")
% hold on
% scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3), "filled")



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
