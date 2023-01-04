tejido = load("/Users/alex/Desktop/resultados/tejido_12_30_20_5.npz.mat").sum_F_t;
cel = load("/Users/alex/Desktop/resultados/reflectores_11_25_20_31.npz.mat").sum_F_t;

sen_in = tejido + cel/4.3498e+03;

imagesc(pos,linspace(z0,z0+(2048-1)/fs/2*ca, 2048)-0.15, abs(sen_in))
daspect([1 1 1])
colormap('jet')

limit = 0.1;
rng(133993219)

% Ruido
SNR=80;
R=10^(-SNR/20)*(rand(size(sen_in))-0.5);
   
% Limitación
sen_lim=real(sen_in)+R*0;
sen_lim(sen_lim>limit)=limit;
sen_lim(sen_lim<-limit)=-limit;

% Filtrado en frecuencia
fs=150;
Rp=0.1;
Rs=60;
[Nfil,Wp]=ellipord([10 30]/fs*2,[5 35]/fs*2,Rp,Rs); %%% Frecuencias en MHz
[bf,af]=ellip(Nfil,Rp,Rs,Wp);     %%%% Pasabanda
sen_filt=filtfilt(bf,af,sen_lim);

% Envolvente
sen_out=reshape(abs(hilbert(sen_filt(:))),size(sen_filt,1),size(sen_filt,2));

imagesc(pos,linspace(z0,z0+(2048-1)/fs/2*ca, 2048)-0.15, sen_out)
daspect([1 1 1])
colormap('jet')