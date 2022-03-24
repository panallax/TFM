function [SL2]=mod_reflectores_Alex(rof,c,rop,cpl,cps,a,frecs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% REFLECTORES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Modelo simplificado: se asume una difracción puramente radial desde el
%%%%%%%%%% dispersor, con el valor correspondiente a dispersión directa (¿o
%%%%%%%%%% retrodispersión? a comprobar)


%%%%%%%%%% Los parámetros de entrada son:
%%%%%%%%%% rof y rop, densidades del fluido y el dispersor en kg/mm^3
%%%%%%%%%% c,cpl y cps son la velocidad del fluido, la velocidad
%%%%%%%%%% longitudinal del dispersor y la velocidad transversal del
%%%%%%%%%% dispersor en mm/s
%%%%%%%%%% a es el radio del dispersor en mm
%%%%%%%%%% frecs es el vector de las frecuencias

%%%%%%%%%% El parámetro de salida SL2 es la presión dispersada compleja en
%%%%%%%%%% función de la frecuencia.




%% Modelo celular (leucocito) ¡¡¡ATENCION, ESTOS TRES PARÁMETROS VAN FUERA DE LA FUNCIÓN, EN EL PROGRAMA PRINCIPAL!!!
% rof = 998
% c = 1500
% rop = 1064
% cpl = 1585
% cps = 0
% a = 5e-6
% frec = [0:75] MHz
% -----------
% rop=1064*1e-9;
% cpl=1.585;
% cps=0*0.270*(1-1i*8e-10*frecs);

  
%%%%%% Magnitudes auxiliares
w=2*pi*frecs;
sigma=0.5*(1-2*(cps./cpl).^2)./(1-(cps./cpl).^2);
if sigma==0.5
  sigma=0.499999999;
  cps=0.000000001;
end
 
k=w./c;
k1=w./cpl;
k2=w./cps;

x=k*a;
x1=k1*a;
x2=k2*a;

nfo=5;         %%% Orden máximo de la serie de polinomios.
clear ln2
SL2=zeros(1,length(frecs));
for n=0:nfo
    jn_1=sqrt(pi/2 ./x1).*besselj(n+0.5,x1);
    jnmas1_1=sqrt(pi/2 ./x1).*besselj(n+1.5,x1);
    jnp_1=n./x1.*jn_1-jnmas1_1;
    jns_1=(n^2-n-x1.^2)./x1.^2 .*jn_1+2*jnmas1_1./x1;

    jn_2=sqrt(pi/2 ./x2).*besselj(n+0.5,x2);
    jnmas1_2=sqrt(pi/2 ./x2).*besselj(n+1.5,x2);
    jnp_2=n./x2.*jn_2-jnmas1_2;
    jns_2=(n^2-n-x2.^2).*jn_2./x2.^2+2*jnmas1_2./x2;


    jn=sqrt(pi/2 ./x).*besselj(n+0.5,x);
    jnmas1=sqrt(pi/2 ./x).*besselj(n+1.5,x);
    jnp=n./x.*jn-jnmas1;

    nn=sqrt(pi/2 ./x).*bessely(n+0.5,x);
    nnmas1=sqrt(pi/2 ./x).*bessely(n+1.5,x);
    nnp=n./x.*nn-nnmas1;

    hn=jn+1i*nn;
    hnp=jnp+1i*nnp;
    
    
    Mn=2*(jn_1-x1.*jnp_1)./((n^2+n-2)*jn_2+jns_2.*x2.^2);
    A11=-rof./x.*(x1.*jnp_1+Mn*n*(n+1).*jn_2);
    A12=-hnp;
    A21=-rop./(1-sigma).*(sigma.*jn_1+(1-2*sigma).*(-jns_1+n*(n+1)*Mn.*(jn_2-x2.*jnp_2)./x1.^2));
    A22=-hn;
    ln2(n+1,:)=(-1)^n*(2*n+1)*(A11.*jn-A21.*jnp)./(A11.*A22-A21.*A12)./k;
    SL2=SL2+ln2(n+1,:);

end
