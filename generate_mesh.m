function [Nodes] = generate_mesh(r,h,d)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generador de los nodos de conforman  %%% 
    %%% la malla de un casquete esférico con %%%
    %%% una separación mínima d entre nodos  %%%
    %%% r [radio]                            %%%
    %%% h [altura casquete]                  %%%
    %%% d [distancia mínima entre nodos]     %%%
    %%% Nodos = [Matriz nx3 de los nodos]    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    dth = d/r;
    theta = 0;
    Nodes = zeros(1,3);
    M = zeros(1,3);
    z=0;
    while z <= h
        z = r*(1 - cos(theta));
        Ra = r*sin(theta);
        phi = 0;
        for i=1:round(2*pi/(10*(d/Ra)))
            M(i,1) = Ra*cos(phi);
            M(i,2) = Ra*sin(phi);
            M(i,3) = z;
            phi = phi + 10*d/Ra;
        end
        Nodes = [Nodes; M];
        theta = theta + dth;
    end
end
