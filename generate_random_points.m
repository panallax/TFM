function [points] = generate_random_points(x,y,z,zmin,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generador de n puntos uniformemente distribuidos      %%%
    %%% dentro de un espacio poliedrico centrado en el orígen %%%
    %%% de dimensiones dadas.                                 %%%
    %%% x,y,z son las dimensiones del espacio y zmin es la    %%%
    %%% distancia al centro de coordenadas.                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng('default')
    dims = [x y z];
    coords = zeros(2, length(dims));
    for i=2:length(dims)
        coords(1,i) = -dims(i)/2;
        coords(2,i) = dims(i)/2;
    end
    coords(:,1) = [-1 x+1]; 
    coords(:,3) = coords(:,3)+ (zmin + z/2);
    for i=1:length(coords)
        points(:,i) = coords(1,i) + (coords(2,i) - coords(1,i))*rand(n,1);
    end

end 