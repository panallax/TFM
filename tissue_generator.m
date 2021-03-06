function [tissue, n] = tissue_generator(x_t,y_t,z1,z2, ro)
    tissue_coords = [x_t(1) x_t(end); y_t(1) y_t(end); z1 z2]'; 
    n = round(ro*(x_t(end) - x_t(1))*(y_t(end) - y_t(1))*(z2 - z1));
    tissue = zeros(n,3);
    for i=1:length(tissue_coords)
        tissue(:,i) = tissue_coords(1,i) + (tissue_coords(2,i) - tissue_coords(1,i))*rand(n,1);
    end


end