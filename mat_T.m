function [T] = mat_T(nodes, points, k)

    T = zeros(length(points), length(nodes));
    for m=1:size(points,1)
        for n=1:length(nodes)
            dist = norm(points(m,:) - nodes(n,:));
            T(m,n) = exp(-1i*k*dist)/dist; 
        end
    end

end