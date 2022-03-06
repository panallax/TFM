function [T] = mat_T(nodes, points, k)

    T = zeros(length(points), length(nodes));
    for m=1:length(points)
        for n=1:legth(nodes)
            dist = norm(points(m,:) - nodes(n,:));
            T(m,n) = exp^(-1i*k*dist)/dist; 
        end
    end

end