def mat_T(nodes, points, k):
    dst = pdist2(nodes, points,"euclidean")';
    T = exp(-1i*k.*dst)./dst;
end