function [area,vol,k] = get_polyhedron(ptcloud)

    % This function takes a point cloud, finds its convex hull, and returns
    % the area and volume, along with the indices of the points in each
    % triangle (stored in k)
    
    [k,vol] = convhull(ptcloud); % easy volume/triangulation of polyhedron
    
    
    % find area by calculating the area of each triangle
    n_triangle = size(k,1);
    areas = zeros(n_triangle,1);
    
    for i = 1:n_triangle
        pt1 = ptcloud(k(i,1),:);
        pt2 = ptcloud(k(i,2),:);
        pt3 = ptcloud(k(i,3),:);
        v1 = (pt2-pt1);
        v2 = (pt3-pt1);
        th = acos(dot(v2,v1)/(norm(v2)*norm(v1)));
        areas(i) = 1/2*norm(v2)*norm(v1)*sin(th);
    end
    
    area = real(sum(areas));
end