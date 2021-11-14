function [IQs,Vs,As,Bounds,cell_d,angles,atomlocs,Flagged] = polymer_sphericity_calc(filename,phase)
    
    % This function reads a converged SCFT coordinate grid file of monomer
    % volume fraction, specified by "filename", and also takes the phase of
    % the calculation as an input.
    %
    % The function has the 
    % 
    % File #1: the results file, which gives the isoperimetric quotients, volumes, and areas of
    % each unique particle (.iq)
    %
    % File #2: the particle coordinates file, which gives the boundary
    % points of each particle analyzed in the phase (.bound)
    fprintf('Processing file of name: %s\n',filename)

    % read the converged SCFT solution coordinate grid data
    [R,x,y,z,dim,lattype,cell_d,angles,n_mnr,grid] = read_rgrid(filename); 
    phase = lower(phase);

    % unit cell properties
    %a = cell_d(1); b = cell_d(2); c = cell_d(3); alp = angle(1); bet = angle(2); gam = angle(3);
    %V_unitcell = a*b*c*sqrt(1+2*cos(alp)*cos(bet)*cos(gam)-cos(alp)^2-cos(bet)^2-cos(gam)^2);
    
    basis = get_basis(cell_d,angles);
    atomlocs = get_atomloc(phase); % miller indices locations of particles 
    atomlocs = round(((basis')*(atomlocs'))',7); % multiplied by basis vectors
    n_atoms = size(atomlocs,1);
    
    % initialize matrices for storing data
    IQs = zeros(n_atoms,1);
    Vs = zeros(size(IQs));
    As = zeros(size(IQs));
    Bounds = {};
    Ks = {};
    Radii = {};
    Thetas = {};
    Flags = zeros(n_atoms,1);
    
    % set the interpolation method
    
    t0 = tic;
    if (phase =="a15")||(phase=="sigma")
        r_type = 'spline';
    elseif (phase == "bcc")||(phase=="fcc")||(phase=="c15")
        r_type = 'spline';
    elseif (phase == "c14") || (phase=="z") || (phase=="c36")
        r_type = 'natural';
    else
        error('Phase input not expected.');
    end
    
    % for each "particle", find the particle's isoperimetric quotient, volume,
    % area, and boundary coordinates, among other things.
    
    for i = 1:n_atoms
        tn = tic;
        atomloc = atomlocs(i,:);
        [iq,v,a,bound,k,radii,thetas,flagged] = particle_sphericity(atomloc,R,x,y,z,cell_d,angles,r_type);
        IQs(i) = iq;
        Vs(i) = v;
        As(i) = a;
        Bounds{i} = round(bound,8);
        Ks{i} = k;
        Radii{i} = radii;
        Thetas{i} = thetas;
        Flags(i) = flagged;
        fprintf('%d/%d complete. %0.2f seconds.\n',[i,n_atoms,toc(tn)])
    end
    t = toc(t0);
    fprintf('Calculation took: %0.2f seconds\n',t)
    
    Flagged = any(Flags);
    
    % Code below can be used to recenter unit cell on a certain atom! (or
    % just on a certain location..)
    % for i = 1:n_atoms
    %    Bounds{i} = center_pbc(Bounds{i},atomlocs(1,:),cell_d,angles);
    % end
    
    % Code below can display the isosurface point clouds
    %disp_phase(Bounds,IQs,cell_d,angles)
    
    % clear unnecessary things for.... no reason i guess
    clear iq v a bound k radii thetas flagged    
    
end

function [sum,E] = rad_hist(radii,thetas,plotHist)
    % Calculates a histogram of "radii", or distances of points on the
    % isosurface from the particle center
    % Note that this histogram is reweighted to approximately account for 
    % the "area" represented by each point given that the angles are spaced 
    % evenly, as opposed to treating each point as representing the same
    % area of the particle's surface
    % N is number of bins
    
    h = histogram(radii);
    N = length(h.BinEdges)-1;
    [Y,E] = discretize(radii,N);

    areas = (radii.^2).*sin(thetas);
    sum = zeros(N,1);
    for i = 1:size(areas,1)
        binidx = Y(i);
        sum(binidx) = sum(binidx)+areas(i);
    end
    
    if plotHist
        bar(E(2:N+1),sum,'hist')
    end
    
end

function disp_phase(Bounds,IQs,cell_d,angles)
    % plots point clouds and colors them by isoperimetric quotient
    n = length(Bounds);
    IQmin = min(IQs);
    IQmax = max(IQs);
    
    figure()
    for i = 1:n
        bound = Bounds{i};
        if n~=1
            color = (IQs(i)-IQmin)/(IQmax-IQmin)*[205 0 -205] + [50 0 255];
            plot3(bound(:,1),bound(:,2),bound(:,3),'.','Color',color/255,'MarkerSize',15)
        else
            plot3(bound(:,1),bound(:,2),bound(:,3),'.','MarkerSize',15)
        end
        hold on
    end
    
    draw_lattice(cell_d,angles)
    hold off
    
    set(gca,'PlotBoxAspectRatio',[1,1,1],'DataAspectRatio',[1,1,1])

end

function draw_lattice(cell_d,angles)
    % Drawing the Unit Cell Outline

    axes = gca;
    x_pos(1) = 0;
    x_pos(2) = cos(angles(3))*cell_d(2);
    x_pos(3) = cell_d(1)+ cos(angles(3))*cell_d(2);
    x_pos(4) = cell_d(1);
    x_pos(5) = 0 + cos(angles(1))*cell_d(3);
    x_pos(6) = cos(angles(3))*cell_d(2) + cos(angles(1))*cell_d(3);
    x_pos(7) = cos(angles(3))*cell_d(2) + cos(angles(1))*cell_d(3) + cell_d(1);
    x_pos(8) = cell_d(1) + cos(angles(1))*cell_d(3);

    x_start = x_pos;
    for i =9:12
        x_start(i)=x_pos(i-8);
    end
    x_end = zeros(length(x_start),1)';

    for i =1:4
        if i == 4
            x_end(i)=x_start(1);
        else
            x_end(i)=x_start(i+1);
        end
    end


    for i =5:8
        if i == 8
            x_end(i)=x_start(5);
        else
            x_end(i)=x_start(i+1);
        end
    end

    for i =9:12
        x_end(i)=x_start(i-4);
    end

    X1 = [x_start;x_end];

    y_pos(1) = 0;
    y_pos(2) = sin(angles(3))*cell_d(2);
    y_pos(3) = sin(angles(3))*cell_d(2);
    y_pos(4) = 0;
    y_pos(5) = 0 + cos(angles(2))*cell_d(3);
    y_pos(6) = sin(angles(3))*cell_d(2) + cos(angles(2))*cell_d(3);
    y_pos(7) = sin(angles(3))*cell_d(2) + cos(angles(2))*cell_d(3);
    y_pos(8) = 0 + cos(angles(2))*cell_d(3);

    y_start = y_pos;
    for i =9:12
        y_start(i)=y_pos(i-8);
    end
    y_end = zeros(length(y_start),1)';

    for i =1:4
        if i == 4
            y_end(i)=y_start(1);
        else
            y_end(i)=y_start(i+1);
        end
    end

    for i =5:8
        if i == 8
            y_end(i)=y_start(5);
        else
            y_end(i)=y_start(i+1);
        end
    end

    for i =9:12
        y_end(i)=y_start(i-4);
    end

    Y1 = [y_start;y_end];

    z_pos(1) = 0;
    z_pos(2) = 0;
    z_pos(3) = 0;
    z_pos(4) = 0;
    z_pos(5) = cell_d(3)*sin(angles(1))*sin(angles(2));
    z_pos(6) = cell_d(3)*sin(angles(1))*sin(angles(2));
    z_pos(7) = cell_d(3)*sin(angles(1))*sin(angles(2));
    z_pos(8) = cell_d(3)*sin(angles(1))*sin(angles(2));

    z_start = z_pos;
    for i =9:12
        z_start(i)=z_pos(i-8);
    end
    z_end = zeros(length(z_start),1)';

    for i =1:4
        if i == 4
            z_end(i)=z_start(1);
        else
            z_end(i)=z_start(i+1);
        end
    end


    for i =5:8
        if i == 8
            z_end(i)=z_start(5);
        else
            z_end(i)=z_start(i+1);
        end
    end

    for i =9:12
        z_end(i)=z_start(i-4);
    end

    Z1 = [z_start;z_end];

    line(X1,Y1,Z1,'Color',[0,0,0],'LineStyle','-','LineWidth',1.2)

end

function [IQ,V,A,bound,k,radii,thetas,flagged] = particle_sphericity(atomloc,R,x,y,z,cell_d,angles,r_type)
    % for each particle, find its boundary, and use this to form a convex
    % hull of the particle and calculate the surface area and volume of the
    % particle
    
    [bound, radii, thetas] = find_boundary(atomloc,R,x,y,z,cell_d,angles,r_type);
    [bound, radii, flagged] = check_data(bound,atomloc,radii,cell_d,angles);
    [A,V,k] = get_polyhedron(center_pbc(bound,atomloc,cell_d,angles));
    
    IQ = calc_IQ(A,V);
end

function [bound, radii, flagged] = check_data(bound, atomloc, radii, cell_d, angles)
    
    flagged = false; % flag if there are outliers.
    window = round(length(radii)/4);
    
    outliers = isoutlier(radii,'movmean',window,'ThresholdFactor',4);
    if any(outliers)
        flagged = true; % flag this file for review
        
        if sum(outliers) < 5
            fprintf('Outliers (%d) found in current particle. Replacing.\n',sum(outliers))
            radii = filloutliers(radii,'pchip','movmean',window);
            
            % find the new vector from the center to the boundary point.
            % Must center coords on the particle to get actual vector
            boundTemp = center_pbc(bound,atomloc,cell_d,angles);
            atomlocTemp = center_pbc(atomloc,atomloc,cell_d,angles);
            newVectors = (boundTemp(outliers,:)-atomlocTemp)./vecnorm((boundTemp(outliers,:)-atomlocTemp),2,2).*radii(outliers);

            bound(outliers,:) = wrap_pbc(atomloc+newVectors,cell_d,angles);            
        else
            fprintf('Outliers (%d) found in current particle, but too many to fix. Will fix 5 max. \n',sum(outliers))
        end
    end
    
end

function [boundarycoords, radii, theta_out] = find_boundary(atomloc,R,x,y,z,cell_d,angles,r_type)
   
    
    % WHAT MONOMER IN CORE %
    % find which monomer is in the core of the particle, call it n_max.
    % isolate the volume fractions for that monomer
    % May need to modify this for generality (triblock??)

    cen_sub = find_subs(atomloc,x,y,z);
    r_center = R(cen_sub(1),cen_sub(2),cen_sub(3),:);
    [~,n_max]=max(r_center);
    
    R_max = R(:,:,:,n_max);
    
    % A SOLUTION FOR NON-RECTANGULAR LATTICES %
    % if not all right angles, create the scattered interpolating function! 
    % that way you aren't creating it over and over... 
    
    xdim = length(x(:,1,1));
    ydim = length(y(1,:,1));
    zdim = length(x(1,1,:));
    [xn,yn,zn] = ndgrid(linspace(0,x(xdim,ydim,zdim),xdim),linspace(0,y(xdim,ydim,zdim),ydim),linspace(0,z(xdim,ydim,zdim),zdim));
    
    if ~all(angles==pi/2)
        F = scatteredInterpolant(x(:), y(:), z(:), R_max(:),r_type,'linear');
    else
        F = griddedInterpolant(xn, yn, zn, R_max, r_type); % gridded interpolant time? might speed up! 
    end

    
    % THE MEAT - BIG "PARALLELIZED" FOR LOOP %
    
    % all of this is to make the loop suitable for parallelization through
    % parfor. make meshgrid of thetas and phis.
    
    % just kidding - no longer parallelized, because for some reason the
    % interpolant F has issues with parfor loops
    
    spacing = 0.02*pi;
    
    th_list = spacing:spacing:(pi-spacing);
    ph_list = 0:spacing:2*pi;
    l_th = length(th_list);
    l_ph = length(ph_list);
    [phis,thetas] = meshgrid(ph_list,th_list);
    
    boundarycoords = zeros(l_th*l_ph,3);
    radii = zeros(l_th*l_ph, 1);
    
    % loop over spherical angles to find the points at the
    % monomerinterface, identified by a given volume fraction (50%).This is
    % done using interpolation of the grid points.
    
    for i = 1:(l_th*l_ph)
        th = thetas(i);
        ph = phis(i);
        [bdpt,rbdpt] = find_bdpt(atomloc, th, ph, cell_d, angles, F);
        boundarycoords(i,:) = bdpt;
        radii(i) = rbdpt;
    end
    
    % now, do the edge cases, theta = 0 and theta = pi
    [bdpt0,rbdpt0] = find_bdpt(atomloc, 0, 0, cell_d, angles, F);
    [bdptpi,rbdptpi] = find_bdpt(atomloc, pi, 0, cell_d, angles, F);
    boundarycoords = [bdpt0; boundarycoords; bdptpi];
    radii = [rbdpt0; radii; rbdptpi];
    theta_out = [0; thetas(:); pi];
    
end

function [bdpt,r_bdpt] = find_bdpt(atomloc,theta,phi,cell_d,angles,F)
    % numerically solve for the point where the expansion is equal to 0.50,
    % IE the interface of the particle
    
    frac = 0.50; % the isosurface we want to find. you can change this!
    
    % decide which function to use, depending on whether a cubic grid or a
    % something else grid (triclinic? hexagonal?)
    
    fun = @(r) frac - ray_interp(atomloc, r, theta, phi, F, cell_d, angles);
    
    options = optimset('TolX',1e-6);
    
    % find approximate bounds for the function... hard-coded/engineered in
    r0 = 0.0001;
    r1 = 0.06*norm(cell_d);

    while sign(fun(r1))==sign(fun(r0))
        r0 = r1;
        r1 = r1 + 0.025*norm(cell_d);
        if r1 > cell_d*10^3
            error('Can not find optimal radius for fzero. Likely a broken simulation result, or something wrong with the algorithm.')
        end
    end

    % find the value of r, for this theta and phi, where function "fun" is zero.
    r_bdpt = fzero(fun,[r0,r1],options);

    % get coordinates of the boundary point, bdpt, and wrap everything into
    % the unit cell if it isn't already
    bdpt = atomloc+r_bdpt*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    bdpt = wrap_pbc(bdpt,cell_d,angles);
end

function f = ray_interp(atomloc, r, theta, phi, F, cell_d, angles)
    % evaluates field given a central point (atomloc) and (r,th,phi) for 
    % to describe the point, and other things about the lattice
    % find the point
    coord = atomloc+r*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    coord = wrap_pbc(coord,cell_d,angles); % keep ray within the unit cell.
    
    % evaluate the scattered or gridded interpolant 
    f = F(coord(1),coord(2),coord(3));

end

function coord = center_pbc(coord,atomloc,cell_d,angles)
    % center a set of coordinates on a point, atomloc, considering the
    % constraints of the unit cell.
    basis = get_basis(cell_d,angles);
    center = ((basis')*[1/2; 1/2; 1/2;])';
    for i = 1:size(coord,1)
        delta = coord(i,:) - atomloc;
        coord(i,:) = delta+center; % puts the atom at 1/2,1/2,1/2, not at 0,0,0
    end
    coord = wrap_pbc(coord,cell_d,angles); % bring everything back into the unit cell
end

function coord = wrap_pbc(coord,cell_d,angles)
    % find the position of a coordinate inside of the unit cell
    % coord - the coordinate (coordinates)
    % cell_d - unit cell lattice parameters
    % angle - unit cell angles (alp, bet, gam)
    %coord = round(coord,7);
    
    basis = get_basis(cell_d,angles);
    
    
    for i = 1:size(coord,1)
        coeff = (basis')\((coord(i,:))');
        for k = 1:3
            coeff(k)=coeff(k)-floor(coeff(k));
        end
        coord(i,:) = (basis')*coeff;
    end
    
end

function basis = get_basis(cell_d,angles)
    
    % construct basis vectors in cartesian coordinates based on lattice
    % parameters and angles. Generalized for triclinic.
    
    % by convention, align the first axis, with length a, with the x axis 
    
    a = cell_d(1);
    b = cell_d(2);
    c = cell_d(3);
    alpha = angles(1);
    beta = angles(2);
    gamma = angles(3);

    basis(1,1) = a;
    basis(2,1) = b * cos(gamma);
    basis(2,2) = b * sin(gamma);
    basis(3,1) = c * cos(beta);
    basis(3,2) = c * (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma);
    basis(3,3) = c * sqrt(1 - cos(beta)^2 - ((cos(alpha) - cos(beta)*cos(gamma))/sin(gamma))^2);
    basis = round(basis,8);
    
end

function subs = find_subs(loc,x,y,z)
    
    % find the subscripts of the closest grid point to "loc"
    pts = [x(:) y(:) z(:)];
    
    k = dsearchn(pts,loc); % uses qhull (such a useful package!) to find closest point. the way I have written this assumes only one is found..
    subs = zeros(size(loc));
    
    for i = 1:size(loc,1)
        
        [ix_,iy_,iz_] = ind2sub(size(x),k);
        
        subs(i,1) = ix_; 
        subs(i,2) = iy_; 
        subs(i,3) = iz_;
    end
end

function IQ = calc_IQ(A,V)
    % isoperimetric quotient. Equal to 1 for sphere. equal to 0 for thin
    % sheet
    IQ = 36*pi*V^2/A^3;
end
        
