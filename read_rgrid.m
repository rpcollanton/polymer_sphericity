function [R,x,y,z,dim,lattype,cell_d,angles,n_mnr,grid] = read_rgrid(filename)

    % reads volume fraction coordinate grid data outputted by PSCF
    % outputs a 4-d array R, containing volume fraction values
    % the indices of R correspond to x point, y point, z point, and monomer type
    % x, y, and z contain the grid point coordinates
    % dim is 1, 2, or 3, denoting whether this is a 1D, 2D, or 3D simulation
    % lattype is the crystal system
    % cell_d contains the unit cell side lengths typically called a, b, and c
    % angles contains the three unit cell angles, alpha, gamma, and beta
    % n_mnr is the number of monomers
    % grid is the grid size (ie: [36, 36, 36] for 36 points in all directions)
    
    C = textread(filename, '%s','delimiter', '\n');

    ic = 0;
    start_row = 1000;
    
    % read in data about crystal from preamble.
    while ic <= start_row
        ic = ic +1;
        

        if strcmp(strrep(char(C(ic)),' ', ''),'dim')==1
            dim = str2double(C{ic+1});          % Reads the grid dimensions
        elseif strcmp(strrep(char(C(ic)),' ', ''),'crystal_system')==1
            lattype = strrep(C(ic+1), '''', '');   % Reads the system type
        elseif strcmp(strrep(char(C(ic)),' ', ''),'cell_param')==1
            param = sscanf(C{ic+1},'%f')';            % Reads the cell parameters
        elseif strcmp(strrep(char(C(ic)),' ', ''),'N_monomer')==1
            n_mnr = str2double(C{ic+1});        % Reads the number of monomers
        elseif strcmp(strrep(char(C(ic)),' ', ''),'ngrid')==1
            grid = sscanf(C{ic+1},'%f')';            % Reads the grid size
            start_row = ic+2;       % record the row in which the volume fractions start
        end
    end

    % record the row in which the supplementary info ends
    end_info = start_row - 1;

    % Reading the values from the file, but not yet organized in grid form..
    A = zeros(length(C) - end_info,n_mnr);
    for i = start_row:length(C)
        A(i - end_info,:) = sscanf(C{i},'%f')';
    end

    % Get x,y,z grid points and unit cell parameters
    [x,y,z,cell_d,angles] = gen_xyz(lattype,param,grid);

    % place points from A on grid, for each monomer. R is a 4D array.
    % 3 dims are for the x,y,z grid. The 4th dim is for the monomer type.
    R = rearrangePts(A, grid, dim);

end

function R = rearrangePts(A, grid, dim)

    n_mnr = size(A,2);
    R = zeros([grid n_mnr]); % store values of volume fraction on grid
    counter = 0;
    
    for iz=1:grid(3)+1
        for iy=1:grid(2)+1
            for ix=1:grid(1)+1
                counter = counter + 1;
                for in = 1:n_mnr
                    if ix == grid(1)+1
                        R(grid(1)+1,:,:,in) = R(1,:,:,in); % periodic bc's
                        counter = counter - (1/n_mnr); 
                    elseif iy == grid(2) + 1
                        R(:,grid(2)+1,:,in) = R(:,1,:,in);
                        counter = counter - (1/n_mnr);
                    elseif iz == grid(3) + 1
                        R(:,:,grid(3)+1,in) = R(:,:,1,in);
                        counter = counter - (1/n_mnr);
                    else
                        R(ix,iy,iz,in) = A(round(counter),in);
                    end
                end
            end
            if(dim==1)
                counter = 0;
            end
        end
        if(dim==2)
            counter=0;
        end
    end
end

function [x,y,z,cell_d,angle] = gen_xyz(lattype, param, grid)

    % get crystal system dimensions and angles
    if strcmp(lattype,'hexagonal') == 1
        angle = [pi/2 pi/2 (2*pi)/3];
        cell_d = param;
    elseif strcmp(lattype,'cubic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'tetragonal') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'orthorhombic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'triclinic') == 1
        angle = [param(4) param(5) param(6)];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'monoclinic') == 1
        angle = [pi/2 param(4) pi/2];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'trigonal') == 1
        angle = [param(2) param(2) param(2)];
        cell_d = [param(1)];
    elseif strcmp(lattype,'lamellar') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    else
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    end
    
    % extract unit cell dimensions, such that there are 3 for each system.
    % PSCF will output only 1 if they are all the same, for example. 
    if(length(cell_d)==1) % Cubic crystals
        new_cell(1:3) = cell_d;             
    elseif(length(cell_d)==2) % Tetragonal crystals
        new_cell(1:2) = cell_d(1);            
        new_cell(3)   = cell_d(2);
    else % Orthorhombic crystals
        new_cell = cell_d;                    
    end
    clear cell_d; cell_d = new_cell; % replace cell_d
    
    % for 2D and 1D crystals, "fill in" the grid.
    if(length(grid)==1) % 3D grid for 1D crystals
        grid(2) = grid(1); 
        grid(3) = grid(1);
    elseif(length(grid)==2) % 3D grid for 2D crystals
        grid(3) = grid(1); 
    end
    
    % matrices for grid coords
    x = zeros(grid); % x coords on grid
    y = zeros(grid); % y coords on grid
    z = zeros(grid); % z coords on grid
    
    nround = 10;
    % calculate the coordinates in x,y,z.
    for iz=1:grid(3)+1
        for iy=1:grid(2)+1
            for ix=1:grid(1)+1
                % calulate coordinate values using angles and cell parameters
                xtemp = cell_d(1) * (ix-1)/grid(1) + (cos(angle(3)))*( cell_d(2) * (iy-1)/grid(2)) + ((iz-1)/grid(3))*(cos(angle(1))*cell_d(3));;
                ytemp = cell_d(2) * (iy-1)/grid(2) * sin(angle(3)) + ((iz-1)/grid(3))*cos(angle(2))*cell_d(3);
                ztemp = cell_d(3) * (iz-1)/grid(3) * sin(angle(1)) * sin(angle(2));
                % round to rectify numerical nonsense.
                x(ix,iy,iz) = round(xtemp,nround);
                y(ix,iy,iz) = round(ytemp,nround);
                z(ix,iy,iz) = round(ztemp,nround);
            end
        end
    end
end