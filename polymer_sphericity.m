function polymer_sphericity(filelist_txt)
    % "filelist" is a text file containing a list of file paths. They are 
    % assumed to have a directory structure that includes the following
    % folders in no particular order:
    %
    % .../eps_Y/...
    % .../chiN_X/...
    % .../PHASE/...
    %
    % Y gives the value of epsilon, the conformational asymmetry
    % PHASE gives the type of phase being analyzed - A15,BCC,etc
    % X gives the value of chiN
    %
    % The file type is "rho_rgrid", which is a volume fraction field in
    % coordinate grid format storing the results of a converged SCFT
    % simulation done with Polymer Self-Consistent Field (PSCF).
    % 
    % It is required that the accompanying "out" file is in the same folder 
    % as the "rho_rgrid" file
    
    % BEGIN CODE %
    
    % open list of file paths.
    % this list can be made using the unix "find" command. 
    
    filelistID = fopen(filelist_txt, 'r');
    filelist = strsplit(fscanf(filelistID,'%c'),'\n');
    
    % number of listed files
    nFiles = length(filelist);
    % array to store paths of files with "failed" analyses, where an error
    % was thrown
    failed = {};
    
    % loop over all rho_rgrid files
    for f = 1:nFiles
        clear IQs Vs As Bounds cell_d angles atomlocs Flagged
		clear phase eps chiN grid blockL
		
		filename = filelist{f};
        % Determine the phase, epsilon value, and chiN value from the path.
        % NOTE %
        % we could incorporate this into the "findParams" function somehow. 
        % This may be preferable long term for robustness against different directory structures than the one in this dataset
        % We can get eps (ratio of kuhn length of monomer 1 to 2) and chiN (directly) from the .out file. Phase isn't in there, only space group, so that HAS to be in the directory somewhere!!! 

		try
	        [phase,eps,chiN] = parsePath(filename);
    	    % skip 'p6mm' (hex) because it is not a particle-forming phase.
        	if strcmp(phase,'p6mm')
            	continue
        	end
        	% skip files that have already been processed
        	if isfile(strrep(filename,'rho_rgrid','iq'))
	            fprintf("File already processed: %s\n",filename)
    	        continue
        	end
        	% From the "out" file, get some relevant parameters like block
        	% length, grid sizes. Could expand easily to get more if we wanted.
        	[grid,blockL] = findParams(filename);
        catch EX
            % if extracting result details fails, mark this file as
            % "failed" and move onto the next file
            failed{end+1} = strcat(filename,'\n');
			fprintf("Prepping of file failed: %s\n",filename)
        	warning(EX.message)
			continue
		end

		% Perform calculation to find surface cloud and isoperimetric quotient. 
        % Involves identifying atom centers, and then finding the 50% volume
        % fraction isosurface around that center.
        try
            [IQs,Vs,As,Bounds,cell_d,angles,atomlocs,Flagged] = polymer_sphericity_calc(filename,phase);
        catch EX
            % if geometric fails for any reason, mark as failed and
            % continue to the next file.
            failed{end+1} = strcat(filename,'\n');
            fprintf("Analysis of file failed: %s\n",filename)
            warning(EX.message)
			continue
        end
        
        % Write calculated data to two files: one IQ file that stores results, one bounds file that stores the boundary points of each sphere
        writeData(filename,IQs,Vs,As,Bounds,atomlocs,phase,eps,chiN,cell_d,angles,grid,blockL,Flagged)
        
    end
	writeToFile(strcat(filelist_txt,'_failed'),failed)
end

function [gridsize,blocklength] = findParams(filepath)
    
    % get path of the .out file. it is in same path as rho_rgrid, by my
    % personal design.
    
    outpath = strrep(filepath,'rho_rgrid','out');
    
    % get ID of file
    outID = fopen(outpath,'r');

    % read the file
    out = strsplit(fscanf(outID,'%c'),'\n');
    
    % declare variables
    blocklength = [];
    gridsize = []; 
    
    % loop through file to find variable values
    for i = 1:length(out)
        if ~isempty(blocklength) && ~isempty(gridsize)
            break
        end
        
        line = strip(out{i});
        if strcmp(line,'block_length')
            blockStr = strip(out{i+1});
            blocklength = str2double(strsplit(blockStr));
            
        elseif strcmp(line,'ngrid')
            gridStr = strip(out{i+1});
            gridsize = str2double(strsplit(gridStr));
            
        end
    end
end

function [phase,eps,chiN] = parsePath(filepath)
    
    % Split string by '/'
    folders = strsplit(filepath,'/');
    
    phases = {'a15','c15','c14','z','sigma','fcc','bcc','p6mm'};
    
    % Declare variables
    chiN = [];
    eps = [];
    phase = [];
    
    % Isolate eps, phase, and chiN:
    for i = 1:length(folders)
        
        % if we already got the values of chiN, epsilon, and phase, stop!
        if ~isempty(chiN) && ~isempty(eps) && ~isempty(phase)
            break
        end
        
        % check for information
        dirname = folders{i};
        
        if contains(dirname,'chiN_')
            chiNstr = strrep(dirname,'chiN_',''); % isolate value as string
            chiN = str2double(strrep(chiNstr,'_','.')); % replace underscore with decimal point, convert to double
        elseif contains(dirname,'eps_')
            epsStr = strrep(dirname,'eps_','');
            eps = str2double(strrep(epsStr,'_','.'));
        else
            % check to see if it is one of the phases
            for j = 1:8
                if strcmp(phases{j},lower(dirname))
                    phase = lower(dirname);
                end
            end
        end            
    end 
end

function writeData(filename,IQs,Vs,As,Bounds,atomlocs,phase,eps,chiN,cell_d,angles,gridsize,blocklength,Flagged)
    
    header = get_header(phase,eps,chiN,cell_d,angles,gridsize,blocklength,Flagged); % in cell format
    nAtoms = length(IQs);
    
    iqPath = strrep(filename,'rho_rgrid','iq');
    
    % Process IQ values and save to a file. %
    iqCell = {};
    iqCell{end+1} = 'filetype\n';
    iqCell{end+1} = 'iq\n\n';
    iqCell = [iqCell header];
    
    iqCell{end+1} = '[atomlocs]\n';
    for i = 1:nAtoms
        iqCell{end+1} = strcat(num2str(atomlocs(i,:),'%10.6f'),'\n');
    end
    iqCell{end+1} = '\n';
    
    iqCell{end+1} = '[isoperimetric_quotients]\n';
    for i = 1:nAtoms
        iqCell{end+1} = strcat(num2str(IQs(i)),'\n');
    end
    iqCell{end+1} = '\n';
    
    iqCell{end+1} = '[volumes]\n';
    for i = 1:nAtoms
        iqCell{end+1} = strcat(num2str(Vs(i)),'\n');
    end
    iqCell{end+1} = '\n';
    
    iqCell{end+1} = '[areas]\n';
    for i = 1:nAtoms
        iqCell{end+1} = strcat(num2str(As(i)),'\n');
    end
    iqCell{end+1} = '\n';
    iqCell{end+1} = 'FINISH';

    % Now process the boundary points and save them to a file. % 
    boundPath = strrep(filename,'rho_rgrid','bound');
    
    boundCell = {};
    boundCell{end+1} = 'filetype\n';
    boundCell{end+1} = 'bound\n\n';
    boundCell = [boundCell header];
    
    for i = 1:nAtoms
        boundCell{end+1} = '[atomnum]\n';
        boundCell{end+1} = strcat(num2str(i),'\n');
        boundCell{end+1} = '[atomloc]\n';
        boundCell{end+1} = strcat(num2str(atomlocs(i,:),'%10.6f'),'\n');
        boundCell{end+1} = '[boundarypoints]\n';
        Bounds_i = Bounds{i};
        nPoints = size(Bounds_i,1);
        for j = 1:nPoints
            boundCell{end+1} = strcat(num2str(Bounds_i(j,:),'%10.6f'),'\n');
        end
        boundCell{end+1} = '\n';
    end
    boundCell{end+1} = '\n';
    boundCell{end+1} = 'FINISH';
    
    fprintf("Writing IQ file: %s\n",iqPath)
    writeToFile(iqPath,iqCell);
    fprintf("Writing bounds file: %s\n",boundPath)
    writeToFile(boundPath,boundCell);

end

function header = get_header(phase,eps,chiN,cell_d,angles,gridsize,blocklength,Flagged) 

    header = {};
        
    header{end+1} = '---------------------------------------------------\n\n';
    
    header{end+1} = 'eps\n';
    header{end+1} = strcat(num2str(eps),'\n\n');
    
    header{end+1} = 'chiN\n';
    header{end+1} = strcat(num2str(chiN),'\n\n');
    
    header{end+1} = 'phase\n';
    header{end+1} = strcat(phase,'\n\n');
    
    header{end+1} = 'block_length\n';
    header{end+1} = strcat(num2str(blocklength),'\n\n');
    
    header{end+1} = 'cell_d\n';
    header{end+1} = strcat(num2str(cell_d),'\n\n');
    
    header{end+1} = 'angles\n';
    header{end+1} = strcat(num2str(angles),'\n\n');
    
    header{end+1} = 'grid\n';
    header{end+1} = strcat(num2str(gridsize),'\n\n');
    
    header{end+1} = 'flagged_for_review\n';
    header{end+1} = strcat(mat2str(Flagged),'\n\n');

    header{end+1} = '---------------------------------------------------\n\n'; 

end

function writeToFile(abspath,cellarray)
    % accepts a 1xn or nx1 cell array, writes each cell as a line in a file
    
    fileID = fopen(abspath,'w');
    nLines = max(size(cellarray));
    for i = 1:nLines
        fprintf(fileID,cellarray{i});
    end
    fclose(fileID);
end
