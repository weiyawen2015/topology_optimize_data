function Top3dSTL_v3(fout, varargin)
%Top3dSTL_v3   A simple STL writter for Top3d by Liu (Apr 2015)
%   Top3dSTL_v3(fout) writes a STL file with name fout using cubic
%   representation and binary file format if xPhys exists in Workspace.
%
%   Top3dSTL_v3(fout, xPhys) writes a STL file with name fout using cubic
%   representation and binary file format
%
%   Top3dSTL_v3(___, Name, Value) writes a STL file with one or more Name,
%   Value pair argments. Use this option with any of the input argument
%   combinations in the prvious syntaxes.
%       FORMAT     - File is written in 'binary' (default) or 'ascii' format.
%       TITLE      - Header text (max 80 characters) written to the STL file.
%       MODE       - Facets are generated using 'cube' (default) or 'iso'.
%       CUTOFF     - Density cutoff value. default: 0.5
%       FACECOLOR  - Face color. default: 'cyan'
%       ALPHA      - Face alpha value. default: 1.
%       UNITLENGTH - Vector of element unit length. default: [1 1 1]
%       PLOT       - Logic flag to display structures. default: true
%
%   Example 1:
%       Top3dSTL_v3('MyTop3d.stl') % when xPhys is in Workspace
%
%   Example 2:
%       Top3dSTL_v3('MyTop3d.stl', density, ...
%       'Format', 'ascii', 'Mode', 'iso', 'FaceColor', 'm', 'Plot', false)
%

% Determine input type
if ~isempty(varargin) && (isnumeric(varargin{1}) || islogical(varargin{1}))
    xPhys = varargin{1};
    options = parseInputs(varargin{2:end});
else
    try
        xPhys = evalin('base', 'xPhys');
    catch ME
        switch ME.identifier
            case 'MATLAB:UndefinedFunction'
                error('xPhys is not input argument nor exist in workspace');
            otherwise
                rethrow(ME)
        end
    end
    options = parseInputs(varargin{:});
end

% Generate faces and verts
if strcmp(options.mode, 'cube')
    [faces, verts] = getCube(xPhys, options);
else
    [faces, verts] = getISO(xPhys, options);
end
% Facets
facets = single(verts);
facets = reshape(facets(faces',:)', 3, 3, []);
% facets: (:,:,1) --> Vertices of face 1,
% facets(:,1,1)   --> First vertice of face 1
V1 = squeeze(facets(:,2,:) - facets(:,1,:));
V2 = squeeze(facets(:,3,:) - facets(:,1,:));

% Normal vectors
normals = cross(V1, V2);
clear V1 V2
% Normal vectors normalization
normals = bsxfun(@times, normals, 1 ./ sqrt(sum(normals .* normals, 1)));

facets = cat(2, reshape(normals, 3, 1, []), facets);
clear normals

% Write STL
if strcmp(options.format, 'ascii')
    writeAscii(facets, fout, options.title);
else
    writeBinary(facets, fout, options.title);
end

end

function options = parseInputs(varargin)
OP = inputParser;
defaultFormat    = 'binary';
expectedFormat   = {'ascii', 'binary'};
defaultMode      = 'cube';
expectedMode     = {'cube', 'iso'};
defaultTitle     = sprintf('Created by Top3dSTL.m %s',datestr(now));
defaultCutoff    = 0.5;
defaultFcolor    = 'c';
defaultAlpha     = 1;
defaultUnitLegth = [1, 1, 1];
defaultPlot      = true;

OP.addParamValue('format', defaultFormat, ...
    @(x) any(validatestring(x,expectedFormat)))
OP.addParamValue('mode', defaultMode, ...
    @(x) any(validatestring(x,expectedMode)))
OP.addParamValue('title', defaultTitle, @ischar);
OP.addParamValue('cutoff', defaultCutoff, @isnumeric)
OP.addParamValue('facecolor', defaultFcolor, @ischar)
OP.addParamValue('alpha', defaultAlpha, @isnumeric)
OP.addParamValue('unitlength', defaultUnitLegth, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
OP.addParamValue('plot', defaultPlot, @islogical)

OP.parse(varargin{:});
options = OP.Results;
end

function [faces, verts] = getCube(xPhys, options)
% Generate Mesh
nx = options.unitlength(1);
ny = options.unitlength(2);
nz = options.unitlength(3);
[nely, nelx, nelz] = size(xPhys);
nele    = nelx*nely*nelz;
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids, size(nodeidz))+repmat(nodeidz, size(nodeids));
enodVec = nodeids(:)+1;
enodMat = repmat(enodVec,1,8)+ ...
    repmat([0 nely + [1 0] -1 ...
    (nely+1)*(nelx+1)+[0 nely + [1 0] -1]],nele,1);

% Faces connectivities
enodidx = [...
    1 3 2; 1 4 3; ... % back
    5 6 7; 5 7 8; ... % front
    1 5 8; 1 8 4; ... % left
    6 2 3; 6 3 7; ... % right
    8 7 3; 8 3 4; ... % up
    1 2 6; 1 6 5];    % down

faces = [];

% Filter out Low density
xPhys(xPhys < options.cutoff)  = 0;
xPhys(xPhys >= options.cutoff) = 1;

for f = 1:size(enodMat,1)
    if xPhys(f) == 0
        continue;
    end
    
    % Coordinates
    [j, i, k] = ind2sub([nely, nelx, nelz], f);
    eFace   = enodidx;      % element faces connectivities
    idx     = [];           % element faces to be deleted
    % Neighbor on back
    if (k ~= 1 && xPhys(j, i, k - 1) == 1)
        idx = [idx; [1 2]];
    end
    % Neighbor on front
    if (k ~= nelz && xPhys(j, i, k + 1) == 1)
        idx = [idx; [3 4]];
    end
    % Neighbor on left
    if (i ~= 1 && xPhys(j, i - 1, k) == 1)
        idx = [idx; [5 6]];
    end
    % Neighbor on right
    if (i ~= nelx && xPhys(j, i + 1, k) == 1)
        idx = [idx; [7 8]];
    end
    % Neighbor on up
    if (j ~= 1 && xPhys(j - 1, i, k) == 1)
        idx = [idx; [9 10]];
    end
    % Neighbor on down
    if (j ~= nely && xPhys(j + 1, i, k) == 1)
        idx = [idx; [11 12]];
    end
    
    eFace(idx, :) = [];
    tmp   = enodMat(f,:);
    faces = cat(1, faces, tmp(eFace));

end

% Vertices
[xx, yy, zz] = meshgrid(0:nx:nelx*nx, ...
    0:ny:nely*ny, ...
    0:nz:nelz*nz);
verts = [xx(:) nely-yy(:) zz(:)];

% Visualization
if options.plot
    dverts = [xx(:) zz(:) yy(:)];
    cla, hold on, view(30,30), rotate3d on, axis equal
    axis([0 nelx*nx 0 nelz*nz 0 nely*ny]), box
    set(gca,'YDir', 'reverse', 'ZDir', 'reverse', 'ZtickLabel', flipud(get(gca, 'Ztick')'));
    patch('faces', faces, 'vertices', dverts, 'FaceColor', options.facecolor, 'FaceAlpha', options.alpha)
    xlabel('x'), ylabel('z'), zlabel('y')
end
end

function [faces, verts] = getISO(xPhys, options)
nx = options.unitlength(1);
ny = options.unitlength(2);
nz = options.unitlength(3);
[nely, nelx, nelz] = size(xPhys);

aux = zeros(nely+2, nelx+2, nelz+2);
aux(2:end-1,2:end-1,2:end-1) = xPhys;

[X,Y,Z] = meshgrid(0:nx:nx*(nelx+1), ...
    0:ny:ny*(nely+1), ...
    0:nz:nz*(nelz+1));

[faces, verts] = isosurface(X-0.5, Z-0.5, Y-0.5, aux, options.cutoff);

% Visualization
if options.plot
    cla, hold on, view(30,30), rotate3d on, axis equal
    axis([0 nx*nelx 0 nz*nelz 0 ny*nely]), box
    set(gca, 'YDir', 'reverse', 'ZDir', 'reverse', 'ZtickLabel', flipud(get(gca, 'Ztick')'));
    
    patch('Faces', faces, 'Vertices', verts,...
        'FaceColor', options.facecolor, 'EdgeColor', 'none', 'FaceAlpha', options.alpha);
    camlight, lighting gouraud;
    xlabel('x'), ylabel('z'), zlabel('y')
    drawnow
end
end

function writeAscii(facets, fout, title)
% Write ASCII STL file
%{
FORMAT:

solid name
    facet normal ni nj nk
        outer loop
            vertex v1x v1y v1z
            vertex v2x v2y v2z
            vertex v3x v3y v3z
        endloop
    endfacet
end solid name

%}
fid = fopen(fout, 'wb+');
fprintf(fid, [title, '\r\n']);
fprintf(fid,[...
    'facet normal %.7E %.7E %.7E\r\n' ...
    'outer loop\r\n' ...
    'vertex %.7E %.7E %.7E\r\n' ...
    'vertex %.7E %.7E %.7E\r\n' ...
    'vertex %.7E %.7E %.7E\r\n' ...
    'endloop\r\n' ...
    'endfacet\r\n'], facets);
fprintf(fid, ['end ', title, '\r\n']);
fclose(fid);
fprintf('Wrote %d facets to %s\n',size(facets, 3), fout);
end

function writeBinary(facets, fout, title)
% Write Binary STL file
%{
FORMAT:

UINT8[80] ? Header
UINT32 ? Number of triangles

foreach triangle
    REAL32[3] ? Normal vector
    REAL32[3] ? Vertex 1
    REAL32[3] ? Vertex 2
    REAL32[3] ? Vertex 3
    UINT16 ? Attribute byte count
end

%}
fid = fopen(fout, 'wb+');
fprintf(fid, '%-80s', title);                  % Title
fwrite(fid, size(facets, 3), 'uint32');        % Number of facets
facets = typecast(facets(:), 'uint16');        % Convert to unit16
facets = reshape(facets, 12*2, []);
facets(end+1, :) = 0;                          % Add color(0) to the end of each facet
fwrite(fid, facets, 'uint16');
fclose(fid);
fprintf('Wrote %d facets to %s\n',size(facets, 2), fout);
end

% =========================================================================
% === This code was written by K Liu, Dept. of Mechanical Engineering   ===
% === Indiana University-Purdue University Indianapolis,                ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === This code is an extension for Top3D program that helps users to   ===
% === generate STL file for 3D printing or post processing              ===
% === ----------------------------------------------------------------- ===
% === This code as well as Top3d program can be downloaded freely from  ===
% === the website: http://www.top3dapp.com/                             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, and =
% === they shall not be liable in any event caused by the use of the code.=
% =========================================================================