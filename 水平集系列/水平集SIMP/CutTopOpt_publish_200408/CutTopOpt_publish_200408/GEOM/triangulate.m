%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [triangles, cellcase]=triangulate(phi)
% Triangulate using the pattern shown below
% phi is the level set values,

% Description of the cell cases
% 0 denotes void x=0
% 1 denotes solid x=1
%
% 4--3 Node numbering
% |  |
% 1--2
%
% [8 4 2 1] * [n1 n2 n3 n4]' = cellcase
%
% +--+ 1--+ +--1 1--1 +--+ 1--+ +--1 1--1
% | 0| | 1| | 2| | 3| | 4| | 5| | 6| | 7|
% +--+ +--+ +--+ +--+ +--1 +--1 +--1 +--1
%
% +--+ 1--+ +--1 1--1 +--+ 1--+ +--1 1--1
% | 8| | 9| |10| |11| |12| |13| |14| |15|
% 1--+ 1--+ 1--+ 1--+ 1--1 1--1 1--1 1--1

cellcase=[8 4 2 1]*(0 < phi)';
% For each cell case create a corresponding number of lines
switch cellcase
    case {0, 15}
        triangles=[1 2 3; 1 3 4];
    case {1, 14}
        triangles=[1 2 6; 2 5 6; 2 3 5; 6 5 4];
    case {2, 13}
        triangles=[1 2 5; 1 5 6; 1 6 4; 5 3 6];
    case {3, 12}
        triangles=[1 2 5; 1 5 6; 6 5 3; 6 3 4];
    case {4, 11}
        triangles=[1 5 4; 5 6 4; 4 6 3; 5 2 6];
    case {5, 10}
        triangles=[1 5 8; 5 2 6; 6 3 7; 8 7 4; 5 6 7; 5 7 8;];
    case {6, 9}
        triangles=[1 5 6; 1 6 4; 5 2 3; 5 3 6;];
    case {7, 8}
        triangles=[1 5 6; 5 2 3; 5 3 6; 6 3 4;];
    otherwise
%         triangles = [];
        error('Case not valid');
end