%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mma = MMAInit(mesh,opt)

mma.m=1;
mma.n=mesh.nnodes;
mma.xold1=zeros(mma.n,1);
mma.xold2=zeros(mma.n,1);
mma.a0=1;
mma.a=ones(mma.m,1)*0;
mma.c=ones(mma.m,1)*100;
mma.d=ones(mma.m,1);
mma.df0dx2=zeros(mma.n,1);
mma.dfdx2=zeros(mma.n,1);
mma.U=zeros(mma.n,1);
mma.L=zeros(mma.n,1);
mma.asyminit=opt.asyminit;
mma.asymdecrease=opt.asymdecrease;
mma.asymincrease=opt.asymincrease;

