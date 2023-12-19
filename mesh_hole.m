function [materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect...
    ,nelnodes,elident,nfix,fixnodes,ndload,dloads]=mesh_hole_part1(meshn)

%        Variables read from input file;
%        materialprops(i)    List of material parameters
%        ncoord              No. spatial coords (2 for 2D, 3 for 3D)
%        ndof                No. degrees of freedom per node (2 for 2D, 3 for 3D)
%                            (here ndof=ncoord, but the program allows them to be different
%                            to allow extension to plate & beam elements with C^1 continuity)
%        nnode               No. nodes
%        coords(i,j)         ith coord of jth node, for i=1..ncoord; j=1..nnode
%        nelem               No. elements
%        maxnodes            Max no. nodes on any one element (used for array dimensioning)
%        nelnodes(i)         No. nodes on the ith element
%        elident(i)          An integer identifier for the ith element.  Not used
%                            in this code but could be used to switch on reduced integration,
%                            etc.
%        connect(i,j)        List of nodes on the jth element
%        nfix                Total no. prescribed displacements
%        fixnodes(i,j)       List of prescribed displacements at nodes
%                            fixnodes(1,j) Node number
%                            fixnodes(2,j) Displacement component number (1, 2 or 3)
%                            fixnodes(3,j) Value of the displacement
%        ndload              Total no. element faces subjected to tractions
%        dloads(i,j)         List of element tractions
%                            dloads(1,j) Element number
%                            dloads(2,j) face number
%                            dloads(3,j), dloads(4,j)
%                            Components of traction in x and y direction

%lp: a temporary parameter to help creat mesh
lp=meshn+1;

E = 1000000;  % Young's modulus = 10e11 Pa
%nu = 0.4999;
nu=0.3; %Poission ratio

%shear modulus G /miu
materialprops(1)=E/2/(1+nu);

%Poisson ratio
materialprops(2)=nu;
materialprops(3)=1;

ncoord=2;
ndof=2;
nnode=(meshn+1)^2;
nelem=meshn^2;
maxnodes=4;
nelnodes=zeros(1,nelem);
nelnodes(1:nelem)=4;
elident=zeros(1,nelem);

%line 55 - 61: fixed nodes (displacement boundary condition) for Cook's membrane. 
nfix=lp*2;
fixnodes=zeros(3,2*lp);
for i = 1:lp
    fixnodes(1,i) = i;
end
for i = lp+1:nfix
    fixnodes(1,i) = (lp*(lp-2))+i;
end
fixnodes(2,1:lp) = 2;
fixnodes(2,lp+1:nfix) = 1;
fixnodes(3,:)=0;


%line 73 - 84: node coornidates for infinite plate with circular hole
lengthSquare = 0.004;
radiusHole = 0.001;

x=[]; y=[];
for i=1:lp

    innerStartx = radiusHole*cos(0.5*pi*(i-1)/(meshn));
    innerStarty = radiusHole*sin(0.5*pi*(i-1)/(meshn));
    if(i <= round(lp/2))
        outerStartx = lengthSquare;
        outerStarty = ((lengthSquare*2)/meshn)*(i-1);
    else
        outerStartx = (lengthSquare*2) - (((lengthSquare*2)/meshn)*(i-1));
        outerStarty = lengthSquare;
    end

    xp=linspace(innerStartx,outerStartx,lp);
    yp=linspace(innerStarty,outerStarty,lp);
    
    x=[x xp];
    y=[y yp];
end
coords=[x;y];


%line 64 - 69: traction (force boundary condition) for Cook's membrane. 
ndload=meshn;
dloads=zeros(4,meshn);


for i = 1:meshn
    nodeNumberLoads1 = (i*meshn)+i;
    nodeNumberLoads2 = (i*meshn)+i+(meshn+1);
    xcoordLoads = 0.5*(coords(1,nodeNumberLoads1)+coords(1,nodeNumberLoads2));
    ycoordLoads = 0.5*(coords(2,nodeNumberLoads1)+coords(2,nodeNumberLoads2));
    rLoads = sqrt(xcoordLoads^2 + ycoordLoads^2);
    thetaLoads = atan(ycoordLoads/xcoordLoads);

    stressxx = 1 - ((radiusHole^2)/(rLoads^2))*(1.5*cos(2*thetaLoads) + cos(4*thetaLoads)) + (1.5*((radiusHole^4)/(rLoads^4))*cos(4*thetaLoads));
    stressyy = -((radiusHole^2)/(rLoads^2))*(0.5*cos(2*thetaLoads) - cos(4*thetaLoads)) - (1.5*((radiusHole^4)/(rLoads^4))*cos(4*thetaLoads));
    stressxy = -((radiusHole^2)/(rLoads^2))*(0.5*sin(2*thetaLoads) + sin(4*thetaLoads)) + (1.5*((radiusHole^4)/(rLoads^4))*sin(4*thetaLoads));

    dloads(1,i) = i*meshn;
    dloads(2,i) = 2;
    if i <= round(meshn/2)
        dloads(3,i)=stressxx;
        dloads(4,i)=stressxy;
    else
        dloads(3,i)=stressxy;
        dloads(4,i)=stressyy;
    end
end

% line 86 to 97 connectivity for Cook's membrane
rowcount = 0;
connect=zeros(4,meshn^2);
for elementcount = 1:meshn^2
    connect(1,elementcount) = elementcount + rowcount;
    connect(2,elementcount) = elementcount + 1 + rowcount;
    connect(3,elementcount) = elementcount + (lp + 1) + rowcount;
    connect(4,elementcount) = elementcount + (lp) + rowcount;
    if mod(elementcount,lp-1) == 0
        rowcount = rowcount + 1;
    end
end

end



