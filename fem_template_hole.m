function fem_template_hole

% ==================== Read data from the input file ===========================
%
close all;
clear all;
meshn = 20;

[materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect...
    ,nelnodes,elident,nfix,fixnodes,ndload,dloads]=mesh_hole(meshn);


% Plot the initial mesh as a check
figure
plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g');
hold on
for oi = 1:size(coords,2)
    plot(coords(1,oi),coords(2,oi),'k.',"MarkerSize",20);
    text(coords(1,oi)+0.5,coords(2,oi)+0.5,sprintf('%d',oi),'Fontsize',15);
end
title('initial mesh and node number');
%
%============================ MAIN FEM ANALYSIS PROCEDURE ========================
%
%   dofs        Nodal displacements.  Let u_i^a be ith displacement component
%               at jth node.  Then dofs contain (u_1^1, u_2^1, u_1^2, u_2^2....) for 2D
%               and (u_1^1, u_2^1, u_3^1, u_1^2, u_2^2, u_3^2....) for 3D
%
%   K           Global stiffness matrix.  Stored as (K_1111 K_1112  K_1121  K_1122...
%                                                    K_1211 K_1212  K_1221  K_1222...
%                                                    K_2111 K_2112  K_2121  K_2122...)
%               for 2D problem and similarly for 3D problem
%   r           Force vector.  Currently only includes contribution from tractions
%               acting on element faces (i.e. body forces are neglected)
%
dofs = zeros(ndof*nnode,1);

K = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs);

r = globaltraction(ncoord,ndof,nnode,ndload,coords,nelnodes,elident,connect,dloads,dofs);

%
%  Fix constrained nodes.
%
for n = 1:nfix
    rw = ndof*(fixnodes(1,n)-1) + fixnodes(2,n);
    for cl = 1:ndof*nnode
        K(rw,cl) = 0;
    end
    K(rw,rw) = 1.;
    r(rw) = fixnodes(3,n);
end
%
% Solve for the displacements
%

dofs = K\r;

%
%================================= POST-PROCESSING =================================
%
% Create a plot of the deformed mesh
%

defcoords = zeros(ndof,nnode);
scalefactor = 1.0;
for i = 1:nnode
    for j = 1:ndof
        defcoords(j,i) = coords(j,i) + scalefactor*dofs(ndof*(i-1)+j);
    end
end

figure
hold on
plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g');
plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'r');
title('deformed and undeformed mesh');
hold off

%
%================================= POST-PROCESSING FOR STRESS =================================
%
sss=numberofintegrationpoints(ncoord,nelnodes,elident);
sigma=zeros(ncoord,ncoord,sss,nelem);

lmncoord = zeros(ncoord,maxnodes);
displacements = zeros(ndof,maxnodes);

%
%   Loop over all the elements
%
for lmn = 1:nelem


    %
    %   Extract coords of nodes, DOF for the current element
    %
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
            displacements(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
    end
    n = nelnodes(lmn);
    ident = elident(lmn);

    npoints = numberofintegrationpoints(ncoord,n);
    dNdx = zeros(n,ncoord);
    dxdxi = zeros(ncoord,ncoord);
    strain = zeros(ndof,ncoord);
    xi = zeros(ncoord,1);
    x = zeros(ncoord,1);
    %
    %  Set up integration points
    %
    xilist = integrationpoints(ncoord,n,npoints);
    %
    %  Loop over the integration points
    %
    for intpt = 1:npoints

        %     Compute shape functions && derivatives wrt local coords
        %
        for i = 1:ncoord
            xi(i) = xilist(i,intpt);
        end
        N = shapefunctions(n,ncoord,ident,xi);
        dNdxi = shapefunctionderivs(n,ncoord,ident,xi);
        %
        %     Compute the coords of the integration point
        %
        for i = 1:ncoord
            x(i) = 0.;
            for a = 1:n
                x(i) = x(i) + lmncoord(i,a)*N(a);
            end
        end
        %
        %     Compute the jacobian matrix && its determinant
        %
        for i = 1:ncoord
            for j = 1:ncoord
                dxdxi(i,j) = 0.;
                for a = 1:n
                    dxdxi(i,j) = dxdxi(i,j) + lmncoord(i,a)*dNdxi(a,j);
                end
            end
        end

        dxidx = inv(dxdxi);
        %
        %     Convert shape function derivatives:derivatives wrt global coords
        %
        for a = 1:n
            for i = 1:ncoord
                dNdx(a,i) = 0.;
                for j = 1:ncoord
                    dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
                end
            end
        end
        %
        %     Compute the (infinitesimal) strain by differentiating
        %     displacements, similar as Exercise 8 on Page 107, Hughes book
        %
        %------------IMPORTANT------------------------
        %For Bbar method: need to change dNdx here to B_bar

        for i = 1:ncoord
            for j = 1:ncoord
                strain(i,j) = 0.;
                for a = 1:n
                    strain(i,j) = strain(i,j) + 0.5*(displacements(i,a)*dNdx(a,j)+displacements(j,a)*dNdx(a,i));
                end
            end
        end

        lmnstress = materialstress(ndof,ncoord,strain,materialprops);
        sigma(:,:,intpt,lmn) = lmnstress;

    end


end

%L2 projection of stress post-processing
[L2_LHS, L2_RHS]=L2_global(ncoord,ndof,nnode,coords, ...
    nelem,maxnodes,elident,nelnodes,connect,squeeze(sigma(1,1,:,:))  );  %sigma_22
%     nelem,maxnodes,elident,nelnodes,connect,squeeze(-1/ncoord*(sigma(1,1,:,:)+sigma(2,2,:,:)) ) );  %hydrostatic pressure


stressplot=L2_LHS \ L2_RHS;

%displaying corner stress concentration
disp(squeeze(sigma(1,1,:,(meshn*(meshn-1))+1)));

%L2 energy norm
radiusHole = 0.001;
C = materialstiffness(ndof,ncoord,0,materialprops); %calculating material stiffness matrix
D = zeros(3,3);
D(1,1) = C(1,1,1,1); D(2,1) = C(2,2,1,1); D(1,2) = C(1,1,2,2); D(2,2) = C(2,2,2,2); D(3,3) = C(1,2,1,2);

stressElementsFEM = zeros(3,nelem);
for i = 1:nelem
    stressElementsFEM(1,i) = sigma(1,1,1,i);
    stressElementsFEM(2,i) = sigma(2,2,1,i);
    stressElementsFEM(3,i) = sigma(1,2,1,i);
end

strainElementsFEM = zeros(3,nelem);
for i = 1:nelem
    strainElementsFEM(:,i) = D\stressElementsFEM(:,i);
end

stressElementsExact = zeros(3,nelem);
for i = 1:nelem
    nodeofElement = connect(1,i);
    xcoordNode = coords(1,nodeofElement);
    ycoordNode = coords(2,nodeofElement);
    r = sqrt(xcoordNode^2 + ycoordNode^2);
    theta = atan(ycoordNode/xcoordNode);

    stressxx = 1 - ((radiusHole^2)/(r^2))*(1.5*cos(2*theta) + cos(4*theta)) + (1.5*((radiusHole^4)/(r^4))*cos(4*theta));
    stressyy = -((radiusHole^2)/(r^2))*(0.5*cos(2*theta) - cos(4*theta)) - (1.5*((radiusHole^4)/(r^4))*cos(4*theta));
    stressxy = -((radiusHole^2)/(r^2))*(0.5*sin(2*theta) + sin(4*theta)) + (1.5*((radiusHole^4)/(r^4))*sin(4*theta));
    stressElementsExact(1,i) = stressxx;
    stressElementsExact(2,i) = stressyy;
    stressElementsExact(3,i) = stressxy;
end

strainElementsExact = zeros(3,nelem);
for i = 1:nelem
    strainElementsExact(:,i) = D\stressElementsExact(:,i);
end

for i = nelem
    l2norm = sqrt(0.5*(stressElementsFEM(1,i)-stressElementsExact(1,i))*(strainElementsFEM(1,i)-strainElementsExact(1,i)));
end
disp(l2norm)


figure
for ooe=1:nelem

    XX=zeros(1,5); YY=zeros(1,5); dd=zeros(1,5);

    for ooq=1:5
        tempooq=mod(ooq,4);
        if tempooq==0
            tempooq=4;
        end
        XX(1,ooq)=coords(1,connect(tempooq,ooe));
        YY(1,ooq)=coords(2,connect(tempooq,ooe));
        dd(1,ooq)=stressplot(connect(tempooq,ooe));
    end
    fill(XX,YY,dd);
    hold on;
    shading interp;
    colorbar;
    hold on;

end

title('Stress_{11} on undeformed mesh');
% title(['hydrostatic pressure on undeformed mesh']);
axis equal
%     xlim([0 44])
%     ylim([0 60])
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
set(gca,'XColor','white')
set(gca,'YColor','white')
set(gca,'Color','white')
end

function [M,R] = L2_global(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,target)

%
%   Assemble the global mass matrix

npoints = numberofintegrationpoints(ncoord,nelnodes,elident);


M = zeros(nnode,nnode);
R = zeros(nnode,1);
lmncoord = zeros(ncoord,maxnodes);
lmntarget = zeros(1,npoints);
%
%   Loop over all the elements
%
for lmn =  1 : nelem
    %
    %   Extract coords of nodes, DOF for the current element
    %
    for a =  1 : nelnodes(lmn)
        for i = 1 : ncoord

            lmncoord(i,a) = coords(i,connect(a,lmn));

        end
    end
    for ooi=1:npoints
        lmntarget(ooi)=target(ooi,lmn);
    end
    n = nelnodes(lmn);
    ident = elident(lmn);
    [mel,rel] = L2_el(ncoord,ndof,n,ident,lmncoord,lmntarget);
    %
    %   Add the current element stiffness to the global stiffness
    for a = 1 : nelnodes(lmn)
        rw = (connect(a,lmn));
        R(rw) = R(rw) + rel((a));
    end

    for a = 1 : nelnodes(lmn)
        for b = 1 : nelnodes(lmn)
            rw = connect(a,lmn);
            cl = connect(b,lmn);
            M(rw,cl) = M(rw,cl) + mel(a,b);

        end
    end
end

end




%================= ELEMENT MASS MATRIX ====================================
%
function [mel,rel] = L2_el(ncoord,ndof,nelnodes,elident,coord,target)
%
%  Assemble the (consistent) element mass matrix
%
%    Arguments:
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords[i,a]        ith coord of ath node
%      materialprops      Material properties passed on to constitutive procedures
%      displacement[i,a]  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi[i,inpt]         ith local coord of integration point no. intpt
%      w[intpt]           weight for integration point no. intpt
%      N[a]               Shape function associated with ath node on element
%      dNdxi[a,i]         Derivative of ath shape function wrt ith local coord
%      dNdx[a,i]          Derivative of ath shape function wrt ith global coord
%      dxdxi[i,j]         Derivative of ith global coord wrt jth local coord
%      dxidx[i,j]         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian


npoints = numberofintegrationpoints(ncoord,nelnodes,elident);

xi = zeros(ncoord,1);
dxdxi = zeros(ncoord,ncoord);
mel = zeros(nelnodes,nelnodes);
rel = zeros(nelnodes,1);
%
%  Set up integration points and weights
%
xilist = integrationpoints(ncoord,nelnodes,npoints,elident);
w = integrationweights(ncoord,nelnodes,npoints,elident);
%
%  Loop over the integration points
%
for intpt = 1 : npoints

    %     Compute shape functions and derivatives wrt local coords
    %
    for i = 1 : ncoord
        xi(i) = xilist(i,intpt);
    end
    N = shapefunctions(nelnodes,ncoord,elident,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi);
    %
    %     Compute the jacobian matrix and its determinant
    %
    for i = 1 : ncoord
        for j = 1 : ncoord
            dxdxi(i,j) = 0.;
            for a = 1 : nelnodes
                dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
            end
        end
    end
    dt = det(dxdxi);

    for a = 1 : nelnodes
        for b = 1 : nelnodes
            for i = 1 : ndof
                mel(a,b) = mel(a,b) + N(b)*N(a)*w(intpt)*dt;
            end
        end
    end

    for a = 1 : nelnodes
        rel(a) = rel(a) + N(a)*target(intpt)*w(intpt)*dt;
    end

end


end

%
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%
function n = numberofintegrationpoints(ncoord,nelnodes,elident)

if (ncoord == 1)
    n = nelnodes;
elseif (ncoord == 2)

    if (nelnodes == 4)
        n = 4;
    end


end
end
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
function xi = integrationpoints(ncoord,nelnodes,npoints,elident)

xi = zeros(ncoord,npoints);
%
%  1D elements
%
if (ncoord == 1)
    if (npoints==1)
        xi(1,1) = 0.;
    elseif (npoints == 2)
        xi(1,1) = -0.5773502692;
        xi(1,2) = -xi(1,1);
    elseif (npoints == 3)
        xi(1,1) = -0.7745966692;
        xi(1,2) = 0.0;
        xi(1,3) = -xi(1,1);
    end
    %
    %  2D elements
    %
elseif (ncoord == 2)
    if ( nelnodes==4  )

        if (npoints == 1)
            xi(1,1) = 0.;
            xi(2,1) = 0.;
        elseif (npoints == 4)
            xi(1,1) = -0.5773502692;
            xi(2,1) = xi(1,1);
            xi(1,2) = -xi(1,1);
            xi(2,2) = xi(1,1);
            xi(1,3) = xi(1,1);
            xi(2,3) = -xi(1,1);
            xi(1,4) = -xi(1,1);
            xi(2,4) = -xi(1,1);

        end
    end
end
end

%
%================= INTEGRATION WEIGHTS ==================================
%
%   Defines integration weights w_i
%
function w = integrationweights(ncoord,nelnodes,npoints,elident)

w = zeros(npoints,1);

%
%  1D elements
%
if (ncoord == 1)
    if (npoints == 1)
        w(1) = 2.;
    elseif (npoints == 2)
        w = [1.,1.];
    elseif (npoints == 3)
        w = [0.555555555,0.888888888,0.555555555];
    end
    %
    %  2D elements
    %
elseif (ncoord == 2)
    %

    %    Rectangular element
    %
    if ( nelnodes==4 )

        if (npoints == 1)
            w(1) = 4.;
        elseif (npoints == 4)
            w = [1.,1.,1.,1.];
        elseif (npoints == 9 )
            w1D = [0.555555555,0.888888888,0.55555555555];
            for j = 1:3
                for i = 1:3
                    n = 3*(j-1)+i;
                    w(n) = w1D(i)*w1D(j);
                end
            end
        end
    end

end
end

%
%================= SHAPE FUNCTIONS ==================================
%
%        Calculates shape functions for various element types
%
function N = shapefunctions(nelnodes,ncoord,elident,xi)


N = zeros(nelnodes,1);
%
%  1D elements
%
if (ncoord == 1)
    if (nelnodes==2)
        N(1) = 0.5*(1.+xi(1));
        N(2) = 0.5*(1.-xi(1));
    end
    %
    %  2D elements
    %
elseif (ncoord == 2)
    %    Rectangular element
    %
    if ( nelnodes == 4 )
        N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
        N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
        N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
        N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
    end

end

end

%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
function dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi)

dNdxi = zeros(nelnodes,ncoord);
%
% 1D elements
%
if (ncoord == 1)
    if (nelnodes==2)
        dNdxi(1,1) = 0.5;
        dNdxi(2,1) = -0.5;
    elseif (nelnodes == 3)
        dNdxi(1,1) = -0.5+xi(1);
        dNdxi(2,1) =  0.5+xi(1);
        dNdxi(3,1) = -2.*xi(1);
    end
    %
    %  2D elements
    %
elseif (ncoord == 2)
    %    Rectangular element
    %
    if ( nelnodes == 4 )
        dNdxi(1,1) = -0.25*(1.-xi(2));
        dNdxi(1,2) = -0.25*(1.-xi(1));
        dNdxi(2,1) = 0.25*(1.-xi(2));
        dNdxi(2,2) = -0.25*(1.+xi(1));
        dNdxi(3,1) = 0.25*(1.+xi(2));
        dNdxi(3,2) = 0.25*(1.+xi(1));
        dNdxi(4,1) = -0.25*(1.+xi(2));
        dNdxi(4,2) = 0.25*(1.-xi(1));
    end

end
end

%====================== No. nodes on element faces ================
%
%   This procedure returns the number of nodes on each element face
%   for various element types.  This info is needed for computing
%   the surface integrals associated with the element traction vector
%
function n = nfacenodes(ncoord,nelnodes,elident,face)
if (ncoord == 2)
    if ( nelnodes == 4)
        n = 2;
    end

end
end
%======================= Lists of nodes on element faces =============
%
%    This procedure returns the list of nodes on an element face
%    The nodes are ordered so that the element face forms either
%    a 1D line element or a 2D surface element for 2D or 3D problems
%
function list = facenodes(ncoord,nelnodes,elident,face)

i4 = [2,3,4,1];

list = zeros(nfacenodes(ncoord,nelnodes,face),1);

if (ncoord == 2)
    if (nelnodes==4)
        list(1) = face;
        list(2) = i4(face);
    end
end
end
%
%====================== ELEMENT DISTRIBUTED LOAD VECTOR ==============
%
function r = eldload(ncoord,ndof,nfacenodes,elident,coords,traction)

npoints = numberofintegrationpoints(ncoord-1,nfacenodes);
xi = zeros(ncoord-1,1);
dxdxi = zeros(ncoord,ncoord-1);
r = zeros(ndof*nfacenodes,1);

xilist = integrationpoints(ncoord-1,nfacenodes,npoints);
w = integrationweights(ncoord-1,nfacenodes,npoints);

for intpt = 1:npoints

    for i = 1:ncoord-1
        xi(i) = xilist(i,intpt);
    end

    N = shapefunctions(nfacenodes,ncoord-1,elident,xi);
    dNdxi = shapefunctionderivs(nfacenodes,ncoord-1,elident,xi);
    %
    %     Compute the jacobian matrix && its determinant
    %
    for i = 1:ncoord
        for j = 1:ncoord-1
            dxdxi(i,j) = 0.;
            for a = 1:nfacenodes
                dxdxi(i,j) = dxdxi(i,j) + coords(i,a)*dNdxi(a,j);
            end
        end
    end
    if (ncoord == 2)
        dt = sqrt(dxdxi(1,1)^2+dxdxi(2,1)^2);
    end
    for a = 1:nfacenodes
        for i = 1:ndof
            row = ndof*(a-1)+i;
            r(row) = r(row) + N(a)*traction(i)*w(intpt)*dt;
        end
    end
end
end
%
%====================== Assemble the global stiffness matrix =================
%
function Stif = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs)
%
%   Assemble the global stiffness matrix
%

Stif = zeros(ndof*nnode,ndof*nnode);
lmncoord = zeros(ncoord,maxnodes);
lmndof = zeros(ndof,maxnodes);
%
%   Loop over all the elements
%
for lmn = 1:nelem
    %
    %   Extract coords of nodes, DOF for the current element
    %
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
            lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
    end
    n = nelnodes(lmn);
    ident = elident(lmn);
    kel = elstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof);
    %
    %   Add the current element stiffness:the global stiffness
    %
    for a = 1:nelnodes(lmn)
        for i = 1:ndof
            for b = 1:nelnodes(lmn)
                for k = 1:ndof
                    rw = ndof*(connect(a,lmn)-1)+i;
                    cl = ndof*(connect(b,lmn)-1)+k;
                    Stif(rw,cl) = Stif(rw,cl) + kel(ndof*(a-1)+i,ndof*(b-1)+k);
                end
            end
        end
    end
end
end
%
%===================== Assemble the global traction vector =============
%
function r = globaltraction(ncoord,ndof,nnodes,ndload,coords,nelnodes,elident,connect,dloads,dofs)

r = zeros(ndof*nnodes,1);
traction = zeros(ndof,1);

for load = 1:ndload
    %
    %     Extract the coords of the nodes on the appropriate element face
    %
    lmn = dloads(1,load);
    face = dloads(2,load);
    n = nelnodes(lmn);
    ident = elident(lmn);
    nfnodes = nfacenodes(ncoord,n,ident,face);
    nodelist = facenodes(ncoord,n,ident,face);
    lmncoord = zeros(ncoord,nfnodes);
    for a = 1:nfnodes
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(nodelist(a),dloads(1,load)));
        end
        for i = 1:ndof
            lmndof(i,a) = dofs(ndof*(connect(nodelist(a),dloads(1,load))-1)+i);
        end
    end
    %
    %    Compute the element load vector
    %
    for i = 1:ndof
        traction(i) = dloads(i+2,load);
    end

    rel = eldload(ncoord,ndof,nfnodes,ident,lmncoord,traction);
    %
    %    Assemble the element load vector into global vector
    %
    for a = 1:nfnodes
        for i = 1:ndof
            rw = (connect(nodelist(a),dloads(1,load))-1)*ndof+i;
            r(rw) = r(rw) + rel((a-1)*ndof+i);
        end
    end

end
end


function plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,color)
% Function to plot a mesh.
f2D_4 = [1,2,3,4];
hold on
if (ncoord==2)  % Plot a 2D mesh
    for lmn = 1:nelem
        for i = 1:nelnodes(lmn)
            x(i,1:2) = coords(1:2,connect(i,lmn));
        end
        %         scatter(x(:,1),x(:,2),'MarkerFaceColor','r');
        if (nelnodes(lmn)==4)
            patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor',color);
        end
    end
end
axis equal
hold off
end


%================= Material Stress ==================================
%
%   Computes stress sigma_{ij} given strain epsilon_{ij}
%
function stress = materialstress(ndof,ncoord,strain,materialprops)

C = materialstiffness(ndof,ncoord,strain,materialprops);
stress = zeros(ndof,ncoord);
for i = 1 : ndof
    for j = 1 : ncoord
        for k = 1 : ndof
            for l = 1: ncoord
                stress(i,j) = stress(i,j) + C(i,j,k,l)*strain(k,l);
            end
        end
    end
end
end

%================= Material Stiffness ==================================
%
%    Computes elasticity tensor C_{ijkl} = shear modulus and Poissons ratio
%    Currently coded either for plane strain, plane stress or general 3D.
%
function C = materialstiffness(ndof,ncoord,strain,materialprops)

mu = materialprops(1);
nu = materialprops(2);

C = zeros(ndof,ncoord,ndof,ncoord);


%  planestrain = 0 => plane stress, planestrain = 1 => plane strain
planestrain = materialprops(3);

for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                if (planestrain==1)
                    if (i==j && k==l)
                        C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-2*nu);
                    end
                else
                    if (i==j && k==l)
                        C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-nu);
                    end
                end
                if (i==l && k==j)
                    C(i,j,k,l) = C(i,j,k,l)+mu;
                end
                if (i==k && j==l)
                    C(i,j,k,l) = C(i,j,k,l)+mu;
                end
            end
        end
    end
end


end
%
%================= ELEMENT STIFFNESS MATRIX ================================
%
function kel = elstif(ncoord,ndof,nelnodes,elident,coord,materialprops,displacement)
%
%  Assemble the element stiffness
%
%    Arguments;
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures
%      displacement(i,a)  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness
%
%
npoints = numberofintegrationpoints(ncoord,nelnodes,elident);
dNdx = zeros(nelnodes,ncoord);
dxdxi = zeros(ncoord,ncoord);
strain = zeros(ndof,ncoord);
kel = zeros(ndof*nelnodes,ndof*nelnodes);
%
%  Set up integration points && weights
%
xilist = integrationpoints(ncoord,nelnodes,npoints,elident);
w = integrationweights(ncoord,nelnodes,npoints,elident);
%
%  Loop over the integration points
%
for intpt = 1:npoints

    %     Compute shape functions && derivatives wrt local coords
    %
    for i = 1:ncoord
        xi(i) = xilist(i,intpt);
    end
    N = shapefunctions(nelnodes,ncoord,elident,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi);


    %
    %     Compute the jacobian matrix && its determinant
    %
    for i = 1:ncoord
        for j = 1:ncoord
            dxdxi(i,j) = 0.;
            for a = 1:nelnodes
                dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
            end
        end
    end

    dxidx = inv(dxdxi);
    dt = det(dxdxi);
    %
    %     Convert shape function derivatives:derivatives wrt global coords
    %
    for a = 1:nelnodes
        for i = 1:ncoord
            dNdx(a,i) = 0.;
            for j = 1:ncoord
                dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
            end
        end
    end
    %
    %     Compute the (infinitesimal) strain by differentiating displacements
    %     This step is not really necessary for linear elasticity calculations
    %     where stiffness is independent of strain.  It is included:allow
    %     extension:nonlinear materials later.
    %
    for i = 1:ncoord
        for j = 1:ncoord
            strain(i,j) = 0.;
            for a = 1:nelnodes
                strain(i,j) = strain(i,j) + 0.5*(displacement(i,a)*dNdx(a,j)+displacement(j,a)*dNdx(a,i));
            end
        end
    end
    %
    %     Compute the material tangent stiffness (d stress/d strain)
    %     ds/de is just C_ijkl for linear elasticity - this notation is used
    %     to allow extension to nonlinear problems
    %
    dsde = materialstiffness(ndof,ncoord,strain,materialprops);
    %
    %     Compute the element stiffness
    %
    for a = 1:nelnodes
        for i = 1:ndof
            for b = 1:nelnodes
                for k = 1:ndof
                    row = ndof*(a-1)+i;
                    col = ndof*(b-1)+k;
                    for j = 1:ncoord
                        for l = 1:ncoord
                            kel(col,row) = kel(col,row) + dsde(i,j,k,l)*dNdx(b,l)*dNdx(a,j)*w(intpt)*dt;
                        end
                    end
                end
            end
        end
    end
end

end



