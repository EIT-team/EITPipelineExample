%% Dependencies

%see readme.md for more info

%iso2mesh
%mesher
%supersolver
%reconstruction
%paraview to look at .vtks

%% Create Mesh Geometry

% Taking Neonatal scalp as an example, we have the .stl file created in
% solidworks, and the desired electrode positions in the stl coordinates

%THIS IS IMPORTANT! This is the resolution of the INR file,i.e. the number
%of pixels/voxels per mm. A finer resolution will result in a bigger .inr
%file and take longer to load into the mesher, but it will represent the
%geometery better. The pixel scale must be carried through to the stl2inr
%and the mesher parameter file, otherwise the mesh and electrodes will all
%be in different coordinates
pixel_scale=1;

% These were defined inside the solidworks file and saved to text file
% directly. As long as elec_pos is Nx3 (x,y,z)
elec_pos=dlmread('resources\NNelecposorig.txt');

% create INR from STL - double check alignment of electrodes
[ full_mask,elec_pos_new_sc ] = stl2inr( 'resources\NNscalp.stl',pixel_scale,elec_pos );

%set parameters
P=getmesherparam;
P.pixel_scale_mm=pixel_scale;
P.facet_distance_mm=2;
P.cell_fine_size_mm=3;
P.cell_coarse_size_mm=8;
P.electrode_radius_mm=4;
P.cell_size_electrodes_mm=2;

%set optimisations
P.opt.exude_opt=1;
P.opt.lloyd_opt=1;
P.opt.odt_opt=1;
P.opt.perturb_opt=1;

P.save.save_nodes_tetra=1; %do save mesh so we can load it into matlab after
writemesherparam('resources/NNscalp_param.txt',P);

% Run mesher
runmesher('resources/NNscalp.inr','resources/NNscalp_elecINRpos.txt',...
    'resources/NNscalp_param.txt','output/','NNexample')

%check the mesh by loading NNexample.vtu

% save mesh as matlab format for future use once this is done you do not
% need to run the mesher again unless you change geometry or electrode
% locations (only if you have refined at the electrode sites)
%%
M=loadmesh('output/NNexample');
%DisplayBoundaries(M)

%% Run forward solver

%get protocol. each row is a separate measurement. with the columns CS+ CS-
%V+ V-
Prt=dlmread('resources/NN_Prt_full.txt');

cond_back=0.35; %background conductivity

%create conductivity vector - currently it is all set to 1 to identify the
%tissue number
mat_ref=zeros(size(M.mat_ref));
mat_ref(M.mat_ref==1) = cond_back;
% mat_ref(M.mat_ref==2) = 0.5; %for example if you had a second tissue and
% you wanted to specify a differerent conductivity to it

% get the structures with defaults
[Mesh,Fem,Fwd,Inv,Sol] = supersolver_init(M.Nodes,M.Tetra,mat_ref,M.elec_pos,M.gnd_pos,Prt);
%change the electrode diameters etc.
Fem.current                         =   300e-6; % uA
Fem.elec_diam                       =   8e-3;% electrode diameter in meters
Fem.zc                              =   200 ; % contact impedance (could be a vector at the length of the number of electrodes);

%Setup system
[Mesh,Fem,Fwd,Inv,Sol] = supersolver_setup(Mesh,Fem,Fwd,Inv,Sol);
% run fwd
[Mesh,Fem,Fwd,Inv,Sol,Data] = supersolver_runfwd(Mesh,Fem,Fwd,Inv,Sol);
% get jac
[Mesh,Fem,Fwd,Inv,Sol,J] = supersolver_solve(Mesh,Fem,Fwd,Inv,Sol);

%% Create inverted jacobian
% create a hexagonised jacobian - this evens out the sensitivity across the
% mesh. Should be at least an order of magnitude less hexes than tetras but
% you dont want any hole in it
[Mesh_hex,J_hex] = convert_fine2coarse(M.Tetra,M.Nodes,J,0.006);

% ###########################
% The code before this only needs to be done once per mesh and protocol, so
% the remainder of the code is creating simulated data and making the
% images
% ######

%% Create perturbation data

% now we have our Jacobian [J], and our baseline simulated data
% [data.bnd_v], but we need our difference data to make an image. So we
% need to create a new set of boundary voltages 

pert_loc= [0.04 0.05, 0.06]; %xyz of centre of perturbation
pert_rad=0.01; % radius
pert_cond=0.01; %conductivity of blob

% find the distance of the centre of all tetrahedra from location
cnts_tetra=(M.Nodes(M.Tetra(:,1),:)+M.Nodes(M.Tetra(:,2),:)+M.Nodes(M.Tetra(:,3),:)+M.Nodes(M.Tetra(:,4),:))/4;
dist_tetra=cnts_tetra-repmat(pert_loc,length(cnts_tetra),1);
dist_tetra=sum(dist_tetra.^2,2).^0.5;

%create index of those inside the blob -1 is normal, 2 is pert
mat_ref_idx=M.mat_ref;
mat_ref_idx(dist_tetra<pert_rad)=2;

% create new conductivity vector
mat_ref_pert=zeros(size(M.mat_ref));
mat_ref_pert(mat_ref_idx==1) = 0.35;
mat_ref_pert(mat_ref_idx==2) = pert_cond;

%make mesh for display purposes
Mpert=M;
Mpert.Tetra=M.Tetra((mat_ref_idx==2),:);
Mpert.Nodes=M.Nodes;
figure

hold on
h=DisplayBoundaries(Mpert);
set(h,'facecolor','r')
set(h,'facealpha',0.5)
set(h,'EdgeAlpha',0);

h2=DisplayBoundaries(M);
set(h2,'facealpha',0);
set(h2,'EdgeColor',[0.5,0.5,0.5]);
set(h2,'EdgeAlpha',0.3);

hold off
title('Perturbation location');

% re run forward with this new conducivity vector
Solp=Sol;
Solp.ref=mat_ref_pert;


[Meshp,Femp,Fwdp,Invp,Solp] = supersolver_setup(Mesh,Fem,Fwd,Inv,Solp);
% run fwd
[Meshp,Femp,Fwdp,Invp,Solp,Datap] = supersolver_runfwd(Meshp,Femp,Fwdp,Invp,Solp);

%p.s. you could just cheat and multiply J*mat_ref_pert but that assumes
%linearity

%% Create difference voltages and inverse jacobian 

dV=(Datap.bnd_v-Data.bnd_v)'; % this is the data you would get from the EIT system, it needs to look like this

%% Make Image

[U,S,V] = svd(J_hex,'econ'); % do this here to avoid doing it over and over again
n=3e-6*randn(500,length(dV)); % estimate of the noise in the measurments, we dont have any here so just make it up
lambda = logspace(-20,-1,200); % the range of regularisation parameters to try

% reconstruct image
[Sigma,X,sv_i]=eit_recon_tik0(dV,J_hex,'output/EITrec',n,U,S,V,lambda);

%write vtk files to view in paraview
writeVTKcell_hex('output/EITrec_Tik0',Mesh_hex.Hex,Mesh_hex.Nodes,X); % no noise based correction
writeVTKcell_hex('output/EITrec_Tik0_correction',Mesh_hex.Hex,Mesh_hex.Nodes,Sigma); % with correction added

% NOTE THE MESH SIZES AND PARAMETERS CHOSEN ARE FOR DEMONSTRATION PURPOSES
% ONLY!

%% doing it without cv


X2=V*(diag(1./sv_i)*U'*dV');

UtNoise = U'*n';

SD2 = std(V*(diag(1./sv_i)*UtNoise),0,2);


Sigma2=X2./SD2;

writeVTKcell_hex('output/EITrec_Tik0_2',Mesh_hex.Hex,Mesh_hex.Nodes,X2); % no noise based correction
writeVTKcell_hex('output/EITrec_Tik0_correction_2',Mesh_hex.Hex,Mesh_hex.Nodes,Sigma2); % with correction added
%%
X3lhs=V*(diag(1./sv_i)*U');

X3=X3lhs*dV';

Sigma3=(X3lhs*dV')./SD2;
Sigma4=(X3lhs./SD2)*dV';

SigLhs=(X3lhs./SD2);
tic
Sigmafast=SigLhs*dV';
toc
%% get slice

    cnts=zeros(length(Mesh_hex.Hex),3);
    for i=1:8
        cnts=cnts+Mesh_hex.Nodes(Mesh_hex.Hex(:,i),:)/8;
    end
    
    %convert centers to xyz sub - not sure how to plot incomplete slices in
    %labview
    
    
    




