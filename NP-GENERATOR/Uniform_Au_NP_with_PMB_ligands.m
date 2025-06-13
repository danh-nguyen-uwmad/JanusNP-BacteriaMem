clear all
close all
clc

tname = 'Janus_core.data';
m  = 0;            % MUS length
%%%%%%%%%%%%%%%%%%%%%%%%%%%
box = 100;  % define the box length Angstrom

%% Model for NP surface
%n_density = 0.0476; % density of the alkyl chain 0.0476 -- originial number
n_density = 0.015; % change number here for the PMB attachment
r_p = 50; % radius of NP (Angstrom)
n_p = round(4*pi*r_p^2*n_density);

%% Genetate the sphere
index        = [linspace(0,n_p-1,n_p)]';
golden_angle = pi*(3-sqrt(5));
index_phobic = find(mod(index,2)==0);
index_philic = find(mod(index,2)==1);
theta        = golden_angle*index_phobic;
z            = [linspace(1-1/n_p,1/n_p-1,n_p)]';
z_phobic = z(index_phobic);
z_philic = z(index_philic);

radius = (1-z_phobic.^2).^(1/2);
x = radius.*cos(theta);
y = radius.*sin(theta);
z_phobic = z_phobic;
XYZ=r_p*[x y z_phobic];

n_surf=size(XYZ,1);
number=1:n_surf;
NP_surf=zeros(n_surf,5);
NP_surf(:,1)=(1:n_surf)';
NP_surf(:,2)=2;
NP_surf(:,3:5)=XYZ;

%% Location of Sulfur beads
% Add the core
lineID=9;
fid = fopen('./Au_FCC_matrix.lammpstrj');
Fccmatrix = textscan(fid, '%f %f %f %f','HeaderLines',lineID);
fclose(fid);
MM=cell2mat(Fccmatrix);

xyz_p=MM(:,2:4);
center=mean(xyz_p);
xyz_p(:,1)=xyz_p(:,1)-center(1);
xyz_p(:,2)=xyz_p(:,2)-center(2);
xyz_p(:,3)=xyz_p(:,3)-center(3);
rescale=1.0;
%rescale=1.82;
xyz_p=xyz_p*rescale;
NP_new=[];
N1=0;

for i=1:size(xyz_p,1)
    dist=norm(xyz_p(i,:));
    if dist<(r_p-2)
        N1=N1+1;
        NP_new=[NP_new; xyz_p(i,:)];
    end
end


XYZ = [XYZ; NP_new];
display(['The S number on surface is ' num2str(n_surf)]);

NP_inner_core = zeros(N1,5);
NP_inner_core(:,1) = (1:N1);
NP_inner_core(:,2) = 1;
NP_inner_core(:,3:5) = NP_new;
NP = [NP_surf;NP_inner_core];
display(['The Au number of the core is ' num2str(size(NP_inner_core,1))]);

% The whole NP core
N = size(NP,1);
n_g=n_surf;
rng(10);
a=randperm(n_surf);
r=a(1:n_g);
r=sortrows(r,1);
NP_surf(r,2)=3;  % grafting sites: type information

% Store the NP as the lammps position
atom=zeros(N,7);
for i=1:N
    atom(i,1)=i;   % id
    atom(i,2)=1;   % mol
    atom(i,3)=1;   % type
    atom(i,4)=0.0;   % charge
    atom(i,5:7)=NP(i,3:5); % position
end

atom(NP_surf(:,1),3)=NP_surf(:,2);
[n_g,tempN]=size(find(atom(:,3)==3));
display([n_surf n_g]);
r=find(atom(:,3)==3);
rr=find(atom(:,3)==2);

% generate the chains on the graft sites
chain=zeros( m*n_g,7);
bond =zeros(m*n_g ,4);  % bond list
angle=zeros((m-1)*n_g,5); % angle list
dihedral=zeros((m-3)*n_g,6); % dihedral list

for i=1:n_g
    id=r(i);         % the atom ID of the target bead in the NP core
    atom(id,3)= 2;   % change the graft site back to the original site
    site=atom(id,5:7);
    r1=norm(site);
    r2=acos(site(3)/r1);
    r3=atan2(site(2),site(1));
    
    % write the position
    for j=1:m
        no = j+(i-1)*m;
        chain(no,1)=no+N;
        chain(no,2)=1+i;
        chain(no,3)=3;
        chain(no,4)=0.0;
        chain(no,5)=(r1+j*4.7)*sin(r2)*cos(r3); % 4.7 must be the polymer bond length
        chain(no,6)=(r1+j*4.7)*sin(r2)*sin(r3);
        chain(no,7)=(r1+j*4.7)*cos(r2);
        % 		if j == m
        % 		chain(no,3)= 4; % determine tomorrow
        %         end
    end
    
    % write the bonds
    for j=1:m
        no=j+(i-1)*m;
        bond(no,1)=no;
        bond(no,2)=1;
        if j==1
            head = id;
            tail = j+(i-1)*m+N;
            bond(no,2) = 2;
        else
            head=j-1+(i-1)*m+N;
            tail=head+1;
        end
        bond(no,3)=head;
        bond(no,4)=tail;
    end
    
    % write the angles
    for j=1:(m-1)
        no=j+(i-1)*(m-1);
        angle(no,1)=no;
        angle(no,2)=1;
        if j <=m-1
            if j==1
                head=id;
                middle=j+(i-1)*m+N;
                tail=middle+1;
            else
                head=j-1+(i-1)*m+N;
                middle=head+1;
                tail=middle+1;
            end
        else
            head = j+2+(i-1)*m+N;
            middle = head + 1;
            tail = middle +1;
        end
        angle(no,3)=head;
        angle(no,4)=middle;
        angle(no,5)=tail;
        
    end
    
    % should ad the lipid angle information here
    % write the dihedrals # dihedral is the not useful in this step
    
    for j=4:m
        no=j-3+(i-1)*(m-3);
        dihedral(no,1)=no;
        dihedral(no,2)=1;
        head=j-3+(i-1)*m+N;
        middle1=head+1;
        middle2=middle1+1;
        tail=middle2+1;
        dihedral(no,3)=head;
        dihedral(no,4)=middle1;
        dihedral(no,5)=middle2;
        dihedral(no,6)=tail;
    end
end
% Add one line here for total of alkyl beads
n_alkyl = size(chain,1);

%% Hydrophilic ligands (PMB)
z_philic = z(index_philic);
golden_angle2 = pi*(3-sqrt(5));
theta2        = golden_angle2*index_philic;
radius2 = (1-z_philic.^2).^(1/2);
x2 = radius2.*cos(theta2);
y2 = radius2.*sin(theta2);
z_philic = z_philic;
XYZ2 = r_p*[x2 y2 z_philic];
n_p2               = size(XYZ2,1);

PMB_Pos=zeros(n_p2,7);
PMB_Pos(:,1) = (1+m*n_g+N:n_p2+m*n_g+N);
PMB_Pos(:,2) = 1;
PMB_Pos(:,3) = 4;
PMB_Pos(:,4) = 1.0;
PMB_Pos(:,5:7) = XYZ2;

all_data=[atom; chain; PMB_Pos];
nu_atom = size(all_data,1);
all_data(:,5)=all_data(:,5); %+box/2;
all_data(:,6)=all_data(:,6); %+box/2;
all_data(:,7)=all_data(:,7); %+box/2;

%% Convert to GROMACS
%% This section is to add PMB to the Janus NP
% Read AuNP coordinate, then save xyz for Au
AuNP_core_xyz = NP_inner_core(:,3:5)/10;
Num_Au = size(AuNP_core_xyz, 1);
AuNP_atom = zeros(N,7);

% Center of mass of AuNP core
AuNP_COM_x = mean(AuNP_core_xyz(:,1));
AuNP_COM_y = mean(AuNP_core_xyz(:,2));
AuNP_COM_z = mean(AuNP_core_xyz(:,3));

AuNP_COM = [ AuNP_COM_x AuNP_COM_y AuNP_COM_z ]

XYZ3 = PMB_Pos(:,5:7) * 1.089;
n_p             = size(XYZ3,1);

% Rescale the XYZ from Angstrong to nm
XYZ3 = XYZ3/10;

% PLot the points to make sure they are good
scatter3(XYZ2(:,1)/10,XYZ2(:,2)/10,XYZ2(:,3)/10);
hold on
scatter3(XYZ3(:,1),XYZ3(:,2),XYZ3(:,3));

%% Save the coordinate for all connecting point C1 = index from 1 to n_p2
C1_Pos=zeros(n_p2,7);
C1_Pos(:,1) = (1: n_p2);
C1_Pos(:,2) = 1;
C1_Pos(:,3) = 4;
C1_Pos(:,4) = 1.0;
C1_Pos(:,5:7) = XYZ3;

%% Add PMB at each point of PMB_pos:
% 1. Read PMB xyz file;
% 2. Rotate the PMB x, y, z by the angle between vector of (COM-point) and PMB vector;
% 3. Translocate the C1 of PMB to the connecting point from C1 atom (PMB)to
% that connecting point.

% Read PMB coordinate, then save xyz for PMB
PMB = readtable('PMB.xlsx');
PMB_CG_name = PMB(:,1);
PMB_CG_name = PMB_CG_name{:,:};
PMB_xyz = PMB(:,3:5);
PMB_xyz = PMB_xyz{:,:};
% Vector connecting C1-C1 (atom 1 - atom 21) of PMB
vec_PMB = PMB_xyz(21,:) - PMB_xyz(1,:);

% Number Of PMBs to add:
num_PMB = n_p2;
num_atom_per_PMB = size(PMB_xyz,1);
Total_atom_PMB = num_PMB * num_atom_per_PMB;

PMB_atom = zeros(Total_atom_PMB,7); % To save coordinates of all PMBs
PMB_atom(:,1) = (1 + nu_atom : Total_atom_PMB + nu_atom); %id of PMB atoms
PMB_atom(:,3) = 2; % type of PMB (AuNP is type 1)
PMB_atom(:,4) = 0.0; % charge
PMB_CG_type = []; % bead type of PMB

% Rotate, translate and save PMB molecule coordinates
for i = 1 : num_PMB
        vec_COM_point = C1_Pos(i,5:7) - AuNP_COM;
        
        p0 = vec_PMB';
        p1 = vec_COM_point';
        %https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
        CosTheta = max(min(dot(p0,p1)/(norm(p0)*norm(p1)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));

        u = cross(p1,p0);
        u = u/norm(u);
        %https://www.mathworks.com/matlabcentral/answers/500030-rotate-a-3d-data-cloud-to-align-with-one-axis
        rotMat=eye(3)*cosd(ThetaInDegrees)+sind(ThetaInDegrees)*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(ThetaInDegrees))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];
        for k = 1:size(PMB_xyz, 1)
             PMB_xyz_rot(k,1:3) = PMB_xyz(k,1:3)*rotMat;
             %PMB_xyz_rot(i,1:3)/norm(PMB_xyz_rot(i,1:3));
        end 

        % Move the rotated PMB to the connecting point
        PMB_C1_point = PMB_xyz_rot(1,1:3); % first atom of PMB = C1
        Connect_point = C1_Pos(i,5:7);
        
        dx = PMB_C1_point(1) - Connect_point(1);
        dy = PMB_C1_point(2) - Connect_point(2);
        dz = PMB_C1_point(3) - Connect_point(3);
        % Translation vector
        ds = [dx dy dz];

        % Translate the rotated PMB to the connecting point
        PMB_xyz_final = PMB_xyz_rot  - ds;
        
        PMB_atom((i-1)*num_atom_per_PMB+1 : i*num_atom_per_PMB, 2) = i+1;   % mol
        PMB_atom((i-1)*num_atom_per_PMB+1 : i*num_atom_per_PMB,5:7) = PMB_xyz_final; % position
        
        for j = 1 : num_atom_per_PMB
            temp = string(PMB_CG_name{j});
            PMB_CG_type  = vertcat(PMB_CG_type, temp);
        end 
        
end 

%% Plot to check the rotation
scatter3(PMB_xyz_rot(:,1), PMB_xyz_rot(:,2), PMB_xyz_rot(:,3),'x');
hold on
scatter3(PMB_xyz(:,1), PMB_xyz(:,2), PMB_xyz(:,3),'o');
%scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3), 'x'); scatter3(AuNP_COM_x,
%AuNP_COM_y, AuNP_COM_z);

%% Sum all atoms

all_data_combine = [all_data; PMB_atom];
nu_atom_combine = size(all_data_combine,1);

%% write the .gro file for Gromacs

fid = fopen('Janus_PMB.gro','w');
fprintf(fid,'Coordination of one JANUS_PMB Particle\n');
fprintf(fid,'%7d\n', nu_atom_combine);

s1=0;
s2=0;
s3=0;
s4=0;
s5=0;
for i=1:size(all_data,1)
      
    if all_data(i,3)==1
        s1=s1+1;
        fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',1,'MOL',strcat('C5'),all_data(i,1), all_data(i,5:7)/10);
    elseif all_data(i,3)==2
        s2=s2+1;
        fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',1,'MOL',strcat('N0'),all_data(i,1), all_data(i,5:7)/10);
    elseif all_data(i,3)==3
        s3=s3+1;
        fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',1,'MOL',strcat('C1'),all_data(i,1), all_data(i,5:7)/10);
    elseif all_data(i,3)==4
        s4=s4+1;
        fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',1,'MOL',strcat('N0'),all_data(i,1), all_data(i,5:7)/10);
    elseif all_data(i,3)==5
        s5=s5+1;
        fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',1,'MOL',strcat('P4'),all_data(i,1), all_data(i,5:7)/10);
    end
    
end



% Write for PMB molecules
for i=1:size(PMB_atom,1)
    
fprintf(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',PMB_atom(i,2),'PMB',PMB_CG_type(i), PMB_atom(i,1), PMB_atom(i,5:7));
    
end

box = 20;

fprintf(fid,'    %10.5f%10.5f%10.5f\n', box, box, box);
fclose(fid);

%% This part is for finding the bonding between Au-Au, Au-S and Au-PMB

% For Au-Au bonding
Au_inner_total = size(NP_inner_core,1);
Au_first = [];
Au_second = [];
Au_Au = [];

for i = (n_surf+1) : (n_surf + Au_inner_total-1)
    for j = i+1 : (n_surf + Au_inner_total)
        if norm(all_data(j,5:7)-all_data(i,5:7)) < 2.93
            Au_first = [Au_first; i];
            Au_second = [Au_second; j];
        end
    end
    
end

Au_Au = [Au_first Au_second];

% For Au-S bonding in the hydrophobic part
Au_with_S1_index = [];
S1_index = [];
Au_S1 = [];
for i = 1 : n_surf
    PQ = all_data(i,5:7);
    P = all_data(n_surf+1 : n_surf + Au_inner_total,5:7);
    [k,distance] = dsearchn(P,PQ);
    S1_index = [S1_index; i];
    Au_with_S1_index = [Au_with_S1_index; n_surf+k]; %k is counted from 1, so the real Au_index is n_surf+k
end

Au_S1 = [S1_index  Au_with_S1_index];


% For Au-S bonding in the hydrophilic part
Au_with_S2_index = [];
S2_index = [];
Au_S2 = [];
for i = (n_surf+Au_inner_total+n_alkyl+1) : nu_atom
    PQ = all_data(i,5:7);
    P = all_data(n_surf+1 : n_surf + Au_inner_total,5:7);
    [k,distance] = dsearchn(P,PQ);
    S2_index = [S2_index; i];
    Au_with_S2_index = [Au_with_S2_index; n_surf+k]; %k is counted from 1, so the real Au_index is n_surf+k
end

Au_S2 = [S2_index  Au_with_S2_index];


%% S2 - C1 (PMB)
C1_index = PMB_atom(1:25:end,1); % 1st C1 and every 25th element thereafter
S2_PMB = [ S2_index  C1_index ];

%% writing Janus_PMB.itp
fid = fopen('Janus_PMB.itp','w');
fprintf(fid,';JANUS Particle model\n');
fprintf(fid,'\n');
fprintf(fid,'[ moleculetype ]\n');
fprintf(fid,'; molname  	nrexcl\n'); 
fprintf(fid,'MOL   1\n');
fprintf(fid,'\n');


fprintf(fid,'[ atoms ]\n');
fprintf(fid,';nr  type  resnr  residu   atom    cgnr   charge  mass\n');

s1=0;
s2=0;
s3=0;
s4=0;
s5=0;
for i=1:nu_atom
    % nr label of atom
    fprintf(fid,'%6d',all_data(i,1));
    
    % type of bead B1 to B4
    if all_data(i,3)==1
        fprintf(fid,'%5s',strcat('C5'));
    elseif all_data(i,3)==2
        fprintf(fid,'%5s',strcat('N0'));
    elseif all_data(i,3)==3
        fprintf(fid,'%5s',strcat('C1'));
    elseif all_data(i,3)==4
        %fprintf(fid,'%5s',strcat('Qd'));
        fprintf(fid,'%5s',strcat('N0'));
    elseif all_data(i,3)==5
        fprintf(fid,'%5s',strcat('P4'));
    end
    
    %     fprintf(fid,'%5s',strcat('B',num2str(all_data(i,3))));
    
    % resnr label of molecule
    fprintf(fid,'%5d',1);
    
    % residu name of the molecule
    fprintf(fid,'%5s','MOL');
    
    % atom name of the bead
    %     fprintf(fid,'%5s',strcat('B',num2str(all_data(i,3))));
    if all_data(i,3)==1
        s1=s1+1;
        fprintf(fid,'%5s',strcat('C5'));
    elseif all_data(i,3)==2
        s2=s2+1;
        fprintf(fid,'%5s',strcat('N0'));
    elseif all_data(i,3)==3
        s3=s3+1;
        fprintf(fid,'%5s',strcat('C1'));
    elseif all_data(i,3)==4
        s4=s4+1;
        %fprintf(fid,'%5s',strcat('Qd'));
        fprintf(fid,'%5s',strcat('N0'));
    elseif all_data(i,3)==5
        s5=s5+1;
        fprintf(fid,'%5s',strcat('P4'));
    end
    
    % cgnr charge group number
    fprintf(fid,'%7d',all_data(i,1));
    
    % charge
    fprintf(fid,'%8.4f',all_data(i,4));
    
    % mass
    fprintf(fid,'%8.2f',72.0);
    
    fprintf(fid,'\n');
end
fprintf(fid,'\n');

% Add atoms of PMB here
PMB_a = readtable('PMB_atom.xlsx');
PMB_CG_charge = PMB_a(:,7);
PMB_CG_mass = PMB_a(:,8);

PMB_CG_nr = PMB_atom(:,1);
PMB_CG_type = PMB_CG_type;
PMB_CG_resnr = PMB_atom(:,2);
PMB_atom_name = PMB_CG_type;
PMB_CG_cgnr = PMB_CG_nr;
PMG_all_charge = [];
PMG_all_mass = [];

for i = 1 : num_PMB       
    PMG_all_charge = [PMG_all_charge; PMB_CG_charge];
    PMG_all_mass = [PMG_all_mass ; PMB_CG_mass];
end

% Convert from table data to double data
PMG_all_charge = PMG_all_charge{:,:};
PMG_all_mass = PMG_all_mass{:,:};

for i=1:size(PMB_atom, 1)
    % nr label of atom
    fprintf(fid,'%6d',PMB_CG_nr(i,1));
    
    % type of PMB bead
    fprintf(fid,'%5s',PMB_CG_type(i,1));

    % resnr label of molecule
    fprintf(fid,'%5d',PMB_CG_resnr(i,1));
    
    % residu name of the molecule
    fprintf(fid,'%5s','PMB');
    
    % atom name of the bead
    fprintf(fid,'%5s',PMB_atom_name(i,1));
    
    % cgnr charge group number
    fprintf(fid,'%7d',PMB_CG_cgnr(i,1));
    
    % charge
    fprintf(fid,'%8.4f',PMG_all_charge(i,1));
    
    % mass
    fprintf(fid,'%8.2f',PMG_all_mass(i,1));
    
    fprintf(fid,'\n');
end
fprintf(fid,'\n');


% For bonding

fprintf(fid,'[ bonds ]\n');
fprintf(fid,';  ai    aj    funct    const.\n');

for i=1:m*n_g
    if bond(i,2)==1
        fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',bond(i,3),bond(i,4),1,0.47,1250); % !!!!!__PLEASE ADJUST THE BOND PARAMETERS__!!!!!
    elseif bond(i,2)==2
        fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',bond(i,3),bond(i,4),1,0.445,1250); % !!!!!__PLEASE ADJUST THE BOND PARAMETERS__!!!!!
    end
end

% For Au-Au bonding
for i=1:size(Au_Au,1)
    fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',Au_Au(i,1),Au_Au(i,2),1,0.2885,10000);
end

% for i=1:size(Au_Au,1)
%     fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',Au_Au(i,1),Au_Au(i,2),1,0.577,10000);
% end

% For Au-S1 bonding
for i=1:size(Au_S1,1)
    fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',Au_S1(i,1),Au_S1(i,2),1,0.24,6400);
end

% For Au-S2 bonding
for i=1:size(Au_S2,1)
    fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',Au_S2(i,1),Au_S2(i,2),1,0.24,6400);
end

% For S2-PMB bonding
for i=1:size(S2_PMB,1)
    fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',S2_PMB(i,1), S2_PMB(i,2),1,0.445,1250);
end

fprintf(fid,'\n');

% Bonds of PMBs
PMB_bond = readtable('PMB_bond.xlsx');
bond_first = PMB_bond(:,1);
bond_second = PMB_bond(:,2);
bond_length = PMB_bond(:,4);
bond_const = PMB_bond(:,5);

bond_first = bond_first{:,:};
bond_second = bond_second{:,:};
bond_length = bond_length{:,:};
bond_const = bond_const{:,:};
bond_all = [];
PMB_bond = [];

for i = 1 : num_PMB
    first = num_atom_per_PMB*(i-1) + bond_first + nu_atom;
    second = num_atom_per_PMB*(i-1) + bond_second + nu_atom;
    bond_all = [ first  second ];
    PMB_bond = [PMB_bond; bond_all];
end

bond_all_length = [];
bond_all_const = [];
for i = 1 : num_PMB       
    bond_all_length = [bond_all_length; bond_length];
    bond_all_const = [bond_all_const; bond_const];
end

for i = 1 : size(PMB_bond, 1)
fprintf(fid,'%9d%9d%5d%10.3f%10.3f\n',PMB_bond(i,1),PMB_bond(i,2),1,bond_all_length(i,1),bond_all_const(i,1));
end


% End of bonding force field here

fprintf(fid,'\n');
fprintf(fid,'[ angles ]\n');
fprintf(fid,';  ai aj ak funct\n');
fprintf(fid,'\n');
fprintf(fid,'[ angles ]\n');
fprintf(fid,';  ai aj ak funct\n');

for i=1:(m-1)*n_g
    fprintf(fid,'%9d%9d%9d%5d%8.3f%8.3f\n',angle(i,3),angle(i,4),angle(i,5),2,180,25.0); % !!!!!__PLEASE ADJUST THE LAST TWO ANGLE POTENTIAL PARAMETERS__!!!!!
end
fprintf(fid,'\n');

% Angle of PMBs
PMB_angle = readtable('PMB_angle.xlsx');
angle_first = PMB_angle(:,1);
angle_second = PMB_angle(:,2);
angle_third = PMB_angle(:,3);
angle_value = PMB_angle(:,5);
angle_const = PMB_angle(:,6);

angle_first = angle_first{:,:};
angle_second = angle_second{:,:};
angle_third = angle_third{:,:};
angle_value = angle_value{:,:};
angle_const = angle_const{:,:};

angle_all = [];
PMB_angle = [];

for i = 1 : num_PMB
    first = num_atom_per_PMB*(i-1) + angle_first + nu_atom;
    second = num_atom_per_PMB*(i-1) + angle_second + nu_atom;
    third = num_atom_per_PMB*(i-1) + angle_third + nu_atom;
    angle_all = [ first  second  third ];
    PMB_angle = [PMB_angle; angle_all];
end

angle_all_value = [];
angle_all_const = [];
for i = 1 : num_PMB       
    angle_all_value = [angle_all_value; angle_value];
    angle_all_const = [angle_all_const; angle_const];
end

for i = 1 : size(PMB_angle,1)
    fprintf(fid,'%9d%9d%9d%5d%8.3f%8.3f\n',PMB_angle(i,1),PMB_angle(i,2),PMB_angle(i,3),2,angle_all_value(i,1),angle_all_const(i,1));
end


% Diheral of PMBs

PMB_dihedral = readtable('PMB_dihedrals.xlsx');
dihedral_first = PMB_dihedral(:,1);
dihedral_second = PMB_dihedral(:,2);
dihedral_third = PMB_dihedral(:,3);
dihedral_four = PMB_dihedral(:,4);

dihedral_first = dihedral_first{:,:};
dihedral_second = dihedral_second{:,:};
dihedral_third = dihedral_third{:,:};
dihedral_four = dihedral_four{:,:};

dihedral_all = [];
PMB_dihedral = [];

for i = 1 : num_PMB
    first = num_atom_per_PMB*(i-1) + dihedral_first + nu_atom;
    second = num_atom_per_PMB*(i-1) + dihedral_second + nu_atom;
    third = num_atom_per_PMB*(i-1) + dihedral_third + nu_atom;
    four = num_atom_per_PMB*(i-1) + dihedral_four + nu_atom;
    dihedral_all = [ first  second  third  four ];
    PMB_dihedral = [PMB_dihedral; dihedral_all];
end

fprintf(fid,'\n');
fprintf(fid,'[ dihedrals ]\n');
fprintf(fid,'; ai aj ak al fu xi0 kxi funct\n');

for i = 1 : size(PMB_dihedral,1)
    fprintf(fid,'%9d%9d%9d%9d%5d%8.3f%8.3f%5d\n',PMB_dihedral(i,1),PMB_dihedral(i,2),PMB_dihedral(i,3),PMB_dihedral(i,4),1,0,50,1);
end


% Add constraint of PMB here
PMB_constraint = readtable('PMB_constraint.xlsx');
constraint_first = PMB_constraint(:,1);
constraint_second = PMB_constraint(:,2);
constraint_length = PMB_constraint(:,4);

constraint_first = constraint_first{:,:};
constraint_second = constraint_second{:,:};
constraint_length = constraint_length{:,:};

constraint_all = [];
PMB_constraint = [];

for i = 1 : num_PMB
    first = num_atom_per_PMB*(i-1) + constraint_first + nu_atom;
    second = num_atom_per_PMB*(i-1) + constraint_second + nu_atom;
    constraint_all = [ first  second ];
    PMB_constraint = [PMB_constraint; constraint_all];
end

constraint_all_length = [];
for i = 1 : num_PMB       
    constraint_all_length = [constraint_all_length; constraint_length];
end

fprintf(fid,'[ constraints ]\n');
fprintf(fid,';  i  j 	funct 	length\n');

for i = 1 : size(PMB_constraint, 1)
fprintf(fid,'%9d%9d%5d%10.3f\n',PMB_constraint(i,1),PMB_constraint(i,2),1,constraint_all_length(i,1));
end

fclose(fid);




