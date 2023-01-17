%% Distributed Conjugate Gradient Algorithm for Solving Linear Equations
% showing the effectiveness of DCG using network localization as an example
clc;clear;close all;
currentFolder = pwd;
addpath(genpath(currentFolder));
%% network settings
edgenoise = 0.00;  
sradius = 0.8; 
npoints = 10;     
Boxscale = 100;    
randstate = 6; 
nanchors = 3;
anchor_initialization = 'triangle';
[points,anchors,non_anchors,distMatrix,ConnectivityM,neighbors]=generate_network_2d...
    (nanchors, randstate,edgenoise,npoints,sradius,anchor_initialization);
plot_gg(points,anchors,non_anchors,[],ConnectivityM,Boxscale,'network topo')
%% transforming network localization problem to linear equations
% see the paper 'A barycentric coordinate based distributed localization algorithm for sensor networks', 
% IEEE Transactions on Signal Processing 62 (2014) 4760–4771
pos_a = points(1:nanchors,:);
AA = zeros(npoints,npoints);
for i = 1:nanchors
    AA(i,:) = [zeros(1,i-1) 1 zeros(1,npoints-i)];
end
for i = nanchors+1:npoints
    this_A = zeros(1,npoints);
    this_neighbors_num = neighbors(i,1);
    if this_neighbors_num <3
        AA(i,i) = 1;
        continue
    end
    this_neighbors = neighbors(i,2:1+this_neighbors_num);
    if size(this_neighbors,2)<3
        AA(i,i) = 1;
        continue
    end
    combines = nchoosek(this_neighbors,3);
    num_combines = size(combines,1);
    m = 0;
    for j = 1:num_combines
        c = combines(j,:);
        if ConnectivityM(c(1),c(2))==1 && ConnectivityM(c(2),c(3))==1 ...
                && ConnectivityM(c(1),c(3))==1
            
            m = m + 1;
            temp_A = F_bc(distMatrix(i,c(1)),distMatrix(i,c(2)),distMatrix(i,c(3)),...
                distMatrix(c(1),c(2)),distMatrix(c(1),c(3)),distMatrix(c(2),c(3)));
            this_A(c(1)) = this_A(c(1))+temp_A(1);
            this_A(c(2)) = this_A(c(2))+temp_A(2);
            this_A(c(3)) = this_A(c(3))+temp_A(3);
        end
    end
    if m>=1
        this_A = this_A/m;
        AA(i,:) = this_A;
    else
        AA(i,i) = 1;
    end
    
end
G_AA = AA~= 0;
ConnectivityM_GG = zeros(npoints,npoints);
for i=1:npoints
    cnum=0;
    for j=1:npoints
        if G_AA(i,j)==1
            cnum = cnum+1;
            ConnectivityM_GG(i,j) = 1;
        end
    end
end
BB = AA(non_anchors,anchors);
CC = AA(non_anchors,non_anchors);
MM = eye(size(non_anchors,2)) - CC;
n_free_agent = size(non_anchors,2);
%% DCG
max_iteration = 3e5;
max_diff = 1e-11;
max_ind = npoints;
pos_iter = zeros(npoints-nanchors,2,max_ind+1);
pos_s = 0*ones(npoints-nanchors,2);
ind = 1;
pos_iter(:,:,ind) = pos_s;
error_lst_dcg = zeros(npoints-nanchors,max_iteration);

N = cell(npoints,1);
nN = neighbors(:,1);
for i = 1:npoints
    N{i} = neighbors(i,2:1+nN(i));
end
% initialize
b = zeros(n_free_agent,2);
bb = zeros(n_free_agent,2);
d = zeros(n_free_agent,2);
for i = 1:n_free_agent
    b(i,:) = sum(BB(i,1:nanchors)*pos_a,1);
end
for i = 1:n_free_agent
    bb(i,:) = [sum(MM(:,i)'*b(:,1)) sum(MM(:,i)'*b(:,2))];
end
AA = MM'*MM;

bx = bb(:,1);
dx = zeros(n_free_agent,1);
rx = -bx;
pos_x = zeros(n_free_agent,1);

by = bb(:,2);
dy = zeros(n_free_agent,1);
ry = -by;
pos_y = zeros(n_free_agent,1);

rr = -bb;
nN = nN(non_anchors);
N = N(non_anchors);
for i = 1:n_free_agent
    N{i} = intersect(N{i},non_anchors);
    nN(i) = size(N{i},2);
end
denom_x = zeros(n_free_agent,1);
denom_y = zeros(n_free_agent,1);
num_x = zeros(n_free_agent,1);
num_y = zeros(n_free_agent,1);
num_x2 = zeros(n_free_agent,1);
num_y2 = zeros(n_free_agent,1);
alpha_x = zeros(n_free_agent,1);
alpha_y = zeros(n_free_agent,1);
temp_vect = zeros(n_free_agent,2);

for iter = 1:max_iteration
    rr_old = rr;
    rr_dis = cell(n_free_agent,1);
    for i = 1:n_free_agent
        rr_dis{i} = [i rr_old(i,:)];
    end
    for ii = 1:n_free_agent 
        rr_dis_pre = rr_dis;
        for i = 1:n_free_agent
            for j = 1:n_free_agent
                if ConnectivityM(non_anchors(i),non_anchors(j))==1
                    rr_dis{i} = [rr_dis{i};rr_dis_pre{j}];
                end
            end
            rr_dis{i} = unique(rr_dis{i},'rows');
        end
    end
    for i = 1:n_free_agent
        deno_x = rr_dis{i}(:,2);
        deno_y = rr_dis{i}(:,3);
        denom_x(i) = sum(deno_x.*deno_x);
        denom_y(i) = sum(deno_y.*deno_y);
        rr(i,1) = AA(i,:)*pos_x;
        rr(i,1) = rr(i,1) - bx(i);
        rr(i,2) = AA(i,:)*pos_y;
        rr(i,2) = rr(i,2) - by(i);
    end
    rr_dis = cell(n_free_agent,1); 
    for i = 1:n_free_agent
        rr_dis{i} = [i rr(i,:)];
    end
    for ii = 1:n_free_agent 
        rr_dis_pre = rr_dis;
        for i = 1:n_free_agent
            for j = 1:n_free_agent
                if ConnectivityM(non_anchors(i),non_anchors(j))==1
                    rr_dis{i} = [rr_dis{i};rr_dis_pre{j}];
                end
            end
            rr_dis{i} = unique(rr_dis{i},'rows');
        end
    end
    for i = 1:n_free_agent
        num_x(i) = inner_product(rr_dis{i}(:,2),rr_dis{i}(:,2),n_free_agent);
        num_y(i) = inner_product(rr_dis{i}(:,3),rr_dis{i}(:,3),n_free_agent);
    end
    
    if max(num_x)<max_diff
        %['DCG meet the min resudial ' num2str(3*n_free_agent*iter)]
        break
    end
    dd_dis = cell(n_free_agent,1); 
    for i = 1:n_free_agent
        dx(i) = -rr(i,1)+num_x(i)*dx(i)/denom_x(i);
        dy(i) = -rr(i,2)+num_y(i)*dy(i)/denom_y(i);
        dd_dis{i} = [i dx(i) dy(i)];
    end
    for ii = 1:n_free_agent 
        dd_dis_pre = dd_dis;
        for i = 1:n_free_agent 
            for j = 1:n_free_agent
                if ConnectivityM(non_anchors(i),non_anchors(j))==1
                    dd_dis{i} = [dd_dis{i};dd_dis_pre{j}];
                end
            end
            dd_dis{i} = unique(dd_dis{i},'rows');
        end
    end
    for i = 1:n_free_agent
        num_x2(i) = inner_product(dd_dis{i}(:,2),rr_dis{i}(:,2),n_free_agent);
        num_y2(i) = inner_product(dd_dis{i}(:,3),rr_dis{i}(:,3),n_free_agent);
        temp_vect(i,1) = inner_product(AA(i,:),dd_dis{i}(:,2),n_free_agent);
        temp_vect(i,2) = inner_product(AA(i,:),dd_dis{i}(:,3),n_free_agent);
    end
    temp_vect_dis = cell(n_free_agent,1);
    for i = 1:n_free_agent
        temp_vect_dis{i} = [i temp_vect(i,:)];
    end
    for ii = 1:n_free_agent 
        temp_vect_dis_pre = temp_vect_dis;
        for i = 1:n_free_agent 
            for j = 1:n_free_agent
                if ConnectivityM(non_anchors(i),non_anchors(j))==1
                    temp_vect_dis{i} = [temp_vect_dis{i};temp_vect_dis_pre{j}];
                end
            end
            temp_vect_dis{i} = unique(temp_vect_dis{i},'rows');
        end
    end    
    for i = 1:n_free_agent
        denom_2x = inner_product(dd_dis{i}(:,2),temp_vect(:,1),n_free_agent);
        denom_2y = inner_product(dd_dis{i}(:,3),temp_vect(:,2),n_free_agent);
        alpha_x(i) = -num_x2(i)/denom_2x;
        alpha_y(i) = -num_y2(i)/denom_2y;
        pos_x(i) = pos_x(i)+alpha_x(i)*dx(i);
        pos_y(i) = pos_y(i)+alpha_y(i)*dy(i);
    end
    pos_loc = [pos_x pos_y];
    for j = 1:n_free_agent
        error_lst_dcg(non_anchors(j)-3,iter) = F_error_2d(1,pos_loc(j,:),points(non_anchors(j),:));
    end
    error = F_error_2d(n_free_agent,pos_loc,points(non_anchors,:));
    pos_iter(non_anchors-3,:,iter+1) = pos_loc;  
end

%% plot traces
pos_lst = cell(npoints-nanchors,1);
pos_iter = pos_iter(:,:,1:iter+1);
for i = 1:npoints-nanchors
    aaa = pos_iter(i,:,:);
    rrr = reshape(aaa,2,iter+1);
    rrr = rrr';
    rrr = rrr(1:iter,:);
    pos_lst{i} = rrr;
end
plot_trace(points,ConnectivityM_GG,pos_lst,Boxscale,non_anchors)












