function [points,anchors,non_anchors,distMatrix,ConnectivityM,neighbors]=generate_network_2d...
    (nanchors, randstate,edgenoise,npoints,sradius,anchor_distribution)
%output
% points: n*2£¬ground truth locations
% anchors: anum*1 anchor IDs
% non_anchors:(n-anum)*1 non anchors\ IDs
% distMatrix£ºn*n distance matrix
% ConnectivityM£ºconnectivity matrix

%input
%edgenoise:     %noise Level
%npoints        %node number
%sradius        %reception radius
%%
nf1=edgenoise;         
%nf2=anglenoise;
radius=sradius;      %distance of neighbor nodes
Boxscale = 100;  
dim=2;           %dimensionality
%% initialization
rand('seed',randstate);
points=rand(dim,npoints)-0.5;
points = points';
Y=pdist(points);
trueDistMatrix=squareform(Y);
distMatrix=zeros(npoints, npoints);
ConnectivityM=zeros(npoints,npoints);
for i=1:npoints
    for j=1:npoints
        if trueDistMatrix(i,j)>radius
            distMatrix(i,j)=NaN;
            ConnectivityM(i,j)=0;
        elseif i==j 
            distMatrix(i,j)=0;
            ConnectivityM(i,j)=0;
        else
            distMatrix(i,j)=1;
            ConnectivityM(i,j)=1;
        end
    end
end

% Make distance matrix symmetry
for i=1:npoints
    for j=1:npoints
        if ConnectivityM(i,j)==1
            dis=(1+nf1*(rand()-0.5))*trueDistMatrix(i,j);
            %dis=trueDistMatrix(i,j)+nf1*randn();
            if dis>0
                distMatrix(i,j)=dis;
            else
                distMatrix(i,j)=trueDistMatrix(i,j);
            end
        end
    end
end
neighbors = zeros(npoints,1);
for i=1:npoints
    cnum=0;
    for j=1:npoints
        if ConnectivityM(i,j)==1
            cnum = cnum+1;
            neighbors(i,cnum+1)=j;
        end
    end
    neighbors(i,1)=cnum;
end
%enum=(sum(DisNeighbor(:,1))+sum(AngleNeighbor(:,1))+sum(BothNeighbor(:,1)))/2;
if strcmp(anchor_distribution,'triangle')
    anchors(1) = 1;
    anchors(2) = neighbors(1,2);
    comm_nei = intersect(neighbors(1,2:end),neighbors(anchors(2),2:end));
    comm_nei = comm_nei(comm_nei>0);
    anchors(3) = comm_nei(1);
    points([2 anchors(2)],:) = points([anchors(2) 2],:);
    points([3 anchors(3)],:) = points([anchors(3) 3],:);
    Y=pdist(points);
    trueDistMatrix=squareform(Y);
    distMatrix=zeros(npoints, npoints);
    ConnectivityM=zeros(npoints,npoints);
    for i=1:npoints
        for j=1:npoints
            if trueDistMatrix(i,j)>radius
                distMatrix(i,j)=NaN;
                ConnectivityM(i,j)=0;
            elseif i==j 
                distMatrix(i,j)=0;
                ConnectivityM(i,j)=0;
            else
                distMatrix(i,j)=1;
                ConnectivityM(i,j)=1;
            end
        end
    end    
    % Make distance matrix symmetry
    for i=1:npoints
        for j=1:npoints
            if ConnectivityM(i,j)==1
                dis=(1+nf1*(rand()-0.5))*trueDistMatrix(i,j);
                %dis=trueDistMatrix(i,j)+nf1*randn();
                if dis>0
                    distMatrix(i,j)=dis;
                else
                    distMatrix(i,j)=trueDistMatrix(i,j);
                end
            end
        end
    end
    neighbors = zeros(npoints,1);
    for i=1:npoints
        cnum=0;
        for j=1:npoints
            if ConnectivityM(i,j)==1
                cnum = cnum+1;
                neighbors(i,cnum+1)=j;
            end
        end
        neighbors(i,1)=cnum;
    end
end
anchors = 1:nanchors;
non_anchors = nanchors+1:npoints;

end