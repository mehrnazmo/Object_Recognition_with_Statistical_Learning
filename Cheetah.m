%%
%data=imread('cheetah.bmp');
%%clear;
%%close all;
%%clc;
image=[im2double(cdata)];
%image=[image,zeros(256,2)]; 
%image is now double with dimensions 256*272
[m,n]=size(image);
%gimage=rgb2gray(image);
%converting image to 8*8 cells
%imcell = mat2cell (image,repmat(8,[1,32]),repmat(8,[1,34]));

%frimage=cell(32,34);



block = zeros(8,8);
friblock = zeros(8,8);
vectors=zeros((m-8)*(n-8),64);
flag = 0;
index = zeros ((m-8)*(n-8),1);

for i=1:m-8
    for j=1:n-8
        block = image(i:i+7,j:j+7);
        friblock=dct2(block);
        flag = flag + 1;
        vectors(flag,:)= toZigzag(friblock);
        %{
        ind = reshape(1:numel(friblock), size(friblock));         %# indices of elements
        ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
        ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
        ind(ind==0) = [];   %# keep non-zero indices
        
        
        vectors(flag,:)=friblock(ind);
        %}
    end
end

    

%% finding second index of each vector (corresponding to each point)
      
for k=1:(m-8)*(n-8)
    v = vectors(k,:);
    v_sorted = sort(v, 'descend');
    second_max = v_sorted(2);
    indexx = min(find (v == second_max));  %%several places can have 2nd greatest values
    index(k,1)=indexx;
end
hist(index);

%% finding second index of each vector of Foregraound (Cheetah)

indexFG = zeros (250,1);
for k=1:250
    v = TrainsampleDCT_FG(k,:);
    v_sorted = sort(v, 'descend');
    second_max = v_sorted(2);
    indexf = (find (v == second_max));  %%several places can have 2nd greatest values
    indexFG(k,1)=indexf;
end
hist(indexFG);

%% finding second index of each vector of Backgraound (Grass)

indexBG = zeros (1053,1);
for k=1:1053
    v = TrainsampleDCT_BG(k,:);
    v_sorted = sort(v, 'descend');
    second_max = v_sorted(2);
    indexb = (find (v == second_max));  %%several places can have 2nd greatest values
    indexBG(k,1)=indexb;
end
hist(indexBG);



%% convertng background indexes to a vector containing indices and the probability each index occurs

a = unique(indexBG);
bout = [a,histc(indexBG(:),a)];
bout(:,2) = bout(:,2)/sum(bout(:,2));


%% convertng foreground indexes to a vector containing indices and the probability each index occurs

b = unique(indexFG);
fout = [b,histc(indexFG(:),b)];
fout(:,2) = fout(:,2)/sum(fout(:,2));


%% convertng image indexes to a vector containing indices and the probability each index occurs

c = unique(index);
iout = [c,histc(index(:),c)];
iout(:,2) = iout(:,2)/sum(iout(:,2));



foreflag = zeros(m,n);
  




for r=1:m-8
    for q=1:n-8
        p=(n-8)*(r-1)+q;
    ipos = find(iout(:,1)==index(p,1));  %ipos is the location of the index in iout
    Pimage = iout(ipos,2);
    
    if (any(fout(:,1)== index(p,1)))
        if (any(bout(:,1)==index(p,1)))
        gpos= find(bout(:,1)==iout(ipos));
        Pgrass = bout(gpos,2)*0.8081;
        cpos= find(fout(:,1)==iout(ipos));
        Pcheetah = fout(cpos,2)*0.1919;
            if Pcheetah >= Pgrass
            foreflag(r,q) = 1;
            else
            foreflag(r,q) = 0;
            end
        else
            foreflag(r,q) = 1;
        end
    elseif (any(bout(:,1)==iout(ipos)))
        foreflag(r,q) = 0;    
    elseif (~any(bout(:,1)==iout(ipos)))
        [ttb,lb] = min(abs(bout(:,1)-iout(ipos)));
        [ttf,lf] = min(abs(fout(:,1)-iout(ipos)));
        if abs(ttb)<=abs(ttf)
            foreflag(r,q) = 0;
        else
            foreflag(r,q) = 1;  
        end  
    end
    end
end


%result = invzigzag(foreflag,m-8,n-8);
%f=mat2gray(foreflag);
imagesc(foreflag);
%imshow(foreflag,[0,255]);

%colormap(gray(255));

A=xor(foreflag,mask);
error = sum(sum(A,1),2)/(m*n);



        


        
        
%%  
%{
        
%vectors=zeros(32*34,64);

%vectors_sorted=zeros(32*34,64);

%%computing DCT of each cell, the result is called frimage which is a cell
%%of matrixes
for j=1:34
    for i=1:32
        frimage{i,j}=dct(imcell{i,j});
    end
end

%%Converting each 8*8 matrix of this cell to a vector of 64*1
flag=0;
for j=1:34
    for i=1:32
        ind = reshape(1:numel(frimage{i,j}), size(frimage{i,j}));         %# indices of elements
        ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
        ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
        ind(ind==0) = [];                           %# keep non-zero indices
        flag = flag + 1;
        vectors(flag,:)=frimage{i,j}(ind);                        %# get elements in zigzag order
    end
end


%%Finding the indices of second largest values of each vector (each row) 
index = zeros (32*34,1);
for k=1:34*32
    v = vectors(k,:);
    v_sorted = sort(v, 'descend');
    second_max = v_sorted(2);
    indexx = min(find (v == second_max));  %%several places can have 2nd greatest values
    index(k,1)=indexx;
end
hist(index);




%%indexBG = indexBG * 0.8081/sum(indexBG);
%%indexFG = indexFG * 0.1919/sum(indexFG);
%%index2 = index/sum(index);

%}



 