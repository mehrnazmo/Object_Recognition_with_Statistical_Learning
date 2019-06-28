
%% This scripts classifies an image of a cheetah to 2 different classes, Cheetah and Grass

%% Reading the image, and defining blocks and vectors
image=[im2double(cdata)];
[m,n]=size(image);


    
%% Fitting a gaussian function
%BG = zeros (1,64);
uhat = zeros (1,64);
%nrm = zeros (1053,1);
sigmahat = zeros(64,64);


%% Background gaussian fit
for k = 1:64
    %BG(k) = (TrainsampleDCT_BG(:,k));
    
    uhatkb(1,k) = sum(TrainsampleDCT_BG(:,k))/1053;
    nrmkb = TrainsampleDCT_BG(:,k) - uhatkb(1,k)*ones(1053,1);
    
    for h = 1:64
        uhathb(1,h) = sum(TrainsampleDCT_BG(:,h))/1053;
        nrmhb = TrainsampleDCT_BG(:,h) - uhathb(1,h)*ones(1053,1);
        
        sigmahatb(h,k) = ((transpose(nrmhb)) * nrmkb)/(1053);
    end
end
%% Foreground gaussian fit
for k = 1:64
    %BG(k) = (TrainsampleDCT_BG(:,k));
    
    uhatkf(1,k) = sum(TrainsampleDCT_FG(:,k))/250;
    nrmkf = TrainsampleDCT_FG(:,k) - uhatkf(1,k)*ones(250,1);
    
    for h = 1:64
        uhathf(1,h) = sum(TrainsampleDCT_FG(:,h))/250;
        nrmhf = TrainsampleDCT_FG(:,h) - uhathf(1,h)*ones(250,1);
        
        sigmahatf(h,k) = ((transpose(nrmhf)) * nrmkf)/(250);
    end
end
    
 
 %% Plotting gaussian background/foreground
 
 %% foreground sigma
 sf = zeros(64,1);
 for i = 1:64 
     sf(i)=sigmahatf(i,i); 
 end
%% background sigma 
sb = zeros(64,1);
for i = 1:64 
    sb(i)=sigmahatb(i,i); 
end

figure;
 %% Background
 for i = 1:64
     mub = uhathb(i);
     muf = uhathf(i);
     x = min  ((-3*sqrt(sb(i)) + mub) , (-3*sqrt(sf(i)) + muf))   : 0.001*min(sqrt(sb(i)) , sqrt(sf(i))) :  max(3*sqrt(sb(i)) + mub, 3*sqrt(sf(i)) + muf);
     yb = normpdf(x,mub,sqrt(sb(i)));
     yf = normpdf(x,muf,sqrt(sf(i)));
     subplot(8,8,i);
     plot(x,yb,'r');
     title(['Index number',num2str(i)]);
     hold on;
     plot(x,yf,'--');
         
 end
 
 
 best = [32,27,25,40,33,18,41,45];
 worst = [64,63,59,3,5,62,4,60];
 figure;
 for j = 1:8
 
 indb = best(j);
 mub = uhathb(indb);
 muf = uhathf(indb);
 x = min  ((-3*sqrt(sb(indb)) + mub) , (-3*sqrt(sf(indb)) + muf))   : 0.001*min(sqrt(sb(indb)) , sqrt(sf(indb))) :  max(3*sqrt(sb(indb)) + mub, 3*sqrt(sf(indb)) + muf);
 yb = normpdf(x,mub,sqrt(sb(indb)));
 yf = normpdf(x,muf,sqrt(sf(indb)));
 subplot(2,4,j)
     plot(x,yb,'r');
     title(['Index number',num2str(indb)]);
     hold on;
     plot(x,yf,'--');
 end
     
 figure;
 for j = 1:8
 indb = worst(j);
 mub = uhathb(indb);
 muf = uhathf(indb);
 x = min  ((-3*sqrt(sb(indb)) + mub) , (-3*sqrt(sf(indb)) + muf))   : 0.001*min(sqrt(sb(indb)) , sqrt(sf(indb))) :  max(3*sqrt(sb(indb)) + mub, 3*sqrt(sf(indb)) + muf);
 yb = normpdf(x,mub,sqrt(sb(indb)));
 yf = normpdf(x,muf,sqrt(sf(indb)));
 subplot(2,4,j)
     plot(x,yb,'r');
     title(['Index number',num2str(indb)]);
     hold on;
     plot(x,yf,'--');
 end
 
 figure;
 x = min  ((-3*sqrt(sb(1)) + mub) , (-3*sqrt(sf(1)) + muf))   : 0.001*min(sqrt(sb(1)) , sqrt(sf(1))) :  max(3*sqrt(sb(1)) + mub, 3*sqrt(sf(1)) + muf);
 yb = normpdf(x,mub,sqrt(sb(1)));
 yf = normpdf(x,muf,sqrt(sf(1)));
 plot(x,yb,'r');
 hold on;
 plot(x,yf,'--');

 
 %% Part c
% 64-dimensional

A_best = zeros(size(image,1),size(image,2));
block = zeros(8,8);
friblock = zeros(8,8);
vectors=zeros((m-8)*(n-8),64);
flag = 0;
index = zeros ((m-8)*(n-8),1);
d = [];
alpha = [];
flagd=zeros(m-8,n-8);
mu_BG_best = uhatkb(best);
sigma_BG_best = sigmahatb(best,best);
mu_FG_best = uhatkf(best);
sigma_FG_best = sigmahatf(best,best);
for i=1:m-8
    for j=1:n-8
        
        block = image(i:i+7,j:j+7);
        friblock=dct2(block);
        flag = flag + 1;
        vectors(flag,:)= toZigzag(friblock);
        d_G = (vectors(flag,:) - uhatkb) * inv(sigmahatb) * transpose(vectors(flag,:) - uhatkb);
        alpha_G = log (((2*pi).^64)* det(sigmahatb)) - 2*log(0.8081);
        d_C = (vectors(flag,:) - uhatkf) * inv(sigmahatf) * transpose(vectors(flag,:) - uhatkf);
        alpha_C = log (((2*pi).^64)* det(sigmahatf)) - 2*log(0.8081);
        if ((d_G+alpha_G) <= (d_C+alpha_C))
            flagd(i,j) = 0;
        else
            flagd(i,j) = 1;
        end
        
        zigzagform_best = vectors(flag,best);
        d_BG_best = ((zigzagform_best-mu_BG_best))*inv(sigma_BG_best)*((zigzagform_best-mu_BG_best).');
        alpha_BG_best = log(((2*pi)^64)*det(sigma_BG_best))-2*log(0.8081);
        d_FG_best = ((zigzagform_best-mu_FG_best))*inv(sigma_FG_best)*((zigzagform_best-mu_FG_best).');
        alpha_FG_best = log(((2*pi)^64)*det(sigma_FG_best))-2*log(0.1919);
        if(d_BG_best+alpha_BG_best  <= d_FG_best+alpha_FG_best )
            A_best(i,j) = 0;
        else
            A_best(i,j) = 1;
        end
            
    end
end

figure();
imagesc(flagd);
%colormap(gray(255));
%image_read_mask = imread('cheetah_mask.bmp');
image_read_mask = im2double(mask);
[p,q]=size(image_read_mask);
flag_det_rate = 0;
false_alarm = 0;
for counter1 = 1:size(image,1)-8
    for counter2 = 1:size(image,2)-8
        if(image_read_mask(counter1,counter2)==1 && flagd(counter1,counter2)==1)
            flag_det_rate = flag_det_rate+1; 
        end
        if(image_read_mask(counter1,counter2)==0 && flagd(counter1,counter2)==1)
            false_alarm = false_alarm + 1;
        end
    end
end
det_rate = flag_det_rate/sum(sum(image_read_mask,1),2);
false_error_rate = false_alarm/((p*q)-sum(sum(image_read_mask,1),2));
p_error = (false_error_rate*0.8081)+((1-det_rate)*0.1919);
fprintf('Detection Rate for Cheetah: %f \n',flag_det_rate/sum(sum(image_read_mask,1),2));
fprintf('False Alarm Rate: %f \n',false_alarm/((p*q)-sum(sum(image_read_mask,1),2)));
fprintf('Probability of Error for 64 dimensional: %f \n',p_error);



figure();
imagesc(A_best);
%colormap(gray(255));
flag_det_rate = 0;
flag_false_error_rate = 0;
for counter1 = 1:size(image,1)-8
    for counter2 = 1:size(image,2)-8
        if(image_read_mask(counter1,counter2)==1 && A_best(counter1,counter2)==1)
            flag_det_rate = flag_det_rate+1;
        end
        if(image_read_mask(counter1,counter2)==0 && A_best(counter1,counter2)==1)
            flag_false_error_rate = flag_false_error_rate+1;
        end
    end
end
det_rate = flag_det_rate/sum(sum(image_read_mask,1),2);
false_error_rate = flag_false_error_rate/((p*q)-sum(sum(image_read_mask,1),2));
p_error = (false_error_rate*0.8081)+((1-det_rate)*0.1919);
fprintf('Detection Rate for Cheetah for 8 dimensional: %f \n',det_rate);
fprintf('False Error Rate for Background 8 dimensional: %f \n',false_error_rate);
fprintf('Probability of Error 8 dimensional: %f \n',p_error);


%{




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

%% Calculating 0 and 1 for each point of the image showing if it is background or foreground; resulting in showing foreflag as the final image mask
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

%% Displaying the final result
imagesc(foreflag);
%% Calculating Error
A=xor(foreflag,mask);
error = sum(sum(A,1),2)/(m*n);

%}

 