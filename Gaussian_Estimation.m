%% Reading the image, and defining blocks and vectors
image=[im2double(cdata)];
[m,n]=size(image);
block = zeros(8,8);
friblock = zeros(8,8);
vectors=zeros((m-8)*(n-8),64);
flag = 0;
index = zeros ((m-8)*(n-8),1);
[mD1_BG,nD1_BG] = size(D1_BG);
[mD1_FG,nD1_FG] = size(D1_FG);
for i=1:m-8
    for j=1:n-8
        block = image(i:i+7,j:j+7);
        friblock=dct2(block);
        flag = flag + 1;
        vectors(flag,:)= toZigzag(friblock);
    end
end
%% Strategy 1&2
%%it is strategy 1 if strategy is 1 and it is strategy 2 if strategy is 0
%strategy = 1







%% %%%%%%%%%%%In This Part we find mu and sigma BG, Fg 
%% Fitting a gaussian function


uhat = zeros (1,64);
sigmahat = zeros(64,64);
%% Background gaussian fit
for k = 1:64
    uhatkb(1,k) = sum(D1_BG(:,k))/mD1_BG;
    nrmkb = D1_BG(:,k) - uhatkb(1,k)*ones(mD1_BG,1);
    
    for h = 1:64
        uhathb(1,h) = sum(D1_BG(:,h))/mD1_BG;
        nrmhb = D1_BG(:,h) - uhathb(1,h)*ones(mD1_BG,1);
        sigmahatb(h,k) = ((transpose(nrmhb)) * nrmkb)/(mD1_BG);
    end
end
%% Foreground gaussian fit
for k = 1:64
    
    uhatkf(1,k) = sum(D1_FG(:,k))/mD1_FG;
    nrmkf = D1_FG(:,k) - uhatkf(1,k)*ones(mD1_FG,1);
    
    for h = 1:64
        uhathf(1,h) = sum(D1_FG(:,h))/mD1_FG;
        nrmhf = D1_FG(:,h) - uhathf(1,h)*ones(mD1_FG,1);
        sigmahatf(h,k) = ((transpose(nrmhf)) * nrmkf)/(mD1_FG);
    end
end
%% foreground sigma and mu
%sigma_FG_D1_ML = zeros (64,64);
sigma_FG_D1_ML=sigmahatf; 
mu_FG_D1_ML = uhatkf;
 
%% background sigma 
%sigma_BG_D1_ML = zeros (64,64);
sigma_BG_D1_ML = sigmahatb; 
mu_BG_D1_ML = uhatkb;







%% %%%%%%%% In this part we find sigma0
%% Defining W0 matrixes
W0_mat_str1 = zeros (64,64);
W0_mat_str2 = zeros (64,64);
for i = 1:64
    W0_mat_str1 (i,i) = W0_str1 (i);
    W0_mat_str2 (i,i) = W0_str2 (i);
end

%% Sigma0
Sigma0_str1 = zeros(64,64);
Sigma0_str2 = zeros(64,64);

for counter = 1 : length(alpha)
    
    Sigma0_str1 = alpha(counter) * W0_mat_str1;
    Sigma0_str2 = alpha(counter) * W0_mat_str2;
  
    mu_D1_BG = Sigma0_str1 * inv (Sigma0_str1 +  sigma_BG_D1_ML ./ length(D1_BG) ) * (transpose(sum(D1_BG)) ./length(D1_BG)) + (( sigma_BG_D1_ML ./ length(D1_BG)) *  inv (Sigma0_str1 +  sigma_BG_D1_ML ./ length(D1_BG) ) * transpose(mu0_BG_str1));
    sigma_D1_BG = Sigma0_str1 * inv (Sigma0_str1 +  sigma_BG_D1_ML ./ length(D1_BG) ) * sigma_BG_D1_ML ./ length(D1_BG);
    mu_D1_FG = Sigma0_str1 * inv (Sigma0_str1 +  sigma_FG_D1_ML ./ length(D1_FG) ) * (transpose(sum(D1_FG)) ./length(D1_FG)) + (( sigma_FG_D1_ML ./ length(D1_FG)) *  inv (Sigma0_str1 +  sigma_FG_D1_ML ./ length(D1_FG) ) * transpose(mu0_FG_str1));
    sigma_D1_FG = Sigma0_str1 * inv (Sigma0_str1 +  sigma_FG_D1_ML ./ length(D1_FG) ) * sigma_FG_D1_ML ./ length(D1_FG);
    
    MU_predictive_D1_BG = mu_D1_BG;
    MU_predictive_D1_FG = mu_D1_FG;
    
    
    SIGMA_predictive_BG_D1 = sigma_BG_D1_ML + sigma_D1_BG;
    SIGMA_predictive_FG_D1 = sigma_FG_D1_ML + sigma_D1_FG;
    
    
    
    
    
end






