function [num_mol, x_pos, y_pos]=peakfinder_HW(Img, threshold, min_I_peak, mode, dist)
%% Background subtraction
%threshold=median(median(Img));
%Img(Img<threshold)=threshold;

%Img_bksub=uint16(zeros(512,256));
%[TG] img-wiener2(image,[m,n]) : [m n] specifies the size (m-by-n) of the neighborhood used to estimate the local image mean and standard deviation
Img_bksub=uint16(Img-wiener2(Img, [25 25]));
%Img_bksub=Img-medfilt2(Img, [25,25]); % Median filtering

%[TG] Gausian filter size - ?! how to define filter size of Gausian filter
Img_bksub2=imgaussfilt(Img_bksub,0.4,'FilterSize',[5 5]); % sigma=0.5 & 5x5 window gaussian filtering
%H = fspecial('average',[3,3]); %average filter
%Img_bksub2 = Img-imfilter(Img,H);

threshold=mean(mean(Img_bksub2));
Img_bksub2(Img_bksub2<threshold)=threshold; % Low intensity pixel removal by threshold
%Img_bksub2(Img_bksub2<threshold)=0;
%Img_bksub3=imhmax(Img_bksub2,threshold); % Low intensity pixel removal by threshold
%Img_bksub2=Img;

%% Finding regional max
%[TG] imregionalmax(I,conn) : conn - pixel connectivity, 
%[TG] if conn = 8, pixels are connected if their edges or corners touch. The neighborhood of a pixel are the adjacent pixels in the horizontal, vertical, or diagonal direction.
loc_max=imregionalmax(Img_bksub2,8);

%Image crop
%[TG] crop image - image 잘라내기, 보고싶지 않은 부분은 0으로 만듬.
cropmask=true(size(Img));
cropmask(10:end-10,10:end-10)=false;
loc_max(cropmask)=false;


% get x/y-position of selected regional max
[y_pos, x_pos]=find(loc_max); % return [row, col] value
%[TG] y_pos 가 가지고 있는 요소 갯수, 예를들어 y_pos = [1,2,3] 이면 numel(y_pos) = 3
num_mol=numel(y_pos);

%% Verification by centroid
selc_mol=false(num_mol,1);
PSF=zeros(5);

for i=1:num_mol
    PSF=Img_bksub2(y_pos(i)+(-2:2), x_pos(i)+(-2:2)); % generate 5x5 PSF
    if PSF(3,3) > min_I_peak
        CTR = centroid(PSF); % return CTR=[y,x] of centroid by Dr MJS
        if all([all(CTR>2.5) all(CTR<3.5)])
            selc_mol(i)=true;
            x_pos(i)=x_pos(i)+CTR(2)-3;
            y_pos(i)=y_pos(i)+CTR(1)-3;
        end
    end
end

x_pos = x_pos(selc_mol); % Position of verified molecules
y_pos = y_pos(selc_mol);

%% Neighboorhood removal by dist criteria
if mode==1
    [x_pos, y_pos]=near_remove(x_pos, y_pos, dist); 
end
num_mol=length(x_pos); % Number of verified molecules
end
        
    









