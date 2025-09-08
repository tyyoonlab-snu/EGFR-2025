function [cond_img, img_data_collection]=condensate_analyzer_confocal_FUNCTION_reanalysis(cell_selc,background_setting, peak_background_fold,tolerance_parameter,fig_num,ch_name, average_int_cond_thres)
%% Performing FFT with the image
%cell_selc = img_ch_green ;
cell_selc = double(cell_selc) ;
Img_FFT = fft2(cell_selc) ; % Fourier transform image
Img_Fsh = fftshift(Img_FFT) ;  %Centered Fourier transform image
log_Fsh = log(1+Img_Fsh);
f2=figure(2);
app.f2=f2;

fig_num_f2 = fig_num*4 ;
subplot(4,2,1+fig_num_f2) ;
imagesc(abs(log_Fsh));  %logarimathic Fourier transform image

% masking with circle
center_pos = fliplr(size(cell_selc)/2); %CELL_v1 updated
circle_radius = min(size(cell_selc))/3; %CELL_v1 updated
 
circ = drawcircle('Center',center_pos,'Radius',circle_radius);
title(['FFT image in ',ch_name]);

circ_mask = createMask(circ);

%
%(made by Dr.G)blurring filter of circular mask to prevent from occurence of diffraction
%pattern on inverseFFT-image
%windowSize : Filter size: we can change this parameter 
windowSize = 100 ;
kernel = ones(windowSize) / windowSize ^ 2;
circ_mask = conv2(single(circ_mask), kernel, 'same');
 
subplot(4,2,2+fig_num_f2) ;
imshow(log_Fsh.*circ_mask,[])
title(['FFT MASK in ',ch_name]);

filt_Fsh= Img_Fsh.*double(circ_mask);

% Performing Inverse FFT with masked image
filt_FFT= ifftshift(filt_Fsh);
filt_img = real(ifft2(filt_FFT)) ;
subplot(4,2,3+fig_num_f2) ;
imshow(cell_selc);
colormap(jet)
title(['Before FFT-inverseFFT in ',ch_name]);

subplot(4,2,4+fig_num_f2) ;
imshow(filt_img, [10,200]);
colormap(jet)
title(['After FFT-inverseFFT in ',ch_name]);

%%
img_vec = filt_img(cell_selc~=0);
[mu,std]= normfit(img_vec); %normal distribution assumption.
background = mu+background_setting*std ; % Background can be defined by average + 1.5*std
filt_dat = filt_img-background;

filt_dat(filt_dat<0) = 0 ;
filt_ratio = filt_img./background ; %230917 updated
%
mask_surround = zeros(size(filt_dat));
mask_surround(10:end-10, 10:end-10) = 1 ;
filt_dat_exclude_surround = filt_dat.*mask_surround;


locmax = imregionalmax(filt_dat_exclude_surround,8) ;  %230917 updated
[centY_temp, centX_temp] = find(locmax); %230917 updated
num_peak_temp = numel(centX_temp) ; %230917 updated

% 2230917 updated start
ppc = [];
for kk = 1:num_peak_temp
    peak_partition_coefficient = filt_ratio(centY_temp(kk),centX_temp(kk)); 
    ppc = [ppc,peak_partition_coefficient] ;

end

selc_peak = false(num_peak_temp,1);
for kk = 1:num_peak_temp
    if ppc(kk)> peak_background_fold
        selc_peak(kk) =true;
    end
    
end

centX = centX_temp(selc_peak) ; centY = centY_temp(selc_peak);
num_peak = numel(centX);

f3=figure(3);
app.f3 =f3;

fig_num_f3 = fig_num*3 ;
subplot(2,3,1+fig_num_f3)
imshow(filt_img, [10,200]);hold on
plot(centX_temp,centY_temp,'r+')
title(['Before correction in ',ch_name])


subplot(2,3,2+fig_num_f3)
histogram(ppc) ; hold on 
xline(peak_background_fold,'r-') ; xlim([0,4])
title(['Peak-Background fold histogram in ',ch_name])

subplot(2,3,3+fig_num_f3)
imshow(filt_img, [10,200]);hold on
plot(centX,centY,'r+')
title(['After applied with threshold in ',ch_name])

% 2230917 updated end

%%
fig_num_f4 = fig_num*3;
f4=figure(4) ; subplot(2,3,3+fig_num_f4);
app.f4=f4;
imagesc(filt_dat); hold on
plot(centX,centY,'r+')
title(['Background substracted in ',ch_name])


subplot(2,3,1+fig_num_f4)
imshow(filt_img, [10,200]);hold on
plot(centX,centY,'r+')
title([ch_name])


subplot(2,3,2+fig_num_f4)
mesh(filt_img) ; hold on;
centZ =[] ;
for i=1:length(centX)
    temp_Z = filt_img(centY(i),centX(i));
    centZ = [centZ,temp_Z];
end
scatter3(centX,centY,centZ,'o')
title(['3D Intensity plot in ',ch_name])

%
%5x5

[rows,cols] = size(filt_dat);
total_region= zeros(rows,cols);
total_final = zeros(rows,cols);
total_raw = zeros(rows,cols);
for i =1:length(centX)
    cent1 = centX(i) ; cent2 = centY(i) ;
    psf=filt_dat(cent2-5:cent2+5,cent1-5:cent1+5);
    psf_cor = psf(psf~=0) ;
    [psf_mu, psf_std]= normfit(psf_cor) ;

    %
    
    tolerance = psf_std.*tolerance_parameter ;
    [region, finalImage,raw_data] = regionGrowing(filt_dat, cent2,cent1, tolerance, cell_selc);
    total_region = total_region +region;
    total_final = total_final+finalImage ;
    total_raw = total_raw + raw_data;
end
%subplot(1,2,2)
%imagesc(finalImage);
%
dilute_part = cell_selc.*(~total_region);
dilute_average_int = mean(dilute_part(:));
%%

f5=figure(5);%0+fig_num)
app.f5 = f5;

fig_num_f5 = fig_num*2 ;
subplot(2,2,1+fig_num_f5)
imagesc(filt_dat)
title(['Before applied with tolerance parameter in ',ch_name])

subplot(2,2,2+fig_num_f5)
imagesc(total_final)
title(['After applied with tolerance parameter in ',ch_name,'red:count, white:average_intensity_thresed count'])

%
labeledImage = bwlabel(total_region, 4);
Labels= unique(labeledImage) ;

num_cond= max(Labels) ;

%%250903 updated
average_int_cond_thres_mask = zeros(1,num_cond); 

%

size_cond = [] ; sum_int_cond = []; average_int_cond = [];
for i = 1:num_cond
    [xv, yv] = find(labeledImage==i) ;
    size_cond(i) = length(xv) ;
    tmp_int = [] ;
    for j = 1:length(xv)
        tmp_int(j) = total_raw(xv(j),yv(j));
    end
    
    sum_int_cond(i) = sum(tmp_int(:));
    average_int_cond(i)=sum_int_cond(i)/size_cond(i);
    
    %%250903 updated
    if average_int_cond(i) > average_int_cond_thres %250903 updated
        average_int_cond_thres_mask(i) = 1; 
        subplot(2,2,2+fig_num_f5)
        text(yv(1)+5,xv(1)+5,num2str(i),'color', 'white','FontSize',12);
    end
    %

    subplot(2,2,2+fig_num_f5)
    text(yv(1),xv(1),num2str(i),'color', 'red','FontSize',10);
    
end

num_cond_thresed = sum(average_int_cond_thres_mask) ;

Kp = average_int_cond./dilute_average_int;
temper = 310.15 ; %37'C -> Kelvin ;
G_transfer = -8.31446261815324*temper*log(Kp); %J/mol


img_data_collection= struct('cell_selc',cell_selc,'num_cond',num_cond,'average_int_cond',average_int_cond,'dilute_average_int',dilute_average_int,...
'sum_int_cond', sum_int_cond, 'size_cond',size_cond,'Kp',Kp,'G_transfer',G_transfer,'bg',background_setting,'peak_bg_fc',peak_background_fold,'tolerance',tolerance_parameter,'average_int_cond_thres',average_int_cond_thres,'num_cond_thresed',num_cond_thresed);
cond_img = total_final ;



end
