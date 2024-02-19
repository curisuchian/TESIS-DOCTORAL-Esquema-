%% Shadow Detection Source Code
% Shared by Beril Sirmacek
% For Academic & Educational Usage Only
% Please consider citing following reference articles.
% 
% B. Sirmacek and C. Unsalan, "Damaged Building Detection in Aerial Images 
% using Shadow Information", 4th International Conference on Recent Advances 
% in Space Technologies RAST 2009, Istanbul, Turkey, June 2009.
%
% C. Unsalan and K. L. Boyer, "Linearized vegetation indices based on a formal 
% statistical framework," IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 42, pp. 1575-1585, 2004. 

%%
clear all
close all
clc

%% Read Images: se elige una imagen del archivo
[filename, pathname] = uigetfile({'*.*'},'Browse');
name=[pathname,filename];
%imageSegmenter
im = imread(name);
figure, subplot(1,2,1), imshow(im);

figure, subplot(1,2,1), imshow(im);

% NOTE: You might need different median filter size for your test image.
% Filtrado de ruido
r = medfilt2(double(im(:,:,1)), [3,3]); 
g = medfilt2(double(im(:,:,2)), [3,3]);
b = medfilt2(double(im(:,:,3)), [3,3]);


%% Calculate Shadow Ratio: cociente entre la resta del canal azul menos el rojo, sobre la suma de ambos
shadow_ratio = ((4/pi).*atan(((b-r))./(b+r)));
umbral = quantile(shadow_ratio(:),.85);
histograma = hist(shadow_ratio(:),100);
hist_cum = cumsum(histograma);
X = -1:2/99:1; %vector para eje de abscisas (valores posibles del índice invariante de color)
figure, plot(X,hist_cum,X,histograma),yyaxis right, title('shadow ratio');

%shadow_ratio = ((4/pi).*atan(((r-g))./(g+r)));
%shadow_ratio = ((4/pi).*atan(((b-g))./(b+g)));
%figure, imshow(shadow_ratio, []); colormap(jet); colorbar;
%imshow(im); hold on,
figure, imshow(im)
f = figure();
%f.WindowState = 'maximized';
%figure2=figure('Position', [1,1 , 1000, 1000]);
%shadow_mask = zeros(size(shadow_ratio));
% NOTE: You might need a different threshold value for your test image.
% You can also consider using automatic threshold estimation methods.

shadow_mask(:,:) = shadow_ratio>umbral;
%figure, imshow(shadow_mask, []); 

shadow_mask(1:5,:) = 0;
shadow_mask(end-5:end,:) = 0;
shadow_mask(:,1:5) = 0;
shadow_mask(:,end-5:end) = 0;

% NOTE: Depending on the shadow size that you want to consider,
% you can change the area size threshold
shadow_mask = bwareaopen(shadow_mask, 100);
[x,y] = find(imdilate(shadow_mask(:,:),strel('disk',2,0))-shadow_mask(:,:));
% [xm,ym] = find(imdilate(BW180a(:,:),strel('disk',2,0))-BW180a(:,:));
imshow(im), hold on, plot(y,x,'.b'), title(['Umbral= ',num2str(umbral)]);
% figure, imshow(im), hold on, plot(ym,xm,'.r'), title(['Manual']);




%%


