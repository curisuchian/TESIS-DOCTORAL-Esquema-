%% Código para detección de sombras en imágenes aéreas de áreas selváticas
% El siguiente código es una adaptación de otro código producido por 
% B. Sirmacek y C. Unsalan, el cual aplica el método del índice invariante
% de color para realizar la detección de sombras
% 
% Referencias:
% B. Sirmacek and C. Unsalan, "Damaged Building Detection in Aerial Images 
% using Shadow Information", 4th International Conference on Recent Advances 
% in Space Technologies RAST 2009, Istanbul, Turkey, June 2009.
%
% C. Unsalan and K. L. Boyer, "Linearized vegetation indices based on a formal 
% statistical framework," IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 42, pp. 1575-1585, 2004. 
clear all
close all
clc
%% Read Images: se elige una imagen del archivo
[filename, pathname] = uigetfile({'*.*'},'Browse');
name=[pathname,filename];
im = imread(name);
figure, subplot(1,2,1), imshow(im);
figure, subplot(1,2,1), imshow(im);
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
figure, imshow(im)
f = figure();
% Máscara de sombras
shadow_mask(:,:) = shadow_ratio>umbral;
shadow_mask(1:5,:) = 0;
shadow_mask(end-5:end,:) = 0;
shadow_mask(:,1:5) = 0;
shadow_mask(:,end-5:end) = 0;
% Marcado de sombras
shadow_mask = bwareaopen(shadow_mask, 100);
[x,y] = find(imdilate(shadow_mask(:,:),strel('disk',2,0))-shadow_mask(:,:));
imshow(im), hold on, plot(y,x,'.b'), title(['Umbral= ',num2str(umbral)]);





%%


