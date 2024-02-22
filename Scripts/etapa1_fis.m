clc
close all
clear all

%Leo una imagen...
%A = imread('ST1.jpg');
%% Read Images: se elige una imagen del archivo
[filename, pathname] = uigetfile({'*.*'},'Browse');
name=[pathname,filename];
%imageSegmenter
A = imread(name);
figure, imshow(A)
% A = A(size(A,1)/2:size(A,1)-1,size(A,2)/2:size(A,2),:);
% A = A(:,1:size(A,2)/2,:);
n=size(A,1);
% if rem(n,2)==1
%     n=n-1;
% end

m=size(A,2);
% if rem(m,2)==1
%     m=m-1;
% end
% imhist(A(:,:,1));
% figure
% imhist(A(:,:,2));
% figure
% imhist(A(:,:,3));
if size(A,3) ==3        %comprueba si es a color
    %conversión a intensidad
    A=double(A); %debe hacerse esta declaración para evitar que sature el valor int8 cuando suma
    V = (A(:,:,1)+A(:,:,2)+A(:,:,3))/3;
    V=(V/255);
end
%Para probar con bandas separadas... RGB
V = A(:,:,1); %Rojo
% V = A(:,:,2); %Verde
% V = A(:,:,3); %Azul
%
% V = histeq(V);
% Sub-Sistema caracteristico de entrada
% n=size(V,1);
% m=size(V,2);

u=[-1/2+1/n:1/n:1/2];
v=[-1/2+1/m:1/m:1/2];

V(V==0)=0.1;
Vsoma=log(V);
Vuv=fft2(Vsoma);

% Sub-Sistama linear
gh=1;
gl=0;
c=9;
Do=0.1;
K=1;

uu=u.^2;vv=v.^2;
uuu=repmat(uu',1,m);
vvv=repmat(vv,n,1);
DUV=sqrt(uuu+vvv);
Huv=K*(1-((gh-gl)*(1-exp(-c*(DUV.^2)/Do^2))+gl));


if rem(m,2)==1
    Huv=[Huv(n/2:n,m/2:m)     Huv(n/2:n,1:m/2);...
     Huv(1:n/2-1,m/2:m)   Huv(1:n/2-1,1:m/2)];
else Huv=[Huv(n/2:n,m/2:m)     Huv(n/2:n,1:m/2-1);...
     Huv(1:n/2-1,m/2:m)   Huv(1:n/2-1,1:m/2-1)];
end

Suv=Huv.*Vuv;
 
Vo=ifft2(Suv);

%remocao da fase
Vo=abs(Vo);
Vout=exp(Vo);

% imshow(A)
% figure


Vt = Vout/max(1*max(Vout));
figure, imshow(Vt)



m = size(Vt,1);
n = size(Vt,2);
%Deformo matriz para que se adapte como entrada al sistema FIS...
B = reshape(Vt(:,:,1),m*n,1);

%Evalúa según la intensidad de los píxeles
C = evalfis(B, discriminadorint);

C=reshape(C,m,n);
D = C.*V;
[W,K] = bwboundaries(C,'noholes');
regiones=regionprops(K);
intensidades=regionprops(K,Vt,'MeanIntensity');
imshow(C)
%Almaceno cada área en una posicón de un vector
for k = 1:length(W)
    areas(k)=regiones(k).Area;
    int_promedio(k)=intensidades(k).MeanIntensity;
end
%normalizo...bah, multiplico por 100

areas_n=areas*100/prod(size(V));
%Filtro`áreas grandes y chicas...e incorporo selección por intensidad
%promedio de cada área
ind_area = evalfis([areas_n' int_promedio'], segundaetapa);
%copia de binarizada
selectas = C;
k=1;
Wprima=W;
intensidades_actualizadas=intensidades;
regiones_actualizadas=regiones;
while k < length(regiones_actualizadas)
    %si es mayor que cero...
    if ind_area(k)>0
        k = k+1;  
    else
%         selectas()=
        regiones_actualizadas(k)=[];
        ind_area(k)=[];
        Wprima(k)=[];
        intensidades_actualizadas(k)=[];
    end

end

% area_fil=areas.*ind_area';
%visualizar en la imagen las sombras marcadas

figure, imshow(V)
hold on
%Superpongo las sombras detectadas en primera etapa (azul)
for k = 1:length(W)
    boundary = W{k};
    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2)
end

%Superpongo las sombras filtradas en la segunda etapa (se observan en
%violeta)
selectas=selectas*0;
for k = 1:length(Wprima)
    boundary = Wprima{k};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
%     for q = 1:length(boundary)
        xmin = min(boundary(:,1));
        ymin = min(boundary(:,2));
        xmax = max(boundary(:,1));
        ymax = max(boundary(:,2));
        selectas(xmin:xmax, ymin:ymax)=1;
%     end
end
% figure
% imshow((selectas));
% hold on
% for k = 1:length(Wprima)
%     boundary = Wprima{k};
%     plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
% end