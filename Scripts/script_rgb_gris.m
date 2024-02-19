clear all;
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% carga de imagen
%A=imread('PNI2_original.jpg');
%% Read Images: se elige una imagen del archivo
[filename, pathname] = uigetfile({'*.*'},'Browse');
name=[pathname,filename];
%imageSegmenter
im = imread(name);
imshow(im)
A=histeq(im);
A=im;
figure, imshow(A)
n=size(A,1);
% if rem(n,2)==1
%     n=n-1;
% end

m=size(A,2);
% if rem(m,2)==1
%     m=m-1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(A,3) ==3        %comprueba si es a color
    %conversión a intensidad
    A=double(A); %debe hacerse esta declaración para evitar que sature el valor int8 cuando suma
    V = (A(:,:,1)+A(:,:,2)+A(:,:,3))/3;
    V=(V/255);
end

% Sub-Sistema caracteristico de entrada
% n=size(V,1);
% m=size(V,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=[-1/2+1/n:1/n:1/2];
v=[-1/2+1/m:1/m:1/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if rem(n,2)==1
     disp('n es impar')
else
     disp('n es par')
end
 if rem(m,2)==1
    Huv=[Huv(n/2:n,m/2:m)     Huv(n/2:n,1:m/2);...
     Huv(1:n/2-1,m/2:m)   Huv(1:n/2-1,1:m/2)];
 disp('... y m es impar')
else Huv=[Huv(n/2:n,m/2:m)     Huv(n/2:n,1:m/2-1);...
     Huv(1:n/2-1,m/2:m)   Huv(1:n/2-1,1:m/2-1)];
 disp('... y m es par')
end
%%%%%%%%%%%%%%%%- aplicación de filtro homomórfico -%%%%%%%%%%%%%%%%%%
Suv=Huv.*Vuv;
%%%%%%%%%%%%%%%%- transformada inversa de Fourier -%%%%%%%%%%%%%%%%%%%%
Vo=ifft2(Suv);

%%%%%%%%%%%%%%%- remoción de fase -%%%%%%%%%%%%%%%%%%%%%%

Vo=abs(Vo);
Vout=exp(Vo);


figure, imshow(V)
figure, imshow(Vout/max(Vout(:)))
figure
Vt=Vout/max(1*max(Vout));
%%%%%%%%%%%%%%%%%%%%%- BINARIZACIÓN -%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%aplico un threshold
mayor=0.75;
menor=0.75;

Vt(Vt>mayor)=1;
Vt(Vt<menor)=0.1;
imshow(Vt)
IM_BIN = Vt;                            %Copia d eimagen binarizada de paso anterior
[x,y] = find(imdilate(IM_BIN(:,:),strel('disk',2,0))-IM_BIN(:,:));
figure, imshow(im), hold on, plot(y,x,'.b'), title(['asi es...']);
% IM_BIN = Vt(200:350,550:720);

n = 25;                                  %tamaño de ventana de inspección
% n = size(matriz_inspeccion,2);
FILAS = round(size(IM_BIN,1)/n)-1; %1012 filas/ 3 = 337
COLUMNAS = round(size(IM_BIN,2)/n)-1; %1600 columnas/ 3 =533
UMBRAL = 0.45;                          %45%
CONTADOR_AB = 0;
CONTADOR_AV = 0;
MATRIZ_SECUNDARIA = zeros(FILAS,COLUMNAS);
SOMBRAS = zeros(FILAS,COLUMNAS);
k = 1;          %contador de sombras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               INICIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=0:FILAS-1
    for j=0:COLUMNAS-1
        intensidad_a = sum(sum(IM_BIN(1+i*n:n*(1+i),1+j*n:n*(1+j))))/n^2;
        
        if  intensidad_a >= UMBRAL            %Si la cantidad de píxeles de sombra superan el umbral...
            %asignar el valor intensidad a un elemento correspondiente de
            %una matriz
            %Sentencia para agrupar sombras en distintas filas pero que son
            %contiguas...
%             if i > 1 && (sum(sum(IM_BIN(1+i*n-n:n*i,1+j*n:n*(1+j))))/n^2) > UMBRAL
%                 SOMBRAS(1+i*n-n,1+j*n) = SOMBRAS(1+i*n-n,1+j*n) + intensidad_a;
%                 intensidad_b = sum(sum(IM_BIN(1+i*n:n*(1+i),1+j*n+n:n*(2+j))))/n^2;
%           
           MATRIZ_SECUNDARIA(i+1,j+1) = 1;
        else %asignar cero
            MATRIZ_SECUNDARIA(i+1,j+1) = 0;
           
        end         %fin if UMBRAL
        
    end             %fin ciclo columnas
end                 %fin ciclo filas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Enmascaramiento de sombras
for i=0:FILAS-1
    for j=0:COLUMNAS-1
        IMAGEN_ENMASCARADA(i*n+1:n*(i+1),1+j*n:n*(1+j)) = MATRIZ_SECUNDARIA(i+1,j+1)*IM_BIN(1+i*n:n*(1+i),1+j*n:n*(1+j));
            
    end             %fin ciclo columnas
end
%figure, imshow(A)
figure,imshow(IMAGEN_ENMASCARADA)
disp('Cantidad de sombras detectadas:')
size(bwboundaries(IMAGEN_ENMASCARADA),1)

[x,y] = find(imdilate(IM_BIN(:,:),strel('disk',2,0))-IM_BIN(:,:));
figure, imshow(im), hold on, plot(y,x,'.b'), title(['asi es...']);
[x,y] = find(imdilate(IMAGEN_ENMASCARADA(:,:),strel('disk',2,0))-IMAGEN_ENMASCARADA(:,:));
hold on, plot(y,x,'.r'), title(['asi es...']);