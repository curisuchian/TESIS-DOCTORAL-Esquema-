                                                                                                        install.packages("ggplot2")
install.packages("imager")
install.packages("mixtools")
install.packages("baseline")
install.packages("Matrix")

install.packages("rgdal")
install.packages("plotly")
#install.packages("EcoGenetics")

install.packages("BiocManager")
BiocManager::install("EBImage")

library(mixtools)
library(baseline)
library(ggplot2)
library(imager)
library(Matrix)

library("EBImage")

#Load and plot the RGB file
#data.file <- '/mnt/usb-WD_Elements_25A2_575852314531383859503438-0:0-part1/CHB/PASANTIA/imagenes_forestal/data/DJI_0805.jpg'
data.file <- '../imagenes_forestal/data/original_referencia.jpg'
im <- load.image(data.file)

#GPerforms the conversion of the image from RGB to HSL colorspace
im_hsl <- RGBtoHSL(im)

#Extracts the L component from HSL colorspace image
#im_L <- im_hsl[200:500,200:500,3]
im_L <- im_hsl[,,3]

param <- normalmixEM(as.vector(im_L))

entrada3 <- im_L
#values of image that are lower than mean are set to 0
salida3 <- (entrada3>param$mu[1])*entrada3
display(salida3)
hist(entrada3)
hist(salida3)

gridSearch <- function(fun, im, a) { #Función definida para el procesamiento del paso 4 en un rango determinado de uno de los parámetros;
    #un argumento es la función, que es específica del paso 4, otro es la imagen a ser procesada y el tercer argumento es el parámetro, correspondiente al ancho de ventana del filtro rolling ball
  
    len_a <- length(a)
    len_im <- dim(im)
    
    #result_array <- array(len_im,len_a,len_b)
    result_array <- array(0, dim=c(len_im[1],len_im[2],len_a))
    k <- 1 # ka se usa como índice de la imagen resultado
    for (i in c(1:len_a)){
            result_array[,,k] <- fun(im,a[i])
            k <- k + 1
    }
    
    return (result_array)
    }

paso4 <- function(im,ws){
    entrada4 <- im
    #inversion of grayscale image and addition of maximum grayscale value
    im_L_inv <- entrada4*(-1)+max(entrada4)
    window_size <- ws
    #both baselines bline1 and bline2 are computed considering one input as inverted image and the other input as transposed inverted image
    bline1 <- baseline(t(im_L_inv),wm=window_size, ws=window_size, method='rollingBall')
    bline2 <- baseline(im_L_inv,wm=window_size, ws=window_size, method='rollingBall')
    #smooth image
    im_smooth <- pmax(t(bline1@baseline)*(-1),(bline2@baseline)*(-1))
    salida4 <- im_smooth
    return (salida4)
}

imagenes <- gridSearch(paso4,im_L,c(6,9,12,15,18,24))

display(imagenes+max(-imagenes), all = TRUE)


gridSearch2 <- function(fun, im,im2,a) { #Función definida para el procesamiento del paso 5 en un rango determinado de uno de los parámetros;
    #un argumento es la función, que es específica del paso 5, otro es la imagen a ser procesada, el tercer argumento es la imagen de salida del paso 4, y el cuarto argumento correspondiente al radio en píxeles
  
    len_a <- length(a)
    len_im <- dim(im)
    
    #result_array <- array(len_im,len_a,len_b)
    result_array <- array(0, dim=c(len_im[1],len_im[2],len_a))
    k <- 1 # k se usa como índice de la imagen resultado
    for (i in c(1:len_a)){
            result_array[,,k] <- fun(im,im2,a[i])
            k <- k + 1
    }
    
    return (result_array)
    }

paso5 <- function(im,im2,radio){ #dos de los argumentos a ser pasados son imágenes, y uno es el parámetro radio
entrada5 <- im
#Top hat
#Structuring element consists in a circular shape of determined radius
#radio <- 14 #radius of 7 pixels, corresponding to crown diameter; con 14 se da un mejor resultado usando la imagen original de referencia
mask <- px.circle(radio)

abertura <- mopening(as.cimg(entrada5),mask,real_mode = FALSE)
t_hat <- as.cimg(entrada5) - abertura
abertura <- abertura>0

maskara <- abertura[,,1,1]
salida5 <- (maskara*im2)*(-2)



salida5 <- ((!maskara&entrada5)*entrada5+salida5)
    return (salida5)

    }

imagenes2 <- gridSearch2(paso5,salida3,imagenes[,,3],c(3,6,14,20,30,60)) 

#display(-imagenes2+max(imagenes2), all = TRUE)

display(imagenes2, all = TRUE)

paso6 <- function(im,percentil) {
    entrada6 <- im
#a normal distribution (n_gaps) is generated, using the parameters that were found with normalmixEM (eg. the media and standard deviation)
n_gaps <- rnorm(length(entrada6), mean = param$mu[1], sd = param$sigma[1])
noventaynueve <- quantile(n_gaps,percentil)
salida6 <- (!(entrada6[,]<noventaynueve))*entrada6
return(salida6)
}


imagenes6 <- gridSearch(paso6,imagenes2[,,3],c(0.99,0.001,0.01,0.1,0.5,0.9))
display(imagenes6,all=TRUE)

paso7 <- function (entrada7,percentil) {
    #entrada7 <- salida6
    #ti <- proc.time()
    mat_riz<-cbind(0,0,0,entrada7,0,0,0) #se rellenan tres columnas con ceros por izquierda y por derecha
    mat_riz<-rbind(0,0,0,mat_riz,0,0,0) #se rellenan tres filas con ceros por arriba y por abajo
    MNZ <- entrada7*0 #MNZ es una matriz de la misma dimensión que mat_riz completa con ceros

    for (i in 4:dim(entrada7)[1]+2) { #i es el índice que recorre las columnas
       for (j in 4:dim(entrada7)[2]+2) { #j es el índice que recorre las filas
           a <- i-3
           b <- i+3
           c <- j-3
           d <- j+3
           ventana <- mat_riz[a:b,c:d]
           MNZ[i-2,j-2] <- nnzero(ventana)
       }

     }
    setentaycinco <- quantile(MNZ,percentil)
   
    huecos_copas <- (MNZ>setentaycinco)*entrada7
    salida7 <- huecos_copas
   
    mascara7 <- MNZ>setentaycinco
   return (salida7)
    }


imagenes7 <- gridSearch(paso7,imagenes6[,,1],c(0.1,0.3,0.5,0.75,0.8,0.9))
#0.1,0.3,0.5,0.75,0.8,0.9
display(imagenes7,all=TRUE)

paso8 <- function(entrada8,distancia) {




dist_min <- distancia #distancia mínima 7 píxeles
mat_riz<-cbind(0,0,0,entrada8,0,0,0) #se rellenan tres columnas con ceros por izquierda y por derecha
mat_riz<-rbind(0,0,0,mat_riz,0,0,0)


    MNZ <- entrada8*0 #MNZ es una matriz de la misma dimensión que mat_riz completa con ceros

    for (i in 4:dim(entrada8)[1]+2) { #i es el índice que recorre las columnas
       for (j in 4:dim(entrada8)[2]+2) { #j es el índice que recorre las filas
           a <- i-3
           b <- i+3
           c <- j-3
           d <- j+3
           ventana <- mat_riz[a:b,c:d]
           MNZ[i-2,j-2] <- nnzero(ventana)
       }
        }

im_dist <- distmap(MNZ)
for (i in 4:dim(entrada8)[1]+2) { #i es el índice que recorre las columnas
   for (j in 4:dim(entrada8)[2]+2) { #j es el índice que recorre las filas
       a <- i-3
       b <- i+3
       c <- j-3
       d <- j+3
       if (im_dist[i-2,j-2] > dist_min) {
            ventana <- mat_riz[a:b,c:d]
       media <- 0
       for (n in 1:4) {
           media <- max(ventana)/4 + media
           ventana[which.max(ventana)] <- 0
       }
       entrada8[i-2,j-2] <- media
       }
           
      
   }
   
 }
#media
salida8 <- entrada8
return (salida8)
    }

imagenes8 <- gridSearch(paso8,imagenes7[,,4],c(2,4,7,10,15,30))
#0.1,0.3,0.5,0.75,0.8,0.9
display(imagenes8,all=TRUE)

paso9 <- function (entrada9,umbral) {

    #Top hat y bottom hat
    mask <- imfill(6,6,val=1)
    top_hat <- as.cimg(entrada9) - mopening(as.cimg(entrada9),mask)
    bottom_hat <-  mclosing(as.cimg(entrada9),mask) - as.cimg(entrada9)
    im_filt <- as.cimg(entrada9) + bottom_hat - top_hat

    percentil <- quantile(im_filt,umbral)
    salida9 <- (im_filt>percentil)

    #hist(salida9)
}

salida9 <-gridSearch(paso9,imagenes8[,,3],c(0.001,0.01,0.1,0.2,0.5,0.9))
display(salida9,all = TRUE)

display(salida9[,,6])

 mask <- imfill(6,6,val=1)
    top_hat <- as.cimg(imagenes8[,,3]) - mopening(as.cimg(imagenes8[,,3]),mask)
    bottom_hat <-  mclosing(as.cimg(imagenes8[,,3]),mask) - as.cimg(imagenes8[,,3])
    im_filt <- as.cimg(imagenes8[,,3]) + bottom_hat - top_hat
    #0.001,0.01,0.1,0.2,0.5,0.9
hist(im_filt)
    (percentil <- quantile(im_filt,0.9))

entrada10 <- salida9
distancia <- distmap(entrada10)
display(distancia, all = TRUE)



