#########################################
#      Modelo Fijo Interacción I        #
#  Ejemplo con todas las observaciones  #
# Apuntes Genética Estadística (97-111) #
#########################################

# Limpiar el ambiente
rm(list=ls())

#Antes que todo: instalar y/o cargar los paquetes 
library (MASS)
library(readxl)

# 1. Importar el set de datos
DatosRan <- read_excel("C:/Users/lucas/Desktop/GENEST/DatosRan.xlsx")
View(DatosRan)

# 2. Convertir variables categóricas en factores
DatosRan$RAT <- as.factor(DatosRan$RAT)
DatosRan$SIRE <- as.factor(DatosRan$SIRE)

# Para comprobar que sea un factor utilizamos función str
str(DatosRan$RAT)
str(DatosRan$SIRE)

#Construir matriz de incidencia
Tabla1 <- table(DatosRan$RAT,
                DatosRan$SIRE, 
                dnn = c("RAT", "SIRE"))

#Construir matriz Y
Y <- matrix(data = DatosRan$GAN, nrow = 18)
View(Y)

#Generar la matriz de incidencia de efectos fijos (X)
Xb <- matrix(0, length(Y), 2)
Xv <- matrix(0, length(Y), 3)
ib <- cbind(1:length(Y), DatosRan$RAT)
iv <- cbind(1:length(Y), DatosRan$SIRE) #Sin valor alfa
Xb[ib] <- 1
Xv[iv] <- 1
X <- cbind(c(rep(1,length(Y))), Xb, Xv)

#Matriz de incidencia de efectos aleatorios (Z)
Z <- matrix(c(rep(1,2),rep(0,length(Y)),rep(1,5),rep(0,length(Y)),rep(1,1),rep(0,length(Y)),
              rep(1,2),rep(0,length(Y)),rep(1,3),rep(0,length(Y)),rep(1,5)),length(Y),6)

#Matriz identidad
XX <- t(X) %*% X #Se construye X'X
XZ <- t(X) %*% Z #Se construye X'Z
ZX <- t(Z) %*% X #Se construye Z'X
ZZ <- t(Z) %*% Z #Se construye Z'Z
XY <- t(X) %*% Y #Se construye X'Y
ZY <- t(Z) %*% Y #Se construye Z'Y
CFX <- cbind(XX,XZ)
CFZ <- cbind(ZX,ZZ)
LS <- rbind(CFX,CFZ) #Se forma la matriz del lado izquierdo (X'X)
la.qr <- qr(LS) #Rango de matriz LS
RS <- rbind(XY,ZY) #Se forma la matriz del lado derecho (X'y)

# Obtener inversa generalizada de LS
#1. Se calcula con función ginv(T) (Requiere del paquete MASS)
GT <- ginv(LS)
View (GT)

#2. Remover filas y/o columnas para deshacer la dependencia y obtener inversa de submatriz
T1 <- LS[c(-1,-2,-3,-4,-5,-6),c(-1,-2,-3,-4,-5,-6)]
T2 <- solve(T1)
#Añadir nuevas filas y columnas con rbind y cbind (ENCONTRAR FORMA MÁS EFICIENTE)
T3 <- rbind(c(rep(0,6)), c(rep(0,6)),c(rep(0,6)), c(rep(0,6)),c(rep(0,6)), c(rep(0,6)), T2)
T4 <- cbind(c(rep(0,12)), c(rep(0,12)),c(rep(0,12)), c(rep(0,12)),c(rep(0,12)), c(rep(0,12)), T3)
View(T4)

# 8. Obtener (X'y)
View(RS)

# Encontrar una de la soluciones de b mediante b = ginv(X'X)*(X'y)
b <- T4 %*% RS
View(b)

# Diferencias entre raciones y entre padres no son estimables ya que estas
#no están libres de los efectos de interacciones.

#Lo que es posible estimar en este tipo de modelo, por ejemplo, son las diferencias
#entre raciones dentro de un mismo padre. Por ejemplo, diferencias entre raciones 1 y 2
#dentro del PADRE 1
K1 <- matrix(c(0,1,-1,0,0,0,1,0,0,-1,0,0),12,1)
t(K1) %*% b

# Medias Mínimo Cuadráticas
#Para este ejemplo, las medias mínimo cuadráticas (MMC) están dadas por:
# Crear una matriz K
K <- matrix(c(1,rep(1/2,2),rep(1/3,3),rep(1/6,6),
               1,1,0,rep(1/3,6),0,0,0,
               1,0,1,rep(1/3,3),0,0,0,rep(1/3,3),
               1,rep(1/2,2),1,0,0,1/2,0,0,1/2,0,0,
               1,rep(1/2,2),0,1,0,0,1/2,0,0,1/2,0,
               1,rep(1/2,2),0,0,1,0,0,1/2,0,0,1/2,
               1,1,0,1,0,0,1,0,0,0,0,0,
               1,1,0,0,1,0,0,1,0,0,0,0,
               1,1,0,0,0,1,0,0,1,0,0,0,
               1,0,1,1,0,0,0,0,0,1,0,0,
               1,0,1,0,1,0,0,0,0,0,1,0,
               1,0,1,0,0,1,0,0,0,0,0,1),12,12)

#Obtener MMC mediante la fórmula K'b
MMC <- t(K) %*% b
MMC

# Suma de Cuadrados
#SCT
SCT <- t(Y) %*% Y
SCT
#SCR
SCR <- t(b) %*% RS
SCR
#SCE
SCE <- SCT-SCR
SCE
#VARE
VARE <- SCE/(length(Y)-la.qr$rank)
VARE <- as.vector(VARE)
VARE
#SCM
SCM <- length(Y)*(mean(Y)^2)
SCM
#SCRm
SCRm <- SCR-SCM
SCRm

# 12. Calcular valor F para hipótesis del modelo
F <- (SCRm/5)/VARE
F

# Suma de Cuadrados Pertinente (Q)
#Ejemplo: La diferencia entre ración 1 y ración 2 dentro del padre 1
# Hipótesis: H0: r1 - r2 + rp11 - rp21 = 0

#Diferencias entre raciones 1 y 2 dentro del PADRE 1
K1 <- matrix(c(0,1,-1,0,0,0,1,0,0,-1,0,0),12,1)
Dif1 <- t(K1) %*% b
#Suma de Cuadrados apropiada
Q1 <- t(Dif1) %*% solve((t(K1) %*% T4 %*% K1)) %*% Dif1
Q1
#Valor F para esta hipótesis
F1 <- (Q1/1)/VARE
F1

#Ejemplo: La diferencia entre ración 1 y ración 2 dentro del padre 2
# Hipótesis: H0: r1 - r2 + rp12 - rp22 = 0

#Diferencias entre raciones 1 y 2 dentro del PADRE 2
K2 <- matrix(c(0,1,-1,0,0,0,0,1,0,0,-1,0),12,1)
Dif2 <- t(K2) %*% b
Dif2
#Suma de Cuadrados apropiada
Q2 <- t(Dif2) %*% solve((t(K2) %*% T4 %*% K2)) %*% Dif2
Q2
#Valor F para esta hipótesis
F2 <- (Q2/1)/VARE
F2

# Pruebas de hipótesis de diferencias entre lasmedias mínimo cuadráticas (MMC)
#de los efectos principales.

#La siguiente hipótesis prueba la igualdad de las MMCs de las raciones:
K3 <- matrix(c(0,1,-1,0,0,0,rep(1/3,3),rep(-1/3,3)),12,1)
Dif3 <- t(K3) %*% b
Dif3
#Suma de Cuadrados apropiada
Q3 <- t(Dif3) %*% solve((t(K3) %*% T4 %*% K3)) %*% Dif3
Q3
#Valor F para esta hipótesis
F3 <- (Q3/1)/VARE
F3

#La siguiente hipótesis prueba la igualdad de las MMCs de los padres:
K4 <- matrix(c(0,0,0,1,0,-1,1/2,0,-1/2,1/2,0,-1/2,
               0,0,0,0,1,-1,0,1/2,-1/2,0,1/2,-1/2),12,2)
Dif4 <- t(K4) %*% b
Dif4
#Suma de Cuadrados apropiada
Q4 <- t(Dif4) %*% solve((t(K4) %*% T4 %*% K4)) %*% Dif4
Q4
#Valor F para esta hipótesis
F4 <- (Q4/2)/VARE
F4

#La siguiente hipótesis prueba la igualdad de las MMCs de las interacciones:
K5 <- matrix(c(0,0,0,0,0,0,1,0,-1,-1,0,1,
               0,0,0,0,0,0,0,1,-1,0,-1,1),12,2)
Dif5 <- t(K5) %*% b
Dif5
#Suma de Cuadrados apropiada
Q5 <- t(Dif5) %*% solve((t(K5) %*% T4 %*% K5)) %*% Dif5
Q5
#Valor F para esta hipótesis
F5 <- (Q5/2)/VARE
F5

#La siguiente hipótesis prueba la igualdad de las MMCs de los padres 1 y 2:
#H0: MMCp1 = MMCp2
K6 <- matrix(c(0,0,0,1,-1,0,1/2,-1/2,0,1/2,-1/2,0),12,1)
Dif6 <- t(K6) %*% b
Dif6
#Suma de Cuadrados apropiada
Q6 <- t(Dif6) %*% solve((t(K6) %*% T4 %*% K6)) %*% Dif6
Q6
#Valor F para esta hipótesis
F6 <- (Q6/1)/VARE
F6
