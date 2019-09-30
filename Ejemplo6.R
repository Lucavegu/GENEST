#############################################
# Modelos de Clasificación de Doble Entrada #
#    Modelo de Clasificación Cruzada        #
#   Apuntes Genética Estadística (75-93)    #
#############################################

#Ejecutar en R
#!/usr/local/bin/Rscript

# Limpiar el ambiente
rm(list=ls())

#Antes que todo: cargar los paquetes 
library (MASS)
library(readxl)

# 1. Importar el set de datos
Ejemplo_Tarea6 <- read_excel("C:/Users/lucas/Desktop/GENEST/Ejemplo_Tarea6.xlsx")
View(Ejemplo_Tarea6)

# 2. Convertir variables categóricas en factores
Ejemplo_Tarea6$RAT <- as.factor(Ejemplo_Tarea6$RAT)
Ejemplo_Tarea6$SIRE <- as.factor(Ejemplo_Tarea6$SIRE)

#Para comprobar que sea un factor utilizamos función str
str(Ejemplo_Tarea6$RAT)
str(Ejemplo_Tarea6$SIRE)

# 3. Obtener tabla con frecuencias absolutas de aparición según variables
Tabla1 <- table(Ejemplo_Tarea6$RAT,
                Ejemplo_Tarea6$SIRE, 
                dnn = c("RAT", "SIRE"))
Tabla1

# 4. Construir matriz Y
Y <- matrix(data = Ejemplo_Tarea6$GDP, nrow = 25)
View(Y)

# 5. Generar la matriz de diseño X
Xr <- matrix(0, length(Y), 3) #3 indica el número de niveles (columnas) del efecto RAT
Xs <- matrix(0, length(Y), 4) #4 indica el número de niveles (columnas) del efecto SIRE
ir <- cbind(1:length(Y), Ejemplo_Tarea6$RAT)
is <- cbind(1:length(Y), Ejemplo_Tarea6$SIRE)
Xr[ir] <- 1
Xs[is] <- 1
X <- cbind(c(rep(1,length(Y))), Xr, Xs)
View(X)
#Además, es posible construir la matriz de incidencia (N), mediante:
N <- crossprod(Xr, Xs)

#Obtener tabla con sumatoria de observaciones
#PENDIENTE

# 6. Obtener la matriz resultante (T) de la multiplicación de X con la matriz traspuesta de X (t(X))

T <- t(X) %*% X
View(T)

#6.1 Obtener el rango de la matriz T
#Se calcula la descomposición QR de la matriz
la.qr <- qr(T)
#Se extrae el atributo rank de la.qr
print(c("El rango de la matriz es", la.qr$rank), quote = F)

# 7. Obtener la inversa generalizada de la matriz T
#7.1 Se calcula con función ginv(T) (Requiere del paquete MASS)
GT <- ginv(T)
View (GT)

#7.2 Remover filas y/o columnas para deshacer la dependencia y obtener inversa de submatriz
T1 <- T[c(-1,-2),c(-1,-2)]
T2 <- solve(T1)
#Añadir nuevas filas y columnas con rbind y cbind
T3 <- rbind(c(rep(0,6)), c(rep(0,6)), T2)
T4 <- cbind(c(rep(0,8)), c(rep(0,8)), T3)
View(T4)

# 8. Obtener (X'y)
U <- t(X) %*% Y
View(U)

# 9. Encontrar una de la soluciones de b mediante b = ginv(X'X)*(X'y)
#Puede ser utilizando GT o T4 como inversa generalizada de X'X
b <- GT %*% U
View(b)

b1 <- T4 %*% U
View(b1)

# 10. Obtener Y estimado (Xb)
Yest <- X %*% b
View(Yest)

Yest1 <- X %*% b1
View(Yest1)

# 11. Obtener sumas de cuadrados
#SCT
SCT <- t(Y) %*% Y
SCT
#SCR
SCR <- t(b)%*%t(X)%*%Y
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
F <- (SCRm/17)/VARE
F

## Funciones estimables

#Las diferencias entre raciones son estimables. Por ejemplo, la diferencia entre la ración 1 y 2 es:
k <- matrix(c(0,1,-1,0,0,0,0,0))
t(k) %*% b

#La varianza de la diferencia entre ración 1 y ración 2 está dado por:
k <- matrix(c(0,1,-1,0,0,0,0,0))
k1 <- k[c(2,3),]
Vard12 <- t(k1) %*% GT[c(2,3),c(2,3)] %*% k1 * VARE
Vard12
#O bien
t(k) %*% GT %*% k * VARE

#El error estándar de esta diferencia es:
sqrt(Vard12)

#Si deseo estimar en forma simultánea diferencias entre raciones utilizo una matriz adecuada
K <- matrix(c(0,1,0,-1,0,0,0,0,
              0,0,1,-1,0,0,0,0),8,2)
Vard1323 <- t(K) %*% GT %*% K * VARE
Vard1323

#Los errores estándar de estas diferencias son (utilizando los elementos de la diagonal):
sqrt(diag(Vard1323))

## Medias Mínimo Cuadráticas (MMC)

# Crear una matriz K
KM <- matrix(c(1,1,0,0,rep(1/4,4),
               1,0,1,0,rep(1/4,4),
               1,0,0,1,rep(1/4,4),
               1,rep(1/3,3),1,0,0,0,
               1,rep(1/3,3),0,1,0,0,
               1,rep(1/3,3),0,0,1,0,
               1,rep(1/3,3),0,0,0,1),8,7)
View(KM)

#Obtener MMC mediante la fórmula K'b
MMC <- t(KM) %*% b
View(MMC)

#La varianza de estas MMC se calculan de la siguiente forma
VarMMC <- t(KM) %*% GT %*% KM * VARE
View(VarMMC)

#El error estándar de las MMC es la raíz cuadrada de los elementos de la diagonal
sqrt(diag(VarMMC))

## Análisis de Varianza

#Si nuestro interés es detectar diferencias entre raciones es posible construir hipótesis a ser probadas.
#Por ejemplo ->  H0: R1-R2=0 / R2-R3=0

#Recordando las diferencias estimadas entre raciones
D1323 <- t(K) %*% b
D1323

#Es necesario calcular una suma apropiada de cuadrados (tipo III)
Qr <- t(D1323) %*% solve((t(K) %*% GT %*% K)) %*% D1323
Qr
#El valor F para esta hipótesis es
Fr <- (Qr/2)/VARE
Fr

#Si queremos probar diferencias entre los padres la hipótesis a probar sería:
#H0: T1-T4=0 / T2-T4=0 / T3-T4=0
#Se construye una matriz K apropiada
K2 <- matrix(c(0,0,0,0,1,0,0,-1,
               0,0,0,-0,0,1,0,-1,
               0,0,0,0,0,0,1,-1),8,3)
View(K2)

#Se estiman las diferencias entre padres
DP <- t(K2) %*% b
DP

#Es necesario calcular una suma apropiada de cuadrados (tipo III)
Qp <- t(DP) %*% solve((t(K2) %*% GT %*% K2)) %*% DP
Qp
#El valor F para esta hipótesis es
Fp <- (Qp/3/VARE)
Fp
