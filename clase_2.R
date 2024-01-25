setwd("C:/Users/Gustavo/Documents/2023/Curso_DataScience/Clases")
getwd()

var_hola <- "Hola mundo"
(vaar_logical <- T)# El paréntesis ayuda a que se imprima la función
(vector <- c("Hola mundo", "Adios mundo"))# La "c" significa concatenación
class(vector)
typeof(vaar_logical)

a <- c(4, 6, 8, 10,12)
b <- c(3, 5, 7, 9, 5)

length(a)
length(b)

a[3] #Corchetes para extrael "elementos"
b[4]

d <- c(a,b)# la "c" sirve para concatenar los vectores, se crea un vector más grande
d

var_ordenadas <- sort(c(a,b), decreasing = T)

var_ordenadas

5:10 # se crea una sucesión de números. Esto es un vector 
10:1

vector_by <- seq(from = 0, to = 200, by = 10)
vector_by

vector_tres <- seq(from = 2, to = 15, by = 3)# FROM - TO - BY 
vector_tres

vector_by3 <- seq(1, 10, 3)# FROM - TO - BY: Es una forma simplificada
vector_by3

rep(5, times = 6) # Se repite el elemento n veces. n se sustituye
rep(a, 2)

c(1, 2, 3) + c(7, 8, 9, 10,2,6)# Tienen que ser múltiplos

print(a + b)#Es una operación de vectores

a <- c(4, 6, 8, 10,12)
b <- c(3, 5, 7, 9, 5)

a +b
a - b
a/b
a*b
a**2
a1 <- a*0.5 + b**2
a1[5]

sort(a1, decreasing = F)

######## MATRICES########
m <- matrix(1:9, nrow = 3, ncol = 3, byrow = T)
m
m[3,1]# El primer número indica el número de filas y el segundo el número de columna

m[1,] 

suma_vector <- c(1,2) + m #suma el 1 a la primera columna, el 2 a la segunda, la tercera a la primera
suma_vector

n <- matrix(1:24, 4, 6, byrow = T)
n
dim(n)

n[n>10] #Extraer un conjunto
which(n > 10)#Extrae la posición de una matriz

(a <- 2:6)
(b <- 5:9)

length(a)
length(b)

cbind(a,b)#Combina las columnas
rbind(a,b)#Combina las filas

n <- matrix(1:16, 4,4)
n

p <- matrix(10:150, 6, 4)
p

apply(p, 1, mean)# 1 se refeire a las filas. Promedio por filas. Cada una
apply(p,2,sort)# 2 se refiere a columnas. Ordenar por columnas de la matriz 

t(p)# Transposición
diag(p) #Diagonal de la matriz 

############## LISTAS Y DATA FRAMES###########

milista <- list(nombre = "Pepe", no.hijos = 3, edades.hijos = c(4, 7, 9))
milista

otralista <- list(nombre = "Gus", num_mascotas = 3, altura = c(1.5,.6,1.4))
str(otralista)

otralista$nombre

x <- 6:8
y <- c("A", "B", "C")

tabla <- data.frame(edad = x, grupo = y)
tabla

str(tabla)

tabla[1]# Extrae la columna
tabla[,1]# Extrae columnna igual
tabla$edad#Extrae columna igual

for (i in 1:5) print(i*i)

w <- rnorm(1:10)
w

wsq <- 0

for (i in 1:10) {
  wsq[i] <- w[i]**2
  print(wsq[i])
}

wsq

for(i in 1:10){
 print(5*i + 3)
}

count <- 0
while (count > 10) {
  print(count)
  count <- count + 1
}

x <- runif(1,0,5)

if(z > 2){
  y <- TRUE
}else {
  y <- FALSE
}

z;y
x
?runif

z <- runif(1,1,10)
z

wsq <- 0
w <- rnorm(44)
for (i in 1:15) {
  wsq[i] <- w[i]**3 + 12
  print(wsq[i])
}

#### Sesion 2 ######
mean(tabla$edad)

paste("La media de la edad es de:", mean(tabla$edad))
paste("El total de las edades es de:", sum(tabla$edad))#paste se usa para crear una leyenda

summary(tabla)#es como el describe de python
dim(tabla)#Tres filas, dos columnas

tabla$Sexo <- c("M", "M", "H")#crear una nueva columna
tabla


tabla$Sexo <- NULL# Eliminar columna
tabla

############# Descarga y lectura de data sets#########

netflix <- read.csv("netflix_titles.csv")
dim(netflix)

colnames(netflix)

mayor_2015 <- netflix["release_year" > 2015]
estrenos_2015 <- netflix[netflix$release_year > 2015,]
estrenos_2015

write.csv(estrenos_2015, "netflix_2015.csv")


######## Instalar Paquetes########

install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)

## Reto 2

Amazon <- read.csv("best_sellers_with_categories.csv")

colnames(Amazon)

t(amazon)
tAmazon <- as.data.frame(t(Amazon))
tAmazon

colnames(tAmazon)
#names(tamazon) 

libro_max <- which.max(tAmazon["Price", ])
libro_min <- which.min(tAmazon["Price", ])

print(libro_max)
print(libro_min)

tAmazon$V1

?names

libro2_max <- which.max(Amazon$Price)
Amazon$Price[70]
Amazon$Name[70]
libro2_min <- which.min(Amazon$Price)
Amazon$Price[43]

names(tAmazon) <- tAmazon[1,]
names(tAmazon)


######## Loops y pseudocódigo #########

w <- rnorm(20)# Guardar 20 números aleatorios
w
i <- 1
wsq <- w[i]**2
wsq[i]

i <- 2
w[i]

for(i in 1:10){
  wsq[i] <- w[i]**2
  print(wsq[i])
}

w <- rnorm(20)
wsq <- 0

for (i in 1:10) {
wsq[i] <- w[i]**2
print(wsq[i])
}

for(i in 1:10) print(i*i)#Primero se multiplcia 1*1, luego 2*2, 3*3,...10*10

#### Reto 3 ######

wsq
for (i in 1:10) {
  print(5*i + 3)
}

### While ##
count <- 0
while(count < 10){
  print(count)
  count <- count + 1
}

### If ####

x <- runif(1,0,10)
if(x > 5){
  y <- TRUE
}else {
  y <- FALSE
}

 y

w <- rnorm(44)
wsq <- 0

for(i in 1:15){
  wsq[i] <- w[i]**3 + 12
  print(wsq[i])
}

m <- w[1:15]
n <- wsq

d <- data.frame(m, n)
d








############ CLASE 2#############
#################################

### MEDIDAS DE TENDENCIA CENTRAL ####


datos <- read.csv("Top_500_Songs.csv")

colnames(datos)
sort(datos$released)

datos_dos <- read.csv("governors_county.csv")
colnames(datos_dos)

order(datos_dos$current_votes)
sum(datos_dos$current_votes)
length(datos_dos$current_votes)

summary(datos_dos)
min(datos_dos$current_votes)
max(datos_dos$current_votes)

mean(datos_dos$current_votes)
median(datos_dos$current_votes)

library(modeest)

mfv(datos_dos$current_votes)
votos <- quantile(datos_dos$current_votes)
votos[4]

var(datos_dos$current_votes)
library(fmsb)

IQR(datos_dos$current_votes)
sd(datos_dos$current_votes, na.rm = TRUE)

#### BOX PLOT ###

x11(width = 10,height = 5)
graf <- boxplot(datos_dos$current_votes, main = "Votos actuales")

x <- c(29,13,62,4,63,96,1,90,50,46)
quantile(x, 0.25)
sort(x)
quantile(x, c(0.25,0.50,0.75))
quantile(x, seq(0.1,0.9, by = 0.1)) #Calcular los deciles, divide el conjunto de datos en 10 partes iguales

IQR(x)# Se usa cuando hay datos extremos
var(x)
sd(x)

### Reto 1 ####

set.seed(134)
x <- round(rnorm(1000,175,6),1)#Población de 1000, con media de 185 y sd de 6
x
?set.seed

mean(x)
median(x)
Mode(x)

quantile(x, seq(0.1, 0.9, by = 0.1))

IQR(x)
var(x)
sd(x)

str(iris)
colnames(iris)
summary(iris)
summary(1:100)

summary(mtcars)

set.seed(57)
x <- rnorm(35)
e <- rnorm(35)
y <- 5 + 2*x + e

modelo <- lm(y~x)
modelo
summary(modelo)

head(mtcars)
tail(mtcars)

View(iris) # Muestra una hoja de cálculo que muestran los datos

operacion <- function(a, b){
  a*b + 10
}# Esto es similar a crear un def en python

operacion(5, 4)
operacion(13, 3)

x <- c(7,7,7,1,3,9)
table(x)# Esto imprime los valores únicos y su frecuencia

f.abs <- table(x)

moda <- function(x){
f.abs <- table(x)
maxfreq_abs <- max(f.abs)
pos_max <- which(f.abs == maxfreq_abs)# Da la posición donde se encuentra la máxima frecuencia absoluta
print("La(s) moda(s) es(son): ")
print(names(f.abs[pos_max]))
paste("Con una frecuencia de:", unique(as.vector(f.abs[pos_max])))
}

?str_extract_all()

x <- sample(1:100, 100, replace = T)
moda(x)
table(x)

str(mtcars)
summary(mtcars)
head(mtcars)
View(mtcars)

x <- mtcars
x

calcular_mediana <- function(vector){
  ord_datos <- sort(vector)# ordenar los valores
  longitud <- length(ord_datos)
  if(longitud %% 2 == 0){ # ver si el conjunto es par o impar
    mediana <- (ord_datos[longitud/2] + ord_datos[(longitud/2) + 1]) / 2
  } else {
    mediana <- ord_datos[(longitud + 1) / 2]
  }
  print(mediana)
}

y <- c(25,6,84,92,101,45,99,43)
sort(y)

calcular_mediana(y)

#### Ejemplo 3 ###

library(dplyr)

url1 <- "https://data.humdata.org/hxlproxy/data/download/time_series_covid19_confirmed_global_narrow.csv?dest=data_edit&filter01=explode&explode-header-att01=date&explode-value-att01=value&filter02=rename&rename-oldtag02=%23affected%2Bdate&rename-newtag02=%23date&rename-header02=Date&filter03=rename&rename-oldtag03=%23affected%2Bvalue&rename-newtag03=%23affected%2Binfected%2Bvalue%2Bnum&rename-header03=Value&filter04=clean&clean-date-tags04=%23date&filter05=sort&sort-tags05=%23date&sort-reverse05=on&filter06=sort&sort-tags06=%23country%2Bname%2C%23adm1%2Bname&tagger-match-all=on&tagger-default-tag=%23affected%2Blabel&tagger-01-header=province%2Fstate&tagger-01-tag=%23adm1%2Bname&tagger-02-header=country%2Fregion&tagger-02-tag=%23country%2Bname&tagger-03-header=lat&tagger-03-tag=%23geo%2Blat&tagger-04-header=long&tagger-04-tag=%23geo%2Blon&header-row=1&url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv"
url2 <- "https://data.humdata.org/hxlproxy/data/download/time_series_covid19_deaths_global_narrow.csv?dest=data_edit&filter01=explode&explode-header-att01=date&explode-value-att01=value&filter02=rename&rename-oldtag02=%23affected%2Bdate&rename-newtag02=%23date&rename-header02=Date&filter03=rename&rename-oldtag03=%23affected%2Bvalue&rename-newtag03=%23affected%2Binfected%2Bvalue%2Bnum&rename-header03=Value&filter04=clean&clean-date-tags04=%23date&filter05=sort&sort-tags05=%23date&sort-reverse05=on&filter06=sort&sort-tags06=%23country%2Bname%2C%23adm1%2Bname&tagger-match-all=on&tagger-default-tag=%23affected%2Blabel&tagger-01-header=province%2Fstate&tagger-01-tag=%23adm1%2Bname&tagger-02-header=country%2Fregion&tagger-02-tag=%23country%2Bname&tagger-03-header=lat&tagger-03-tag=%23geo%2Blat&tagger-04-header=long&tagger-04-tag=%23geo%2Blon&header-row=1&url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_deaths_global.csv"

download.file(url = url1, destfile = "st19ncov-confirmados.csv", mode = "wb")
download.file(url = url2, destfile = "st19ncov-muertes.csv", mode = "wb")

conf <- read.csv("st19ncov-confirmados.csv")
muert <-read.csv("st19ncov-muertes.csv") 

Sconf <- conf[-1,]
Smu <- muert[-1,]# quitar filas

Sconf <- rename(Sconf, Pais = Country.Region, Fecha = Date, Infectados = Value)

str(Sconf)

Sconf <- mutate(Sconf, Fecha = as.Date(Fecha, "%Y-%m-%d"), Infectados = as.numeric(Infectados))
#cambiar el tipo de datos, por ejemplo de "character" a "fecha"


# Se hace lo mismo para los datos de muertes
Smu <- select(Smu, Country.Region, Date, Value) # Seleccionamos país, fecha y acumulado de muertos
Smu <- rename(Smu, Pais = Country.Region, Fecha = Date, Muertos = Value) # Renombramos
Smu <- mutate(Smu, Fecha = as.Date(Fecha, "%Y-%m-%d"), Muertos = as.numeric(Muertos)) # Transformamos

str(Smu)

Scm <- merge(Sconf, Smu) # Unimos infectados y muertos acumulados para cada fecha


head(Scm)
tail(Scm)

mex <- filter(Scm, Pais == "Mexico") # Seleccionamos sólo a México
mex <- filter(mex, Infectados != 0) # Primer día de infectados

head(mex)

mex <- mutate(mex, Infectados = as.numeric(Infectados), Muertos = as.numeric(Muertos))

str(mex)

a <- c(4,8,12,20)
diff(a)# diff resta un valor menos el anterior

mex <- mutate(mex, NI = c(1, diff(Infectados))) # Nuevos infectados por día
mex <- mutate(mex, NM = c(0, diff(Muertos))) # Nuevos muertos por día
mex <- mutate(mex, Letalidad = round(Muertos/Infectados*100, 1)) # Tasa de letalidad

mex <- mutate(mex, IDA = lag(Infectados), MDA = lag(Muertos)) # Valores día anterior, lag desplaza elementos
mex <- mutate(mex, FCI = Infectados/IDA, FCM = Muertos/MDA) # Factores de Crecimiento
mex <- mutate(mex, Dia = 1:dim(mex)[1]) # Días de contingencia

lag(letters)
lag(a)# el último elemento desaparece y el primero es NA

write.csv(mex, "C19Mexico.csv", row.names = FALSE)


#### Ejemplo 4 ####
cbind(1:10, 11:20, 21:30)# Deben de tener la MISMA LONGITUD (numero de filas)
cbind(1:10, matrix(11:30, ncol =2))# hace lo mismo
cbind(data.frame(x = 1:10, y = 11:20), z = 21:30)# combina por columnas

df1 <- data.frame(x = 1:5, y = 6:10, z = 16:20)
df1

df2 <- data.frame(x = 51:55, y = 101:105, z = 151:155)
df2
df1; df2
rbind(df1, df2)# Deben tener el mismo # de columnas y los mismos nombres (rowbind)


### Ejemplo 5 ####

x <- matrix(1:49, ncol = 7)
x
apply(x, 1, mean) # cálculo de la media para las filas, con 1
apply(x, 2, median) # cálculo de la mediana para las columnas, con 2

u1011 <- "https://www.football-data.co.uk/mmz4281/1011/SP1.csv"
u1112 <- "https://www.football-data.co.uk/mmz4281/1112/SP1.csv"
u1213 <- "https://www.football-data.co.uk/mmz4281/1213/SP1.csv"
u1314 <- "https://www.football-data.co.uk/mmz4281/1314/SP1.csv"

download.file(url = u1011, destfile = "SP1-1011.csv", mode = "wb")
download.file(url = u1112, destfile = "SP1-1112.csv", mode = "wb")
download.file(url = u1213, destfile = "SP1-1213.csv", mode = "wb")
download.file(url = u1314, destfile = "SP1-1314.csv", mode = "wb")


#lista <- lapply(dir(), read.csv) # Guardamos los archivos en lista, se imprtan todos los archivos csv en r, sin cargar uno por uno


lista <- lapply(c(5,5,5), rnorm)# lapply da como resultado una lista
lista

head(lista[[1]]); head(lista[[2]]); head(lista[[3]])

data <- do.call(rbind, lista)# do.call da como resultado un data frame 
head(data)
dim(data)


#### Reto 3 ####

url1 <- "https://www.football-data.co.uk/mmz4281/1718/D1.csv"
url2 <- "https://www.football-data.co.uk/mmz4281/1819/D1.csv"
url3 <- "https://www.football-data.co.uk/mmz4281/1920/D1.csv"
url4 <- "https://www.football-data.co.uk/mmz4281/2021/D1.csv"

download.file(url = url1, destfile = "D1.csv", mode = "wb")
download.file(url = url2, destfile = "D2.csv", mode = "wb")
download.file(url = url3, destfile = "D3.csv", mode = "wb")
download.file(url = url4, destfile = "D4.csv", mode = "wb")

x <- read.csv("D1.csv")
y <- read.csv("D2.csv")
z <- read.csv("D3.csv")
a <- read.csv("D4.csv")

datos_select <- select(x, Date, HomeTeam, AwayTeam, FTHG,FTAG,FTR)
datos_select2 <- select(y, Date, HomeTeam, AwayTeam, FTHG,FTAG,FTR)
datos_select3 <- select(z, Date, HomeTeam, AwayTeam, FTHG,FTAG,FTR)
datos_select4 <- select(a, Date, HomeTeam, AwayTeam, FTHG,FTAG,FTR)


data_final <- do.call(rbind, list(datos_select, datos_select2, datos_select3, datos_select4))

data_final
dim(data_final)


### Ejemplo 6 ####
head(airquality)
library(dplyr)

str(airquality)
bien <- complete.cases(airquality)# Muestra los NA en un data frame
bien
sum(bien)# muestra las filas completas
airquality[bien,]#Muestra todas las filas que NO TIENEN NA's

data <- select(airquality, Ozone:Temp)
apply(data, 2, mean)
apply(data, 2, mean, na.rm = T)# na.rm elimina los NA. Reemplaza con "0" los NA

# `na.omit` devuelve el objeto con casos incompletos eliminados

(m1 <- apply(na.omit(data), 2, mean))# Aquí se eliminan los NA en lugar de llenarlos con "0"
b <- complete.cases(data)
b

(m2 <- apply(data[b,], 2, mean))

identical(m1, m2)


########### Sesión 3/ Análisis exploratorio de datos######

library(ggplot2)
names(mtcars)# ver nombres de columnas
head(mtcars)

ggplot(mtcars, aes(x=wt, y = disp, colour = mpg )) + 
  geom_point()  # Tipo de geometría, intenta utilizar alguna otra, "grafico de dispersion"

names(mtcars)
ggplot(mtcars, aes(x=wt, y = disp, colour = am )) + 
  geom_point() +   
  theme_gray() +  # Temas (inteta cambiarlo)
  facet_wrap("am")  # Lo divide por el núm de cilindros


names(mtcars)
ggplot(mtcars, aes(x=wt, y = disp, colour = am )) + 
  geom_point() +   
  theme_gray() +   # Temas (inteta cambiarlo)
  facet_wrap("am") +  # Lo divide por el núm de cilindros
  xlab('Peso 1000 lbs')+  # Nombre en los ejes
  ylab('Desplazamiento')


box <- read.csv("boxp.csv")
head(box)
names(box)
hist(box$Mediciones, breaks = (seq(0,300, 20)), 
     main = "Histograma de Mediciones",
     xlab = "Mediciones",
     ylab = "Frecuencia")

box <- na.omit(box) #Evitar el Warning de filas con NA´s

library(dplyr) 

box %>% #Este símbolo significa "pipe"
  ggplot() + 
  aes(Mediciones) +
  geom_histogram(binwidth = 20)

box %>%
  ggplot() + 
  aes(Mediciones) +
  geom_histogram(binwidth = 20, col="black", fill = "blue") + 
  ggtitle("Histograma de Mediciones") +
  ylab("Frecuencia") +
  xlab("Mediciones") + 
  theme_void()

library(ggplot2)
library(dplyr)
alumnos <- read.csv("BD_Altura_Alunos.csv", sep = ";")
head(alumnos)

hist(alumnos$Altura, breaks = (seq(0,300, 20)), 
     main = "Histograma de Altura",
     xlab = "Altura",
     ylab = "")


alumnos %>% 
  ggplot() + 
  aes(Altura) +
  geom_histogram(binwidth = 10)#Este es el histograma con ggplot

hist(alumnos$Altura)# Solo el histograma

alumnos %>%
  ggplot() + 
  aes(Altura) +
  geom_histogram(binwidth = 15, col="purple", fill = "yellow") + 
  ggtitle("Histograma de Altura") +
  ylab("Frecuencia") +
  xlab("Alturas") + 
  theme_void()

names(alumnos)


##### Ejemplo 3 #########

(my_scatplot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point())
(my_scatplot <- ggplot(mtcars, aes(x = wt, y = mpg)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = T))  # modelo lineal, cambia el parametro `se`, este hace referencia al intervalo de confianza
# se = "Standard error" crea una banda de confianza del 95% para predicciones

my_scatplot + xlab('Weight (x 1000lbs)') + ylab('Miles per Gallon')
(my_scatplot <- ggplot(mtcars, aes(x = wt, y = mpg, col = cyl)) + geom_point())

my_scatplot + labs(x = 'Weight (x1000lbs)',
                   y = 'Miles per Gallon', 
                   colour = 'Number of\n Cylinders')


my_scatplot + facet_wrap("cyl")# esto separa los gráficos
my_scatplot + facet_grid(am~cyl)


#### Reto 2 ###

jugadores <- read.csv("players_stats.csv")
head(jugadores)

hist(jugadores$MIN)

?geom_vline


jugadores %>%
  ggplot() + 
  aes(MIN) +
  geom_histogram(binwidth = 50, col="purple", fill = "yellow") + 
  geom_vline(aes(xintercept = mean(MIN)), color = "red", linetype = "dashed", size = 1) +
  ggtitle("Histograma total de minutos") +
  ylab("Frecuencia") +
  xlab("Minutos totales") + 
  theme_void()

jugadores <- na.omit(jugadores)

jugadores %>%
  ggplot() + 
  aes(Age) +
  geom_histogram(binwidth = 1, col="black", fill = "blue") + 
  geom_vline(aes(xintercept = mean(Age)), color = "red", linetype = "dashed", size = 1) +
  ggtitle("Histograma de Edad") +
  ylab("Frecuencia") +
  xlab("Edad") + 
  theme_void()

jugadores$Age
length(jugadores$Age)

scater <- ggplot(jugadores, aes(x = Weight, y = Height)) + geom_point(color = "blue") +
labs(x = "Peso", y = "Altura", title = "Peso vs Altura")
scater


##### Ejemplo 4 ######

box <- read.csv("boxp.csv")
head(box)
summary(box)

box <- na.omit(box)

box <- mutate(box, Categoria = factor(Categoria), Grupo = factor(Grupo))

summary(box)
names(box)
head(box)
bien <- complete.cases(data)
data <- data[bien,]
data <- mutate(data, Categoria = factor(Categoria), Grupo = factor(Grupo))
str(box)
View(box)

ggplot(box, aes(x = Categoria, y = Mediciones, fill = Grupo)) + geom_boxplot() +
  ggtitle("Boxplots") +
  xlab("Categorias") +
  ylab("Mediciones")


ggplot(box, aes(x = Categoria, y = Mediciones, fill = Grupo)) + geom_boxplot() +
  scale_fill_discrete(name = "Dos Gps", labels = c("G1", "G2")) + 
  ggtitle("Boxplots") +
  xlab("Categorias") +
  ylab("Mediciones")

#scale_fill_discrete = cambia los nombres

#### Ejemplo 5 #####

library(scales) 

Sys.time()



#### Sesion 4 #####

# Ejemplo 1 #

# VARIABLE CONTINUA- variable donde no hay huecos (estatura)

#DISTRIBUCION NORMAL Y T STUDENT

library(ggplot2) # Utilizaremos estos paquetes para algunas gráficas
library(reshape2)

x <- seq(-4,4,0.01)*6 + 175
y <- dnorm(x, mean = 175, sd = 6)

?dnorm

plot(x, y, type = "l", xlab= "", ylab = "")
abline(v = 175, lwd = 2, lty = 2)# línea que dibuja la media 

pnorm(q = 180, mean = 175, sd = 6)# pnorm- probabilidad acumulada, cuál es la probabilidad de...
par(mfrow = c(2,2))# hace que haya mas de una gráfica
plot(x, y, type = "l", xlab = "", ylab = "")

polygon(c(min(x), x[x<=180], 180), c(0, y[x<=180], 0), col = "red")#el área roja es la probabilidad = 0.7976

pnorm(q = 165, mean = 175, sd = 6)
plot(x, y, type = "l", xlab = "", ylab = "")
polygon(c(min(x), x[x<=165], 165), c(0, y[x<=165], 0), col = "yellow")

pnorm(q = 180, mean = 175, sd = 6) - pnorm(q = 165, mean = 175, sd = 6)
plot(x, y, type = "l", xlab = "", ylab = "")
polygon(c(165, x[x<=165 & x<=180], 180), c(0, y[x<=165], 0), col = "green")


dev.off()# quita el panel para volver a poner una sola gráfica


### cuantil ##

(b <- qnorm(p = 0.75, mean = 175, sd = 6))

set.seed(7563)# se generan números aleatorios
muestra <- rnorm(n = 1000, mean = 175, sd = 6)
length(muestra); mdf <- as.data.frame(muestra)
tail(mdf)

ggplot(mdf, aes(muestra))+
  geom_histogram(colour = "red",
                 fill = "blue",
                 alpha = 0.3,
                 birthwith = 3)+
  geom_density(aes(y = 3*...count...))+
  geom_vline(xintercept = mean(mdf$muestra), linetype = "dashed", color = "black")+
  ggtitle("Histograma")+
  theme_dark()

## t de student ##
#Media es 0, se utiliza para crear intervalos de confianza e hipótesis
# se usa con menos de n = 30
# cuantil es el proceso inverso a la probabilidad

### Reto 1 ###


#pnorm(q = 18, mean = 16, sd = 1, lower.tail = FALSE)
m <- pnorm(q = 18, mean = 16, sd = 1)
m

n <- 1-m
n

### Ejemplo 2 ###

x <- seq(0, 5, 0.02)
plot(x, dexp(x, rate = 2), type = "l", lwd = 2, ylab = "")
title(main = "Función de Densidad Exponencial", ylab = "f(x)",
      sub = expression("Parámetro " ~ lambda == 2))
text(x = 3, y = 1.5, labels = expression(f(x)==2*exp(-2*x) ~ " para x "  >= 0))
text(x = 3, y = 1.3, labels = paste("0 en otro caso"))
text(x = 1, y = 1, labels = expression("E(X) = " ~ 1/lambda == 1/2), col = 2)
text(x = 3, y = 0.5, labels = expression("DE(X) = " ~ 1/lambda == 1/2), col = 4)


set.seed(10) # Para reproducir posteriormente la muestra
(m1.4 <- rexp(n = 4, rate = 2))
mean(m1.4)

set.seed(64) # Para reproducir las muestras en el futuro
(m5.3 <- sapply(X = rep(3, 5), FUN = rexp, 2))# sapply genera una matriz
#rexp, 2 es la tasa, que se aplica a rep


set.seed(465) # Para reproducir las muestras en el futuro
m1000.7 <- sapply(X = rep(7, 1000), FUN = rexp, 2)
media1000.7 <- apply(m1000.7, 2, mean)
mdf <- as.data.frame(media1000.7)
tail(mdf)

ggplot(mdf, aes(media1000.7)) + 
  geom_histogram(colour = 'green', 
                 fill = 'orange',
                 alpha = 0.7) + # Intensidad del color fill
  geom_vline(xintercept = mean(media1000.7), linetype="dashed", color = "black") + 
  ggtitle('Histograma para las 1000 medias') + 
  labs(x = 'medias', y = 'Frecuencia')+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) 

mean(media1000.7); 1/2 # Media de las 1000 medias y media de la población de la cual vienen las 1000 muestras
sd(media1000.7); (1/2)/sqrt(7) # DE de las 1000 medias y DE de la población de la cual vienen las 1000 muestras dividida por la raíz del tamaño de la muestra



pnorm(q = 58, mean = 60, sd = 8/10)

## Ejemplo 4 #### Contraste de hipótesis

set.seed(174376)
muestra <- rexp(n = 56, rate = 4.1); 1/4.1 # media = 1/6 aprox 0.1667 (media de la población)
tail(as.data.frame(muestra))

m2 <- rexp(n = 63, rate = 3.4);1/3.4

#hay dos medias poblacionales
1/4.1-1/1.3# diferencia de medias reales

#HO: muestra-m2 = 0 vs HA: muestra-2 diferente de 0
#(contraste de dos colas) lo determina la HA

z0 <- (mean(muestra)-mean(m2)-0)/sqrt(var(muestra)/56 + var(m2)/63) 
z0#diferencia de promedios dividido entre raiz cuadrada de la varianza. Se está estandarizado 
#Estandarización

z0# Clacular un estadístico de prueba 

#región de rechazo (de dos colas) de un nivel
#nivel de significancia alpha = 0.05 (error tipo 1)


z.025 <- qnorm(p = 0.025, lower.tail = FALSE)
z.025

(pvalue <- 2*pnorm(z0, lower.tail = FALSE))
x <- seq(-4,4,0.01)
y <- dnorm(x)



polygon(c(z0, x[x>=z0], max(x)), c(0, y[x>=z0], 0), col="red")
axis(side = 1, at = z0, font = 2, padj = 1, lwd = 2)

#"When the p-value is low, H0 must go"

set.seed(1776)
s1 <- rnorm(n = 23, mean = 175, sd = 3)
tail(as.data.frame(m1))
s2 <- rnorm(n = 20, mean = 160, sd = 3)
tail(as.data.frame(m2))
175-160


t.test(x = s1, y = s2,
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = FALSE)



### Reto 3 ###

muestra <- c(V1 = 13, V2 = 17, V3 = 20, V4 = 17, V5 = 20, V6 = 20, V7 = 18, V8 = 18, V9 = 16, V10 = 19, 
             V11 = 13, V12 = 17, V13 = 15, V14 = 19, V15 = 16, V16 = 19, V17 = 22, V18 = 10, V19 = 13, V20 = 21)

m1 <- unname(muestra)
m1

muestra

#HO: H = 15
#HA: H > 15

?t.test

t.test(x = m1,
       alternative = "greater",
       mu = 15)


###### Sesión 5 ########
#Modelos Predictivos

nyc <- read.csv("nyc.csv", header = TRUE)

head(nyc)
colnames(nyc)

tail(nyc, 2) 
dim(nyc)
attach(nyc)

pairs(~ Price + Food + Decor + Service, data = nyc, gap = 0.4, cex.labels = 1.5)
#Graficos de dispersión 

m1 <- lm(Price ~ Food + Decor + Service + East)#lm es "Modelo Lineal"/ modelo predictivo
summary(m1)
m2 <- lm(Price ~ Food + Decor + East)
summary(m2)

m2 <- update(m1, ~.-Service)#actualizar el modelo, quitando la varaible
summary(m2)


mfull <- lm(Price ~ Food + Decor + Service + East + 
              Food:East + Decor:East + Service:East)
summary(mfull)

anova(m2,mfull)#Contraste de hipótesis - se acepta la H0

summary(m2)


plot(m2$fitted.values, Price, xlab = "Valores ajustados", ylab = "Price")
abline(lsfit(m2$fitted.values, Price))

pairs(~ Food + Decor, data = nyc, gap = 0.4, cex.labels = 1.5)

StanRes2 <- rstandard(m2)#Estimaciones de los errores del modelo
par(mfrow = c(2, 2))
plot(Food, StanRes2, ylab = "Residuales Estandarizados")
plot(Decor, StanRes2, ylab = "Residuales Estandarizados")
plot(East, StanRes2, ylab = "Residuales Estandarizados")

#Los errores son aleatorios, significa que el modelo ha capturado toda la información del dataset


qqnorm(StanRes2)
qqline(StanRes2)

dev.off()


#### Reto 1 ###

datos <- read.csv("advertising.csv", header = TRUE)
dim(datos)
head(datos)
colnames(datos)
attach(datos)

pairs(~ Sales + TV + Radio + Newspaper + Sales, data = datos, gap = 0.4, cex.labels = 1.5)

m1 <- lm(Sales ~ TV + Radio + Newspaper)
summary(m1)


### Ejemplo 2###
#Clasificación (varaible de respuesta cualitativa)

library(dplyr)
library(e1071)
library(ggplot2)
install.packages("ISLR")
library(ISLR)

?Default
head(Default)
tail(Default)
dim(Default)
str(Default)

ggplot(Default, aes(x = balance, y = income, colour = default)) + 
  geom_point() + facet_wrap('student') + 
  theme_grey() + ggtitle("Datos Default")

set.seed(2020)
train = sample(nrow(Default), 
               round(nrow(Default)/2))#5000 valores de manera aleatoria
tail(Default[train, ])#Filtra 5000 filas del df Default. Train-conjunto de entrnamiento

ggplot(Default[train, ], 
       aes(x = balance, y = income, colour = default)) + 
  geom_point() + facet_wrap('student') + 
  theme_dark() + ggtitle("Conjunto de entrenamiento")# <- 

ggplot(Default[-train, ], #Se filtran los otros 5000 datos restantes
       aes(x = balance, y = income, colour = default)) + 
  geom_point() + facet_wrap('student') + 
  theme_light() + ggtitle("Conjunto de prueba")# <- 

# Ahora utilizamos la función `tune` junto con la función `svm` para 
# seleccionar el mejor modelo de un conjunto de modelos, los modelos 
# considerados son aquellos obtenidos al variar los valores de los 
# parámetros `cost` y `gamma`. Kernel Radial

tune.rad = tune(svm, default~., data = Default[train,], 
                kernel = "radial", 
                ranges = list(
                  cost = c(0.1, 1, 10, 100, 1000), 
                  gamma = seq(0.01, 10, 0.5)
                ) 
)

#el . depués de la tilde ~ se usa para abreviar el poner las varaibles
# Se ha elegido el mejor modelo utilizando *validación cruzada de 10 
# iteraciones*

# summary(tune.rad)

# Aquí un resumen del modelo seleccionado

# summary(tune.rad$best.model)

best <- svm(default~.,  data = Default[train,],
            kernel = "radial",
            cost = 100,
            gamma = 1.51
)

mc <- table(true = Default[-train, "default"], 
            pred = predict(best, 
                           newdata = Default[-train,]))#Predict predice del data frame la varaible default 
mc

# El porcentaje total de aciertos obtenido por el modelo usando el 
# conjunto de prueba es el siguiente

round(sum(diag(mc))/sum(colSums(mc)), 5)

# Ahora observemos las siguientes proporciones

rs <- apply(mc, 1, sum)
r1 <- round(mc[1,]/rs[1], 5)
r2 <- round(mc[2,]/rs[2], 5)
rbind(No=r1, Yes=r2)

### Reto 2 ###

datos <- read.csv("datosclases.csv", header = TRUE)
library(dplyr)
colnames(datos)
head(datos)
summary(datos)
str(datos)

datos <- mutate(datos, y = as.factor(y))

library(ggplot2)
library(e1071)

dim(datos)
tail(datos)


ggplot(datos, aes(x = x.1, y = x.2, colour = y)) + 
  geom_point() + 
  theme_grey() + ggtitle("Datos Clases")


set.seed(120)
train = sample(nrow(datos), 
               round(nrow(datos)/2))
tail(datos[train, ])

ggplot(datos[train, ], 
       aes(x = x.1, y = x.2, colour = y)) + 
  geom_point() +  
  theme_dark() + ggtitle("Conjunto de entrenamiento")

ggplot(datos[-train, ], 
       aes(x = x.1, y = x.2, colour = y)) + 
  geom_point() + 
  theme_light() + ggtitle("Conjunto de prueba")

#### Sesion 6 ####

AP <- AirPassengers# Series de tiempo (ordenados cronológicamente)
AP

class(AP)# "ts"- time serie

start(AP); end(AP); frequency(AP)# En qué año empiza, termina y que frecuencia

summary(AP)

plot(AP, ylab = "Pasajeros (1000's)", xlab = "Tiempo", 
     main = "Reserva de pasajeros aéreos internacionales", 
     sub = "Estados Unidos en el periodo 1949-1960")

# Establecer el directorio de trabajo según corresponda
# setwd()
CBE <- read.csv("cbe.csv", header = TRUE)
CBE[1:4,]
class(CBE)

Elec.ts <- ts(CBE[, 3], start = 1958, freq = 12)#transformar a serie de tiempo
Beer.ts <- ts(CBE[, 2], start = 1958, freq = 12)
Choc.ts <- ts(CBE[, 1], start = 1958, freq = 12)

plot(cbind(Elec.ts, Beer.ts, Choc.ts), 
     main = "Producción de Chocolate, Cerveza y Electricidad", 
     xlab = "Tiempo",
     sub = "Enero de 1958 - Diciembre de 1990")

Global <- scan("global.txt")
Global.ts <- ts(Global, st = c(1856, 1), end = c(2005, 12), fr = 12)
Global.annual <- aggregate(Global.ts, FUN = mean)
plot(Global.ts, xlab = "Tiempo", ylab = "Temperatura en °C", main = "Serie de Temperatura Global",
     sub = "Serie mensual: Enero de 1856 a Diciembre de 2005")
plot(Global.annual, xlab = "Tiempo", ylab = "Temperatura en °C", main = "Serie de Temperatura Global",
     sub = "Serie anual de temperaturas medias: 1856 a 2005")


New.series <- window(Global.ts, start = c(1970, 1), end = c(2005, 12)) 
New.time <- time(New.series)
  

# Se debe elegir entre modelo aditivo o modelo multiplicativo cuando sea razonable suponer la descomposición
Elec.decom.A <- decompose(Elec.ts)#decompose-extraer nueva series de la serie original

plot(Elec.decom.A, xlab = "Tiempo", 
     sub = "Descomposición de los datos de producción de electricidad")

Tendencia <- Elec.decom.A$trend
Estacionalidad <- Elec.decom.A$seasonal
Aleatorio <- Elec.decom.A$random

Tendencia

plot(Elec.ts, 
     xlab = "Tiempo", main = "Datos de Producción de Electricidad", 
     ylab = "Producción de electricidad", lwd = 2,
     sub = "Tendencia con efectos estacionales aditivos sobrepuestos")
lines(Tendencia, lwd = 2, col = "blue")
lines(Tendencia + Estacionalidad, lwd = 2, col = "red", lty = 2)

ts.plot(cbind(Tendencia, Tendencia + Estacionalidad), 
        xlab = "Tiempo", main = "Datos de Producción de Electricidad", 
        ylab = "Producción de electricidad", lty = 1:2, 
        col = c("blue", "red"), lwd = 2,
        sub = "Tendencia con efectos estacionales aditivos sobrepuestos")

Tendencia[20] + Estacionalidad[20] + Aleatorio[20]
Elec.ts[20]


Elec.decom.M <- decompose(Elec.ts, type = "mult")

plot(Elec.decom.M, xlab = "Tiempo", 
     sub = "Descomposición de los datos de producción de electricidad")


Trend <- Elec.decom.M$trend
Seasonal <- Elec.decom.M$seasonal
Random <- Elec.decom.M$random

#Se multiplica la tendencia por la estacionalidad

plot(Elec.ts, 
     xlab = "Tiempo", main = "Datos de Producción de Electricidad", 
     ylab = "Producción de electricidad", lwd = 2,
     sub = "Tendencia con efectos estacionales multiplicativos sobrepuestos")
lines(Trend, lwd = 2, col = "blue")
lines(Trend * Seasonal, lwd = 2, col = "red", lty = 2)

ts.plot(cbind(Trend, Trend * Seasonal), 
        xlab = "Tiempo", main = "Datos de Producción de Electricidad", 
        ylab = "Producción de electricidad", lty = 1:2, 
        col = c("blue", "red"), lwd = 2,
        sub = "Tendencia con efectos estacionales multiplicativos sobrepuestos")

Trend[100]*Seasonal[100]*Random[100]
Elec.ts[100]


### Ejemplo 2 ####

set.seed(1)# Modelos estocásticos
w <- rnorm(100)#Ruido blanco, varaibles aletatorias de la misma distribución con media 0
#100 valores de ruido blanco
plot(w, type = "l", xlab = "")
title(main = "Ruido Blanco Gaussiano", xlab = "Tiempo")
#Varianza constante- hay mucha fluctuación

x <- seq(-3, 3, length = 1000)
hist(w, prob = T, ylab = "", xlab = "", main = "") 
points(x, dnorm(x), type = "l")
title(ylab = "Densidad", xlab = "Valores simulados de la distribución normal estandar",
      main = "Comparación de una muestra con su población subyacente")



acf(w, main = "")#acf es un gráfico- correlaciones lineales entre las varaibles aleatorias
title(main = "Función de Autocorrelación Muestral",#serie de tiempo
      sub = "Valores simulados de la distribución normal estandar")

#K = Lag- correlacion de todas las varaibles que se llevan unidades de tiempo
#las líneas que están dentro de las barras azules, significa que las variables no están correlacionedas 


set.seed(2)
x <- w <- rnorm(1000)
for(t in 2:1000) x[t] <- x[t-1] + w[t]

plot(x, type = "l", main = "Caminata Aleatoria Simulada", 
     xlab = "t", ylab = expression(x[t]), 
     sub = expression(x[t]==x[t-1]+w[t]))

acf(x, main = "")
title(main = "Correlograma para la caminata aleatoria simulada", 
      sub = expression(x[t]==x[t-1]+w[t]))
#Todas las variables están correlacionadas en el tiempo > 0, una depende de la anterior

acf(diff(x), main = "")
title(main = "Correlograma de la serie de diferencias", 
      sub = expression(nabla*x[t]==x[t]-x[t-1]))


## Modelos teóricos

#Modelo AR(p)-Autorregresivo de orden 1
rho <- function(k, alpha) alpha^k

plot(0:10, rho(0:10, 0.7), type = "h", ylab = "", xlab = "")
title(main = "Correlograma para un proceso AR(1)",
      ylab = expression(rho[k] == alpha^k),
      xlab = "lag k",
      sub = expression(x[t]==0.7*x[t-1]+w[t]))

plot(0:10, rho(0:10, -0.7), type = "h", ylab = "", xlab = "")
title(main = "Correlograma para un proceso AR(1)",
      ylab = expression(rho[k] == alpha^k),
      xlab = "lag k",
      sub = expression(x[t]==-0.7*x[t-1]+w[t]))
abline(h = 0)

## Modelo aurregresivo de orden 1
set.seed(1)
x <- w <- rnorm(100)
for(t in 2:100) x[t] <- 0.7 * x[t-1] + w[t]

plot(x, type = "l", xlab = "", ylab = "")
title(main = "Proceso AR(1) simulado",
      xlab = "Tiempo",
      ylab = expression(x[t]),
      sub = expression(x[t]==0.7*x[t-1]+w[t]))


###### CORRELOGRAMA #######
#Permite saber qué variables están correlacionadas
#


#Función de autocorrelación parcial
#Correlación de variables entre un determinado tiempo, quitando las correlaciones de en medio
acf(x, main = "")
title(main = "Correlograma del proceso AR(1) simulado", 
      sub = expression(x[t]==0.7*x[t-1]+w[t]))
pacf(x, main = "")
title(main = "Correlograma Parcial del proceso AR(1) simulado", 
      sub = expression(x[t]==0.7*x[t-1]+w[t]))


x.ar <- ar(x, method = "mle")
x.ar$order
x.ar$ar
x.ar$ar + c(-2, 2)*sqrt(x.ar$asy.var)

x


########## Proyecto ##########
#which.max(datos$Altitude)
#which.min(datos$Altitude)
#datos$Altitude[101]
#datos$Country.of.Origin[101]
############


#### Residuales = ruido blanco, es cero
### Modelos móviles de orden 3
## Un modelo puede dar origen a muchas series

acf(x, main = "")
title(main = "Correlograma para un proceso MA(3) simulado", 
      sub = expression(x[t] == w[t] + 0.8*w[t-1] + 0.6*w[t-2] + 0.4*w[t-3]))

#Wt = ruido blanco 

###### Reto 1 #######

set.seed(1)
x <- w <- rnorm(200)
for(t in 2:200) x[t] <- 0.5 * x[t-1] + w[t]

##Serie de tiempo##
plot(x, type = "l", xlab = "", ylab = "")
title(main = "Proceso AR(1) simulado",
      xlab = "Tiempo",
      ylab = expression(x[t]),
      sub = expression(x[t]==0.5*x[t-1]+w[t]))#Modelo Estacionario, serie distribuida alrededor de una serie horizontal
#constante en la serie y varianza

#Varianza constante = varaiblidad alrededor de un rectángulo 

##Correlograma
pacf(x, main = "")
title(main = "Correlograma del proceso AR(1) simulado", 
      sub = expression(x[t]==0.5*x[t-1]+w[t]))


##Ajuste del modelo autorregresivo

x.ar <- ar(x, method = "mle")
x.ar$order
x.ar$ar
x.ar$ar + c(-2, 2)*sqrt(x.ar$asy.var)

?ar


##### Modelos NO estacionarios ###

# Establecer el directorio de trabajo según corresponda
CBE <- read.csv("cbe.csv", header = TRUE)
Elec.ts <- ts(CBE[, 3], start = 1958, freq = 12)

plot(Elec.ts, xlab = "", ylab = "")
title(main = "Serie de Producción de Electricidad Australiana",
      ylab = "Producción de electricidad (GWh)",
      xlab = "Tiempo")# No estacionaria en la media, no distribuida en lo horizontal
#varainza no constante, tiene forma de embudo 

## Hay oscilaciones en la gráfica 

plot(diff(Elec.ts), xlab = "", ylab = "")#Diferencia (es una transformación) se hace constante constante en la media
#En series NO ESTACIONARIAS se puede transformar a estacionarias. Varianza NO constante porque parece EMBUDO
title(main = "Serie Diferenciada de Producción de Electricidad Australiana",
      xlab = "Tiempo", ylab = "Dif Serie",
      sub = "Gráfica de la serie diferenciada de primer órden")

plot(diff(log(Elec.ts)), xlab = "", ylab = "")##Hacer constante la varianza
title(main = "Serie de log dif de Producción de Electricidad Australiana",
      xlab = "Tiempo", ylab = "Dif log-Serie",
      sub = "Gráfica de la serie log-transformada diferenciada de primer órden")


######
## PACF y ACF MODELOS PARA PROPONER MODELOS cON CORRELOGRAMA
set.seed(1)## Esto simula una serie
x <- w <- rnorm(1000)
for(i in 3:1000) x[i] <- 0.5*x[i-1] + x[i-1] - 0.5*x[i-2] + w[i] + 0.3*w[i-1]## Modelo ARMA 1-1


plot(x, type = "l", 
     main = "Serie simulada de un modelo ARIMA(1, 1, 1)",
     xlab = "Tiempo",
     ylab = expression(x[t]),
     sub = expression(x[t] == 0.5*x[t-1] + x[t-1] - 0.5*x[t-2] + w[t] + 0.3*w[t-1]))

arima(x, order = c(1, 1, 1))## 1 = autoregresivo, 1 = parámetros móviles
## i = diferencia se aplica

#####
x <- arima.sim(model = list(order = c(1, 1, 1), ar = 0.5, ma = 0.3), n = 1000)
arima(x, order = c(1, 1, 1))## Ordenes del modelo, arima da la estimación de la serie de tiempo


plot(Elec.ts, xlab = "", ylab = "")
title(main = "Serie de Producción de Electricidad Australiana",
      ylab = "Producción de electricidad (GWh)",
      xlab = "Tiempo")

plot(log(Elec.ts), xlab = "", ylab = "")#Hace constante la VARAINZAa pero la media no
title(main = "Log de Serie de Producción de Electricidad Australiana",
      ylab = "Log de Producción de electricidad (GWh)",
      xlab = "Tiempo")


Elec.AR <- arima(log(Elec.ts), order = c(1, 1, 0), #numero de enmedio representa las diferencias para hacer constante la MEDIA
                 seas = list(order = c(1, 0, 0), 12))#MODELO CON UN PARAMETRO AUTORREGRESIVO, el 0 es que no hay promedios móviles

Elec.MA <- arima(log(Elec.ts), order = c(0, 1, 1),
                 seas = list(order = c(0, 0, 1), 12))


AIC(Elec.AR)#Estimación de autorregresión
AIC(Elec.MA)#Estimación de promedios móviles 

###¿Cómo se elige el mejor modelo?
## Criterio de AIC, el que tenga el menor valor, es el MEJOR MODELOS 

get.best.arima <- function(x.ts, maxord = c(1, 1, 1, 1, 1, 1)){
  best.aic <- 1e8
  n <- length(x.ts)
  for(p in 0:maxord[1])for(d in 0:maxord[2])for(q in 0:maxord[3])
    for(P in 0:maxord[4])for(D in 0:maxord[5])for(Q in 0:maxord[6])
    {
      fit <- arima(x.ts, order = c(p, d, q),
                   seas = list(order = c(P, D, Q),
                               frequency(x.ts)), method = "CSS")
      fit.aic <- -2*fit$loglik + (log(n) + 1)*length(fit$coef)
      if(fit.aic < best.aic){
        best.aic <- fit.aic
        best.fit <- fit
        best.model <- c(p, d, q, P, D, Q)
      }
    }
  list(best.aic, best.fit, best.model)
}

best.arima.elec <- get.best.arima(log(Elec.ts),
                                  maxord = c(2, 2, 2, 2, 2, 2))

best.fit.elec <- best.arima.elec[[2]]  # Modelo
best.arima.elec[[3]] # Tipo de modelo (órdenes)// ORDENES DE LOS MEJORES MODELOS
best.fit.elec## PARAMETROS ESTIMADOS DEL MEJOR MODELO
best.arima.elec[[1]] # AIC, EL MEJOR ENCONTRADO

acf(resid(best.fit.elec), main = "")## CORRELOGRAMA DE RUIDO BLANCO
title(main = "Correlograma de los residuales del ajuste")

pr <- predict(best.fit.elec, 12)$pred ## 12 Predicciones hacia el futuro
ts.plot(cbind(window(Elec.ts, start = 1981),## se grafica la serie de tiempo original 
              exp(pr)), col = c("blue", "red"), xlab = "")## Exponencial se hace a la inversa del logaritmo
title(main = "Predicción para la serie de producción de electricidad",
      xlab = "Mes",
      ylab = "Producción de electricidad (GWh)")

####### Sesión 7 ##########
