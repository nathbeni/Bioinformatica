# Carga los datos y exáminalos en R. Emplea las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? 2 ¿Cuántos tratamientos? 5

# Cargar los datos del archivo de texto
data <- read.table("/Users/nathbeni/Downloads/datos-trabajoR.txt", header = TRUE)

# Ver las primeras filas de los datos
head(data)

# Resumen de las estadisticas de los datos
summary(data)

# Numero de filas y columnas de los datos
dim(data)

# Estructura de los datos (tipo de datos y estructura)
str(data)

#Haz un boxplot para nuestros datos. Uno para cada variable. Colorea a Variable 1 y a Variable 2 de forma diferente (guarda esos colores para las siguientes gráficas)
#Hacemos un boxplot con las variables
boxplot(Variable1 ~ Tratamiento, data = data, col = "pink", main = "Boxplot de Variable 1 del Tratamiento")
boxplot(Variable2 ~ Tratamiento, data = data, col = "purple", main = "Boxplot de Variable 2 del Tratamiento")

#Haz un gráfico de dispersión con las dos variables. Cada tratamiento debe de ir de un color distinto. ¡Como en la siguiente imagen! Pista: usa col=datos$Tratamiento
#Creamos un gráfico de dispersión usando los colores del vector de tratamiento y dandole un título
plot(x = data$Variable1, y = data$Variable2, col=data$Tratamiento, 
     main = "Gráfico de dispersión", xlab = "Var1", ylab = "Var2")

#Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho. Pista: investiga sobre legend()
#Creamos una leyenda usando el vector de tratamientos eliminando duplicados
legend(x = "bottomright", legend = c(unique(data$Tratamiento)), title = "Tratamiento", fill = c(unique(data$Tratamiento)))

#Haz un histograma para cada variable. Recuerda mantener los colores.
#Creamos dos histogramas manteniendo los colores
hist(data$Variable1, col = "pink", main = "Var1")
hist(data$Variable2, col = "purple", main = "Var2")

#Haz un factor en la columna tratamiento y guárdalo en una variable. Pista: factor(factor$Tratamiento)
#Guardamos el factor en una variable
factor <- factor(data$Tratamiento)

#Calcula la media y la desviación estándar para cada tratamiento. Recomendación: es más fácil si usas aggregate() o tapply(). • aggregate(Variable~factor,datos,función) • tapply(datos$Variable,factor,función)
#Calculamos la media y la desviación típica usando aggregate (Adjuntada captura de pantalla de los resultados)
mediaVariable1 <- aggregate(Variable1~factor, data, mean)
desviacionTipicaVariable1 <- aggregate(Variable1~factor, data, sd)

mediaVariable2 <- aggregate(Variable2~factor, data, mean)
desviacionTipicaVariable2 <- aggregate(Variable2~factor, data, sd)

#Averigua cuántos elementos tiene cada tratamiento. Recomendación: es más fácil si usas table() con el factor
#Usando la función table vemos que cada tratamiento tiene 10 valores
table(factor)

#Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
#Usando la función subset extraemos los valores de la tabla 
dataTratamiento1 <- subset(data, Tratamiento == 1)
dataTratamiento4 <- subset(data, Tratamiento == 4)

print(dataTratamiento1)
print(dataTratamiento4)

#Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales. ¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal. En función del resultado de la prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus resultados.
print(dataTratamiento1)
print(dataTratamiento4)
tratamiento1 <- dataTratamiento1$Variable1
tratamiento4 <- dataTratamiento4$Variable1

#Prueba de Normalidad (Shapiro-Wilk)
shapiro_test_t1 <- shapiro.test(tratamiento1)
shapiro_test_t4 <- shapiro.test(tratamiento1)

#Prueba de Igualdad de Varianzas (Fisher)
var_test <- var.test(tratamiento1, tratamiento4)

#Prueba t de dos muestras (varianzas iguales)
t_test <- t.test(tratamiento1, tratamiento4, var.equal = TRUE)

# Imprimimos los resultados
#Resultados prueba de Normalidad (Shapiro-Wilk)
print(shapiro_test_t1)
print(shapiro_test_t4)
# Como el resultado es muy pequeño podemos decir que no sigue una distribución normal, entonces usariamos el test de varianzas

#Prueba de Igualdad de Varianzas
print(var_test)
# Como el resultado no es igual a 1, haremos el de dos muestras

#Prueba t de dos muestras (varianzas iguales)
print(t_test)

#Según los resultados podemos descartar la hipótesis de que las medias son iguales.
