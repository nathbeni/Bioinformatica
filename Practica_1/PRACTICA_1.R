#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa

## Se requiere la entrega de un Jupyter Notebook con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
Dim(data) #esto mide las dimensiones
head(data)
tail(data)

# Hacemos un primer histograma para explorar los datos
hist(data)
hist(data, col = "gray", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data_log=log2(data) #Pasa de una distribución sesgada a una distribución normal. Para visualizar mejor los datos
hist(data_log)

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?

boxplot(data_log, col=c("blue","blue","blue","orange","orange","orange"), main="GSE5583 - boxplots", las=2)·
#Es un gráfico que muestra un analisis mas detallado de los datos respecto a la distribución de los datos de la muestra.
#Las2 sirve para poner los ejes en vertical, main sirve para el titulo, col para poner el color.

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?

hc = hclust(as.dist(1-cor(data_log))) #Si, es correcta la separación
plot(hc, main="GSE5583 - Hierachical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <- data[,1:3] #matrix array
ko <- data[,4:6] 
class(wt)
head(wt) #la coma es para escoger las columnas(wt siendo 1-3 y ko 4-6), head es para acceder a los datos de arriba

# Calcula las medias de las muestras para cada condición. Usa apply

wt.mean = apply(wt, 1, mean)
head(wt.mean)
ko.mean = apply(ko, 1, mean)
head(ko.mean)

# ¿Cuál es la media más alta?
max(wt.mean) #5375.53
max(ko.mean) #37460.5

#La media mas alta es la de Ko

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO", main = "GSE5583 - Scatter plot")
# Añadir una línea diagonal con abline
abline(0,1, col="red")
abline(h=2,col="pink")
abline(v=5,col="green")
# ¿Eres capaz de añadirle un grid?
grid()

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean)

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? Porque los datos no son fiables y los estoy manipulando.

# ¿Cuántas valores tiene cada muestra? 2 condiciones, 3 valores.  hacemos un t.test para cada fila y un pvalue para cada uno 12488 
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data)) { #Para cada gen
  x = wt[i,] # gene wt numero i
  y = ko[i,] # gene wt numero i
  
  #Hacemos el test
  t = t.test(x,y)
  
  #añadimos el p-value a la lista
  pvalue[i] = t$p.value
  #Añadimos las estadísticas a la lista
  tstat[i] = t$statistic
}

head(pvalue)
length(pvalue)
# Ahora comprobamos que hemos hecho TODOS los cálculos

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10? Se homogeniza mejor para las graficas.
hist(pvalue)
hist(-log10(pvalue),col = "gray")
# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean,-log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico? Los pvalues significativos estarán en la cola, pero el logaritmo lo invierte y estarian arriba.
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios? 426 6
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined, ]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]), col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? Si los valores del ko son mayores que wt no sale una expresión negativa entonces los sobre expresados estan en ese lugar y los reprimidos en el lado positivo
plot(diff.mean, -log(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0], -log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0], -log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")        
#diff.mean = wt.mean - ko.mean
# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap? Hemos agrupado por ko y wt y por otro lado por expresión
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, labRow=FALSE)

#cexcol es el tamaño de letra del eje x
#colV y Rowv son los dendograms


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap

# Lo guardamos en un archivo PDF

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t", quote = FALSE)
