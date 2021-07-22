# lncrna
Final Master's degree (Computational methods in Science)
# Introducción
El cáncer de hígado se ha convertido en el segundo tumor más letal de la actualidad, con una supervivencia a 5 años del 18% de los pacientes. El tipo de tumor primario más común es el hepatocarcinoma, suponiendo un 90% de todos los casos. 

Se suele dar en el contexto de una enfermedad hepática crónica. Esto es debido a la existencia de inflamación hepática que se sostiene en el tiempo, a la vez que aparece fibrosis y una regeneración defectuosa de los hepatocitos.

### Factores de riesgo
Las causas predominantes del hepatocarcinoma difieren según geografía. Mientras que en los países más desarrollados la causa primordial es la obesidad, en aquellos subdesarrollados lo es la infección por los virus de la hepatitis B (HBV) y hepatitis C (HCV).
Otros factores, como el tabaquismo, alcoholismo, etc. son considerados de riesgo para el desarrollo de la enfermedad.

La cirrosis es un factor que predispone a sufrir la enfermedad, si bien no es necesaria para el desarrollo del hepatocarcinoma y no está asociada a un peor pronóstico.

### Diagnóstico
Las técnicas diagnósticas son principalmente no invasivas, como el TAC, la ecografía, resonancia magnética... Sin embargo, es cada vez más común realizar biopsias para obtener una mejor caracterización molecular del tumor.

En la clínica, es común analizar los niveles de alfafetoproteína, proteína sintetizada en el hígado.
Los bebés recién nacidos tienen nivel alto de la alfafetoproteína, que se ve reducida al mínimo una vez transcurrido el primer año de vida. Así pues, puede ser de gran utilidad en el diagnóstico y seguimiento del cáncer combinado con otras pruebas. 

### Tratamiento

El hepatocarcinoma suele ser asintomático hasta sus fases más tardías, por lo que la batería de tratamientos disponibles disminuye considerablemente.

La resección hepática y el trasplante son los tratamientos más utilizados a día de hoy para tratar el HCC. En la actualidad, los tratamientos sistémicos son las terapias convencionales utilizadas en HCC, entre las que se incluyen los inhibidores del immune-checkpoint, inhibidores de la tirosina-quinasa y anticuerpos monoclonales. 

Sin embargo, la supervivencia de estos pacientes a 5 años es baja, y muestran tasas elevadas de recurrencia, por lo que es necesaria la investigación de nuevas terapias. 

### Cribados
Es conocido que la cirrosis incrementa el riesgo de desarrollar hepatocarcinoma. Por eso mismo, se ha realizado un ensayo para ver si un cribado a pacientes con perfil cirrótico  sería de utilidad para una detección precoz del cáncer, aumentando así su supervivencia. 
Los resultados arrojan que este cribado sí es de utilidad en pacientes cirróticos en clase funcional Child-Pugh A y B y ha de realizarse mediante ecografía con una frecuencia de unos 6 meses. Los sujetos en los que el trasplante está contraindicado no entrarían dentro del cribado, ya que una detección precoz no les beneficiaría.

### RNAs largos no codificantes (lncRNA)
Los RNAs largos no codificantes son secuencias de más de 200 nucleótidos que no son traducidas a proteína. 
En el pasado se desconocían sus funciones. De hecho, se llegó a creer que no servían para nada.
Sin embargo, se ha demostrado su importante papel en la regulación génica, como puede ser a nivel de transcripcional, post-transcripcional, epigenético y replicación del ADN, así como en la estabilidad cromosomal.

Cuantitativamente hablando, se expresan mucho menos que los RNA mensajeros (mRNA), pero pueden servir como potenciales biomarcadores, así como targets para nuevas terapias.



# Objetivos

### Objetivos generales
El objetivo principal es realizar una estratificación del hepatocarcinoma en base a los long non-coding RNA (lncRNA).

### Objetivos específicos
Clusterizar los genes y relacionarlos con una mayor o menor supervivencia en el paciente.

# Hipótesis
La desregulación de los lncRNA parece darse de forma más específica para cada enfermedad. Así pues, el estudio de estos puede ser útil para realizar mejores diagnósticos y pronósticos, así como para encontrar nuevas terapias.

# Materiales y metodología
### 1. Obtención de las muestras y control de calidad
Se utilizan los datos de RNA-seq de muestras tumorales de pacientes con hepatocarcinoma. Las muestras han sido obtenidas a través del proyecto NASIR, un ensayo clínico de medicamento realizado a nivel nacional, financiado por la farmacéutica BRISTOL-MYERS SQUIBB PHARMA EEIG y patrocinado por la Clínica Universidad de Navarra.
Una vez obtenidas las muestras, se les realiza un control de calidad.

### 2. Alineamiento y cuantificación de los tránscritos
#### 2.1 Con Kallisto
Se utiliza el programa de pseudoalineamientos Kallisto, ya que ofrece una mayor precisión que los métodos de alineamiento convencionales para los lncRNA ("Benchmark of long non-coding RNA cquantification for RNA sequencing of cancer samples"). A su vez, se ha combinado con una anotación transcriptómica completa ("Homo_sapiens.GRCh38.cdna.all.fa"), tal y como se recomienda. 

#### 2.2 Con Subread
En este caso se utiliza como transcriptoma de referencia el "gencode.v29.annotation.gtf". 

### 3. Normalización y transformación de los counts
Las cuentas se han normalizado mediante DESeq2, una librería de Bioconductor que realiza un de expresión génica diferencial basado en una distribución binomial negativa. A su vez, los lncRNA se expresan menos como norma general, por lo que son susceptibles de sufrir una mayor varianza. Para rectificar este fenómeno, se decide realizar la transformación rlog, disponible en la misma librería.

### 4. Filtrado de genes
Se aplican los siguientes filtros para las cuentas normalizadas:
- Cada gen debe tener un mínimo de 10 counts en al menos 5 pacientes.
- El máximo de cada gen tiene que ser igual o superior a 50
- La desviación estándar de cada gen debe ser igual o superior a 10.

### 5. Clustering de genes
#### 5.1 Clustering de las muestras completas
Se realiza un clustering jerárquico general de los genes que después se utiliza para hacer un análisis de supervivencia.

A continuación, se pretende dividir  en grupos de menos de 50 los genes, por lo que se realizan clusterings jerárquicos sucesivos hasta conseguir dicho propósito.
Para las clusterizaciones se utiliza la función "hclust" del paquete factoextra de R. La métrica de distancia de elección es la euclidea, y el método, "average". 
#### 5.2 Clustering de genes NICO
El grupo de investigación cuenta con una lista de 100 genes que son de especial interés (visto en proyectos anteriores). Por eso mismo, se hace una clusterización específica de estos genes.
Como el grupo de datos es pequeño, se prueba haciendo un clustering jerárquico general y, después, se dividen en clusters de 15 genes o menos.


### 6. Análisis de supervivencia
Se hace un clustering jerárquico de las muestras, dividiéndolas en dos grupos de pacientes (high- y low- risk).




# Resultados




#Discusión

#Conclusión

