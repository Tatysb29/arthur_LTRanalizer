# arthur_LTRanalizer
 
Arthur_LTRanalizer.
Es un  programa que permite  identificar y clasificar retrotransposones LTR en base a dominios proteicos conservados, usando bases de referencia de dominios de perfiles de REXdb  y GyDB. Para la identificación se usa la técnica de perfiles de Modelos Ocultos de Markov (HMM).
Requisitos de Software.
Se puede  usar una máquina virtual  de Linux en Ubuntu o un sistema Linux.
Dependencias.
Python 3.10 o mayor.
Paquete Anaconda.
Librerías de Python como: biopython, hmmer y xopen.

#Instalación.
Una vez que se terminó la aplicación, se creó un archivo ejecutable para Linux para facilitar la instalación al usuario,  de esta forma Arthur_LTRanalizer se puede instalar de dos maneras, la primera opción es directamente desde el ejecutable y la segunda opción es corriendo el script de phython directamente.
##Opción 1 en Linux a través del ejecutable.
1.	Descargar la carpeta “descargaConEjecutable”, esta carpeta contiene una subcarpeta Arthur_LTRAnalizer.
2.	En la terminal ingresar dentro de la carpeta Arthur_LTRAnalizer por  medio del comando cd.
3.	En la terminal añadir una ruta a la variable de entorno del PC Linux hacia la carpeta Arthur_LTRAnalizer con el comando: 
export PATH="$HOME/ miRutaHacialaCarpetaArthur:$PATH"
4.	Correr el programa con el comando: ./Arthur_LTRAnalizer 
Si aparece el error de permiso denegado, se debe habilitar el archivo como ejecutable, dando click derecho en el archivo Arthur_LTRAnalizer- click en propiedades y habilitar archivo ejecutable

##Opción 2 en Linux corriendo el programa .py directamente.
1.	Descargar la carpeta “descargaIndividual”.
2.	Instalar el paquete de Anaconda.
3.	En la terminal crear un ambiente de  trabajo  con el comando: conda create env miAmbiente.
4.	Ingresar al ambiente de trabajo creado con el comando: conda activate miAmbiente.
5.	Dentro de este ambiente configurar los canales de bioconda con los siguientes comandos en el mismo orden:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
6.	Instalar las siguientes librerías:  
biopython con: conda install -c conda-forge biopython # https://anaconda.org/conda-forge/biopython
xopen con: conda install -c bioconda xopen # https://anaconda.org/bioconda/xopen
hmmer con: conda install -c bioconda hmmer # https://anaconda.org/bioconda/hmmer
7.	En la  terminal, estando dentro del ambiente de trabajo elegido, ingresar dentro de la  carpeta Arthur_LTRAnalizer por medio del  comando cd.
8.	Se debe correr el programa con el comando: /Arthur_LTRAnalizer.py

