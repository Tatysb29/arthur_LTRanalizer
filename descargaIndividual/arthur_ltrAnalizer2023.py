
#!/home/taty/anaconda3/envs/ambiente_tesis_PUCE/bin python3.10

# -*- coding: utf-8 -*-
#----Para crear el ambiente por primera vez hacer esto----
# 1.En la terminal de Conda Prompt pegar: para ver versiones disponibles de python en conda: conda search “^python$” 
# 2.En la terminal de  Conda Prompt pegar: para crear un ambiente: conda create -n ambiente_Tesis_BC_PUCE python=3.9.6 
# 3.En la terminal de Conda Prompt pegar: para ver todos los ambientes de trabajo creados: conda info -e
#----Para cargar programas en el nuevo ambiente cada vez hacer esto----
# 1.En la terminal de Conda Prompt pegar: para activar el ambiente de trabajo deseado: conda activate ambiente_Tesis_BC_PUCE
# 2.En la terminal de Conda Prompt pegar: para instalar paquetes en este ambiente: conda install -n ambiente_Tesis_BC_PUCE package
#----Para salir de la terminal
# 1.En la terminal de Conda Prompt pegar: para salir del ambiente: conda deactivate
#----Para eliminar un ambiente de python no de conda, que se creo por error
#En la terminal de Visual Studio: para eliminar un ambiente de pythhon y no conda: rm -r D:\PUCE_MAESTRIA\TesisApp\ambiente_TesisBC_PUCE2

#---Importo librerias--- 
from tkinter import*     #Para realizar interfaz grafica
from tkinter import ttk   #Para crear ventanas en interfaz grafica
from tkinter import messagebox   #Para mostrar mensajes en la interfaz grafica
from tkinter import font         #Para usar variedad de textos y tamanos en la interfaz grafica 
from tkinter import filedialog   #Para abrir la busqueda mipc de archivos
from tkinter import scrolledtext #Para usar cuadros de textos grandes de varias lineas con scroll para subir y bajas, como los que muestran las rutas en la parte derecha de la interfaz
#from typing import Dict, Type
import os                        #Para obtener nombres y rutas de archivos
import itertools
import logging   #Para usar el logging para mostrar mensajes del programa en la terminal
import sys
import re  #Para proporcionar expresiones regulares
import argparse
from io import IOBase
from Bio import SeqIO
from Bio.Data import CodonTable   #Para usar tabla de codon
from collections import OrderedDict, Counter  #Para usar un diccionario ordenado y counter para usar un diccionario con cantidades del    numero de ocurrencias
import subprocess
import multiprocessing #Para el multiproceso
#import uuid

class interfaz_grafica():
    def __init__(self,ventana):
        self.ventana=ventana
        ###------Creacion de ventana principal de la Interfaz grafica---
        #ventana=Tk()
        self.campos_interfaz_grafica()
    def campos_interfaz_grafica(self):
        self.ventana.title("ARTHUR_LTRanalizer VERSION 1.0") #titulo de ventana general
        self.ventana.geometry("1510x755+0+0") #tamanio ventana general
        self.ventana.resizable(True, True)     #Para poder modificar el tamano de pantalla por el usuario
        self.ventana.columnconfigure(0, weight=1)  #Ppara modificar tamano de ventana usando Sizegrip
        self.ventana.rowconfigure(0, weight=1)      #Para modificar tamano de ventana usando Sizegrip
    #ventana.resizable(0, 0)  

        ###--------Creacion del marco para titulo principal, bd es ancho del borde, pady es el espacio vertical desde arriba hasta inicio del texto
        self.marco_superior_titulo=Frame(self.ventana,bd=8,relief=RIDGE,bg="green",padx=2,pady=3) #marco dentro de ventana general, ridge con borde, padx espacio horiz desde borde hasta inicio de nuevo borde
        self.marco_superior_titulo.pack(fill="both", expand=1)     #marco se ubica en el orden que voy poniendo instrucciones, uso pack y no place para modificar tamano de pantalla por usuario, que llene en X
        self.marco_superior_titulo.config(height="60")             #Especifico el alto del marco
        ##---------Creacion del Marco Superior dentro de la Ventana General
        self.marco_superior=Frame(self.ventana,bd=8,relief=RIDGE,bg="green",padx=2,pady=7) #marco dentro de ventana general, ridge con borde, padx espacio horiz desde borde hasta inicio de nuevo borde
        self.marco_superior.pack(fill="both",expand=1)     #ubicacion para que pantalla sea modificable uso pack y no place que no permite modificar tamano de pantalla por usuario
        self.marco_superior.config(height=530)
        #---Creos subtitulos para ingresar datos descriptivos en la parte superior
        subtitulo_definicion=Label(self.marco_superior,font=("arial",14,"bold"),fg="white",bg="green",text="Arthur_LTRanalizer es una herramienta para identificar y clasificar retrotransposones LTR,",padx=2,pady=1)
        subtitulo_definicion.pack(fill="both",expand=1)
        subtitulo_definicion.config(height=1)
        #---Creos subtitulos para ingresar datos descriptivos en la parte superior
        subtitulo_definicion_dos=Label(self.marco_superior,font=("arial",14,"bold"),fg="white",bg="green",text="a nivel de linaje de secuencias de plantas.",padx=2,pady=0)
        subtitulo_definicion_dos.pack(fill="both",expand=1)
        subtitulo_definicion_dos.config(height=1)
        #---Creos subtitulos para ingresar datos descriptivos en la parte superior
        subtitulo_definicion_tres=Label(self.marco_superior,font=("arial",12,"bold"),fg="white",bg="green",text="Nota: La búsqueda se realiza por default con la base de referencia de dominios proteicos REXdb.  El archivo  de entrada debe estar en formato fasta.",padx=2,pady=7)
        subtitulo_definicion_tres.pack(fill="both",expand=1)
        subtitulo_definicion_tres.config(height=1)
        ##--------Creo Submarco Izquierdo, dentro del marco superior para introducir datos personales, labelFrame introduce un titulo incrustado en el marco
        self.submarco_left_superior=Frame(self.marco_superior,bd=8,relief=RIDGE,padx=2)
        self.submarco_left_superior.pack(fill="both",expand=1)
        self.submarco_left_superior.config(height=415)
        ##--------Creo el Sub-submarco Inferior para Boton:escoger archivo, escoger base de referencia e iniciar proceso
        self.subsubmarco_inferior_left_superior=Frame(self.submarco_left_superior,bd=3,relief=RIDGE,padx=2)
        self.subsubmarco_inferior_left_superior.pack(fill="both",expand=1)
        self.subsubmarco_inferior_left_superior.config(height=395) #indico alto del submarco izq
        ##--------Creo el Marco para colocar PUCE
        self.marco_botones=Frame(self.ventana,bd=8,relief=RIDGE,bg="green",padx=2)
        self.marco_botones.pack(fill="both",expand=1)
        self.marco_botones.config(height=47)
        #_________Creo el Marco Inferior donde va a ir en la parte inferior con autor y tutor
        self.marco_inferior=Frame(self.ventana,bd=8,relief=RIDGE,padx=2)
        self.marco_inferior.pack(fill="both",expand=1)
        self.marco_inferior.config(height=90)

        #----Definicion de variables de ventana de interfaz grafica----
        self.nombre_direccion_archivo=()
        self.nombre_direccion_salida=()
        self.nombre_archivo_entrada=()
        #coloco_nombre_archivo_entrada=()
        #coloco_nombre_ruta_salida=()

        self.inicio=()
        self.varmenu_basereferencia_dominios=StringVar()   #variable que almacena la opcion escogida del menu desplegable de archivos de referencia de dominios proteicos
        self.opcionesmenu_basereferencia_dominios=["REXdb","Gydb"] #lista de opciones que se muestran en la interfaz para que usuario escoja el archivo de referencia de dominios proteicos contra el que se va a analizar
        self.varmenu_basereferencia_dominios.set(self.opcionesmenu_basereferencia_dominios[0])   #indico que la opcion Ninguno aparezca al inicio, ocupa la posicion 0 de la lista
        self.varmenu_tipo_secuencia=StringVar()   #variable que almacena la opcion escogida del menu desplegable del tipo de secuencia genomica o TE
        self.opcionesmenu_tipo_secuencia=["TE","Genomicas"] #lista de opciones que se muestran en la interfaz para que usuario escoja el tipo de secuencia genomica o TE 
        self.varmenu_tipo_secuencia.set(self.opcionesmenu_tipo_secuencia[0])   #indico que la opcion Ninguno aparezca al inicio, ocupa la posicion 0 de la lista

        self.varmenu_tipo_secuencia_NuclProt=StringVar()   #variable que almacena la opcion escogida del menu desplegable del tipo de secuencia de  entrada Nucleotidos o Proteinas
        self.opcionesmenu_tipo_secuencia_NuclProt=["Nucleotidos","Proteinas"] #lista de opciones que se muestran en la interfaz para que usuario escoja el tipo de secuencia de entrada Nucleotidos o Proteinas 
        self.varmenu_tipo_secuencia_NuclProt.set(self.opcionesmenu_tipo_secuencia_NuclProt[0])   #indico que la opcion Ninguno aparezca al inicio, ocupa la posicion 0 de la lista


        #----Creo titulo prinicpal
        titulo_principal=Label(self.marco_superior_titulo,font=("arial",17,"bold"),fg="white",bg="green",text="BIENVENIDO a Arthur_LTRanalizer.",padx=2,pady=1)
        titulo_principal.pack(fill="both",expand=1)
        titulo_principal.config(height=1)
        #----Subtitulo para escoger el archivo de secuencias de entrada
        subtitulo_escoger_secuencias=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="1. Escoja el archivo con las secuencias a evaluar ",padx=2,pady=15)
        subtitulo_escoger_secuencias.grid(row=0,column=0,sticky="w")
        #---Subtitulo vacio
        subtitulo_escoger_secuencias_vacio=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="  ",padx=2,pady=3)
        subtitulo_escoger_secuencias_vacio.grid(row=0,column=1,sticky="w")
        #----Subtitulo para escoger el tipo de secuencia si es de genoma o elementos transponibles
        subtitulo_escoger_tiposecuencia=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="2. Escoja el tipo de secuencia de entrada",padx=2,pady=1)
        subtitulo_escoger_tiposecuencia.grid(row=1,column=0,sticky="w")
        #----Subtitulo para escoger el tipo de secuencia de entrada si son nucleotidos o proteinas
        subtitulo_escoger_tiposecuencia_nucleotoProteinas=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="  ",padx=2,pady=5)
        subtitulo_escoger_tiposecuencia_nucleotoProteinas.grid(row=2,column=0,sticky="w")
        #----Subtitulo para escoger el archivo de referencia de dominios proteicos
        subtitulo_escoger_archivo_dominios=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="3. Escoja la base de referencia de dominios proteicos contra la que se analizará",padx=2,pady=35)
        subtitulo_escoger_archivo_dominios.grid(row=3,column=0,sticky="w")
        #----Subtitulo para boton iniciar proceso
        subtitulo_iniciar_proceso=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="4. Escoja la ruta de salida",padx=2,pady=1)
        subtitulo_iniciar_proceso.grid(row=4,column=0,sticky="w")
        #----Subtitulo para archivos de salida
        subtitulo_archivo_salida=Label(self.subsubmarco_inferior_left_superior,font=("arial",14,"bold"),text="5. De click en el boton Iniciar",padx=2,pady=28)
        subtitulo_archivo_salida.grid(row=5,column=0,sticky="w")
        #----Creo titulo final PUCE
        titulo_final=Label(self.marco_botones,font=("arial",17,"bold"),fg="white",bg="green",text="PONTIFICIA UNIVERSIDAD CATOLICA DEL ECUADOR",padx=448,pady=2)
        titulo_final.grid(row=0,column=0)
        #----Creo titulo final autor
        titulo_final_autor=Label(self.marco_inferior,font=("arial",13,"bold"),fg="green",text="Autor: Tatiana Paola Benalcazar Vayas",padx=593,pady=2)
        titulo_final_autor.grid(row=0,column=0,sticky=(N, S, E, W),pady=1)
        #titulo_final_autor.config(height=1)
        #----Creo titulo final director
        titulo_final_director=Label(self.marco_inferior,font=("arial",13,"bold"),fg="green",text="Director: Romain Guyot",padx=50,pady=2)
        titulo_final_director.grid(row=1,column=0,sticky=(N, S, E, W),pady=1)
        #titulo_final_director.config(height=1)
        #----Creo titulo final citar
        titulo_final_citar=Label(self.marco_inferior,font=("arial",11,"bold"),fg="green",text="Citar: Tatiana Paola Benalcazar Vayas",padx=45,pady=2)
        titulo_final_citar.grid(row=2,column=0,sticky=(N, S, E, W),pady=1)
        #titulo_final_citar.config(height=1)
        #-----Creo subtitulos para colocar nombre y direccion de archivo de entrada
        #coloco_nombre_archivo_entrada=Entry(subsubmarco_inferior_left_superior,font=("arial",11,"bold"),bd=0,borderwidth=0,bg="white",relief=RIDGE,width=53,textvariable=nombre_archivo_entrada)
        #coloco_nombre_archivo_entrada.grid(row=0,column=3,sticky="w")
        #coloco_nombre_archivo_entrada.configure(state="disabled")
        self.coloco_nombre_archivo_entrada=scrolledtext.ScrolledText(self.subsubmarco_inferior_left_superior,width=35,height=6,wrap=WORD,bd=0,borderwidth=0,relief=FLAT,bg="gray94",fg="gray40",font=("Arial",11,"bold"),padx=10)
        self.coloco_nombre_archivo_entrada.grid(row=0,column=3,sticky="w")
        self.coloco_nombre_archivo_entrada.configure(state="disabled")
        #-----Creo subtitulos para colocar la ruta de salida en donde se almacenaran loa archivos de salida
        #coloco_nombre_ruta_salida=Entry(subsubmarco_inferior_left_superior,font=("arial",11,"bold"),bd=0,borderwidth=0,bg="white",relief=RIDGE,width=53,textvariable=nombre_direccion_salida)
        #coloco_nombre_ruta_salida.grid(row=3,column=3,sticky="w")
        #coloco_nombre_ruta_salida.configure(state="disabled")
        self.coloco_nombre_ruta_salida=scrolledtext.ScrolledText(self.subsubmarco_inferior_left_superior,width=35,height=3,wrap=WORD,bd=0,borderwidth=0,relief=FLAT,bg="gray94",fg="gray40",font=("Arial",11,"bold"),pady=1,padx=10)
        self.coloco_nombre_ruta_salida.grid(row=4,column=3,sticky="w")
        self.coloco_nombre_ruta_salida.configure(state="disabled")

        #----Definicion de variables para identificar transposones LTR----
        #---Obtener ruta donde esta almacenado el programa---
        ruta_programa=os.path.dirname(os.path.realpath(__file__))  #Se encuentra la ruta donde esta guardado el programa
        print("Esta es la ruta del programa: "+str(ruta_programa))
        #---Se crea un diccionario con los nombres de bases de datos de referencia  de dominios proteicos contra las que se va a analizar
        lista_bdRef_dominios={
            "bdRef_REXdb_plantas":ruta_programa+"/REXdb_protein_database_viridiplantae_v3.0.hmm",  #Coloca la base de datos de referencia que voy a usar que es REXdb
            "bdRef_Gydb":ruta_programa+"/GyDB2.hmm",  #Coloca la base de datos de referencia que voy a usar que es Gydb
        }
        self.nombre_ruta_Ref_REXdb=lista_bdRef_dominios["bdRef_REXdb_plantas"]  #Nombre y ruta de base de datos de dominios proteicos de REXdb
        self.nombre_ruta_Ref_Gydb=lista_bdRef_dominios["bdRef_Gydb"]  #Nombre y ruta de base de datos de dominios proteicos de Gydb

        ##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Escoger archivo de secuencias
        boton_escoger_archivo_entrada=Button(self.subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="ESCOGER ARCHIVO",padx=5,pady=2,width=17,command=self.escoger_archivo_entrada)
        boton_escoger_archivo_entrada.grid(row=0,column=2,sticky="w")
        ##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger tipo de secuencia
        self.tipo_secuencia_menu=OptionMenu(self.subsubmarco_inferior_left_superior,self.varmenu_tipo_secuencia,*self.opcionesmenu_tipo_secuencia) 
        self.tipo_secuencia_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
        self.tipo_secuencia_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
        self.tipo_secuencia_menu.grid(row=1,column=2,padx=0,pady=1,sticky=W)
        ##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger tipo de secuencia Nucleotidos o Proteinas
        self.tipo_secuencia_NuclProt_menu=OptionMenu(self.subsubmarco_inferior_left_superior,self.varmenu_tipo_secuencia_NuclProt,*self.opcionesmenu_tipo_secuencia_NuclProt) 
        self.tipo_secuencia_NuclProt_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
        self.tipo_secuencia_NuclProt_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
        self.tipo_secuencia_NuclProt_menu.grid(row=2,column=2,padx=0,pady=1,sticky=W)
        ##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger archivo de referencia si es REXdb o Gydb
        self.baseref_dominios_menu=OptionMenu(self.subsubmarco_inferior_left_superior,self.varmenu_basereferencia_dominios,*self.opcionesmenu_basereferencia_dominios) 
        self.baseref_dominios_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
        self.baseref_dominios_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
        self.baseref_dominios_menu.grid(row=3,column=2,padx=0,pady=1,sticky=W)
        ##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Inciar el proceso
        boton_inicio=Button(self.subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="INICIAR BUSQUEDA",padx=5,pady=2,width=17,command=self.iniciar_proceso)
        boton_inicio.grid(row=5,column=2,sticky="w")
        ##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Escoger ruta de salida en donde se guardaran loa archivos de salida
        boton_escoger_archivo_salida=Button(self.subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="ESCOGER RUTA",padx=5,pady=2,width=17,command=self.escoger_ruta_salida)
        boton_escoger_archivo_salida.grid(row=4,column=2,sticky="w")

        #----Configuracion para poder cambiar tamano de pantalla usando grid
        self.subsubmarco_inferior_left_superior.columnconfigure(0, weight=1)
        self.subsubmarco_inferior_left_superior.columnconfigure(1, weight=1)
        self.subsubmarco_inferior_left_superior.columnconfigure(2, weight=1)
        self.subsubmarco_inferior_left_superior.columnconfigure(3, weight=1)
        self.subsubmarco_inferior_left_superior.rowconfigure(0, weight=1)
        self.subsubmarco_inferior_left_superior.rowconfigure(1, weight=1)
        self.subsubmarco_inferior_left_superior.rowconfigure(2, weight=1)
        self.subsubmarco_inferior_left_superior.rowconfigure(3, weight=1)
        self.subsubmarco_inferior_left_superior.rowconfigure(4, weight=1)
        self.submarco_left_superior.rowconfigure(0, weight=1)
        self.submarco_left_superior.columnconfigure(0, weight=1)
        self.marco_inferior.rowconfigure(0, weight=1)
        self.marco_inferior.rowconfigure(1, weight=1)
        self.marco_inferior.rowconfigure(2, weight=1)
        self.marco_inferior.columnconfigure(0, weight=1)
        self.marco_botones.rowconfigure(0, weight=1)
        self.marco_botones.columnconfigure(0, weight=1)
        self.marco_superior.rowconfigure(0, weight=1)
        self.marco_superior.columnconfigure(0, weight=1)
        self.ventana.rowconfigure(0, weight=1)

        #------Para poner 4 puntos en la parte inferior derecha de la pantalla para poder modificar tamano de pantalla
        puntos_extender_esquina = ttk.Sizegrip(self.ventana)   #Coloco puntos para  arrastrar diagonalmente y modificar tamano de pantalla
        puntos_extender_esquina.pack(expand = True, fill = BOTH, anchor = SE)


    #----Funcion para escoger archivo de entrada----
    def escoger_archivo_entrada(self):
        #global nombre_direccion_archivo #le hice global porque esta variable almacena info que se quiere usar en otras funciones
        #global nombre_archivo_entrada
        print("Nombre direccion archivo: "+str(self.nombre_direccion_archivo))
        self.nombre_archivo_entrada=()
        #-----Creo subtitulos para colocar nombre y direccion de archivo de entrada
        #coloco_nombre_archivo_entrada=Entry(subsubmarco_inferior_left_superior,font=("arial",11,"bold"),bd=0,borderwidth=0,bg="white",relief=RIDGE,width=53,textvariable=nombre_archivo_entrada)
        #coloco_nombre_archivo_entrada.grid(row=0,column=3,sticky="w")
        #coloco_nombre_archivo_entrada.configure(state="disabled")
        self.coloco_nombre_archivo_entrada=scrolledtext.ScrolledText(self.subsubmarco_inferior_left_superior,width=35,height=6,wrap=WORD,bd=0,borderwidth=0,relief=FLAT,bg="gray94",fg="gray40",font=("Arial",11,"bold"),padx=10)
        self.coloco_nombre_archivo_entrada.grid(row=0,column=3,sticky="w")
        self.coloco_nombre_archivo_entrada.configure(state="disabled")
        try:
            self.nombre_direccion_archivo=filedialog.askopenfilename(title="Escoger Archivo")
        except:
            print("Ingreso por excepcion, no escogio archhivo de entrada")
        try:   # en caso de que no se escoja el archivo de entrada, evito que programa se bloqueo
            self.nombre_archivo_entrada=os.path.basename(self.nombre_direccion_archivo)
        except: #Si se produce la excepcion o el error TypeError entro por este lado y evito que el programa se bloquee
            print("Se produjo  error porque no se escogio archivo de entrada") 
        self.coloco_nombre_archivo_entrada.config(state="normal")
        ##----Borro datos en cuadros de texto de la derecha
        #coloco_nombre_archivo_entrada.delete(0,END)
        self.coloco_nombre_archivo_entrada.delete("0.0",END)
        ##----Ingreso datos personales en cuadros de texto de la derecha
        self.coloco_nombre_archivo_entrada.insert(END,"RUTA:"+"\n")
        #coloco_nombre_archivo_entrada.insert(END,"""
        #""")
        self.coloco_nombre_archivo_entrada.insert(INSERT,str(self.nombre_direccion_archivo))
        self.coloco_nombre_archivo_entrada.insert(END,""+"\n")
        self.coloco_nombre_archivo_entrada.insert(END,"NOMBRE DE ARCHIVO:"+"\n")
        self.coloco_nombre_archivo_entrada.insert(INSERT,str(self.nombre_archivo_entrada))
        self.coloco_nombre_archivo_entrada.configure(state="disabled")

    #----Funcion para escoger ruta de salida----
    def escoger_ruta_salida(self):
        #global nombre_direccion_salida #le hice global porque esta variable almacena info que se quiere usar en otras funciones
        #global coloco_nombre_ruta_salida
        self.nombre_direccion_salida=filedialog.askdirectory(title="Escoger Ruta de Salida")
        self.coloco_nombre_ruta_salida.config(state="normal")
        ##----Borro datos en cuadros de texto de la derecha
        #coloco_nombre_ruta_salida.delete(0,END)
        ##----Ingreso datos en cuadros de texto de la derecha
        #coloco_nombre_ruta_salida.insert(END,str(nombre_direccion_salida))
        #coloco_nombre_ruta_salida.configure(state="disabled")
        self.coloco_nombre_ruta_salida.delete("0.0",END)
        self.coloco_nombre_ruta_salida.insert(INSERT,str(self.nombre_direccion_salida))
        self.coloco_nombre_ruta_salida.configure(state="disabled")
    
    #---Funcion para encerar valores
    def colocar_parametros_vacios(self):
        #global nombre_direccion_archivo
        #global nombre_direccion_salida
        #global inicio
        self.nombre_archivo_entrada=()
        #global coloco_nombre_archivo_entrada
        self.nombre_direccion_archivo=""  #Coloco el nombre vacio para que pueda iniciarse una nueva busqueda
        self.nombre_direccion_salida=""  #Coloco el nombre vacio para que pueda iniciarse una nueva busqueda
        self.inicio="" # Coloco el nombre vacio para que pueda inicarse una nueva busqueda
        print("Estoy dentro de colocar vacios")
        ##----Borro datos en cuadros de texto de la derecha
        self.coloco_nombre_archivo_entrada.config(state="normal")
        self.coloco_nombre_ruta_salida.config(state="normal")
        #coloco_nombre_archivo_entrada.delete(0,END)
        self.coloco_nombre_archivo_entrada.delete("0.0",END)
        #coloco_nombre_ruta_salida.delete(0,END)
        self.coloco_nombre_ruta_salida.delete("0.0",END)
        self.coloco_nombre_ruta_salida.configure(state="disabled")
        self.coloco_nombre_archivo_entrada.configure(state="disabled")

    #---Funcion para dar inicio a la busqueda----
    def iniciar_proceso(self):
        #global nombre_direccion_archivo
        #global nombre_direccion_salida
        #global nombre_archivo_entrada
        #global inicio

        if len(self.nombre_direccion_archivo)!=0 and len(self.nombre_direccion_salida)!=0:
            self.inicio="OK"
            print("Inicio OK")
        else:
            messagebox.showerror(title="LOS CAMPOS NO PUEDEN ESTAR VACIOS",message="ESCOJA TODOS LOS CAMPOS."+"\n"+"\n"+str("El campo de archivo de entrada y/o ruta de salida no se ha definido."))
        print("Este es el nombre y direccion del archivo de secuencias:"+str(self.nombre_direccion_archivo))
        print("Este es el nombre del archivo de secuencias:"+str(self.nombre_archivo_entrada))
        print("Esta es la ruta de salida para archivos:"+str(self.nombre_direccion_salida))
        base_referencia_escogida=self.varmenu_basereferencia_dominios.get()
        tipo_secuencia_escogida=self.varmenu_tipo_secuencia.get()
        tipo_secuencia_escogida_NuclProt=self.varmenu_tipo_secuencia_NuclProt.get()
        print("Este es la base de referencia escogida:"+str(base_referencia_escogida))
        print("Este es el tipo de secuencia escogida:"+str(tipo_secuencia_escogida))
        print("Este es el tipo de secuencia escogida Nucleotidos  o Proteinas:"+str(tipo_secuencia_escogida_NuclProt))
        #---Creo carpeta Archivos_Arthur_LTR en la ruta seleccionada anteriormente
        if self.inicio=="OK":
            ruta_nuevacarpeta_salida=self.nombre_direccion_salida+str("/Archivos_Arthur_LTR")
            print("Esta es la ruta de nueva carpeta creada:"+str(ruta_nuevacarpeta_salida))
            respuesta=messagebox.askokcancel(title="Inicio del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR INICIARÁ."+"\n"+"\n"+str("El archivo de secuencias que se evaluará es:")+"\n"+str(self.nombre_archivo_entrada)+"\n"+"\n"+str("Los archivos de salida se ubicarán en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK para dar inicio y este proceso tardará unos minutos"))
            if respuesta==True:
                os.makedirs(ruta_nuevacarpeta_salida, exist_ok=True)  #Crea la carpeta
                print("La carpeta se creó exitosamente")
                #---Si se va analizar secuencias de TE
                if tipo_secuencia_escogida==self.opcionesmenu_tipo_secuencia[0]:
                    print("El tipo de secuencias es: "+str(self.opcionesmenu_tipo_secuencia[0]))
                    #---Si la base de referencia escogida es REXdb
                    if base_referencia_escogida==self.opcionesmenu_basereferencia_dominios[0]: 
                        print("La base de referencia es: "+str(self.opcionesmenu_basereferencia_dominios[0]))
                        #---Verificar version de programas dependientes
                        comprobar_version_hmmer_bdRef(self.nombre_ruta_Ref_REXdb,nombre_programa='hmmscan') #LLamo a la funcion para comprobar que versiones del programa hmmscan y de base ref Rexdb tengan igual version
                        compresion_indexamiento_bd_ref(self.nombre_ruta_Ref_REXdb)
                        longitud_max_secuencia2=verificar_tamanio_secuencia(self.nombre_direccion_archivo)

                        if longitud_max_secuencia2<1e6 and tipo_secuencia_escogida=="Genomicas":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia2)+str(". Son secuencias de elementos transponibles y  no Genoma. Por favor escoja opcion TE"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de elementos transponibles y no de Genoma. Por favor escoja la opcion TE.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            #respuesta1=messagebox.askokcancel(title="Inicio del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR INICIARÁ."+"\n"+"\n"+str("El archivo de secuencias que se evaluará es:")+"\n"+str(nombre_archivo_entrada)+"\n"+"\n"+str("Los archivos de salida se ubicarán en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK para dar inicio y este proceso tardará unos minutos"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif longitud_max_secuencia2>1e6 and tipo_secuencia_escogida=="TE":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia2)+str(". Es un genoma y no secuencia de elementos transponibles. Por favor escoja opcion Genomicas"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de un genoma y no de elementos transponibles. Por favor escoja la opcion Genomicas.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif (longitud_max_secuencia2>1e6 and tipo_secuencia_escogida=="Genomicas") or (longitud_max_secuencia2<1e6 and tipo_secuencia_escogida=="TE"):
                            #---Escanear con perfiles de modelos ocultos de markov
                            lista_archivos_traducidos2=escanear_hmm(self.nombre_direccion_archivo,ruta_nuevacarpeta_salida,self.nombre_ruta_Ref_REXdb,base_referencia_escogida,self.nombre_archivo_entrada,tipo_secuencia_escogida_NuclProt)
                            #----Proceso de clasificacion
                            palabra4=base_referencia_escogida
                            nombre_ruta_archivo_unif_tabladom2='%s/%s.%s.tabladom' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra4)
                            nombre_ruta_archivo_gff2,nombre_ruta_archivo_seq2=clasificar(lista_archivos_traducidos2,nombre_ruta_archivo_unif_tabladom2,tipo_secuencia_escogida,base_referencia_escogida,self.nombre_direccion_archivo,tipo_secuencia_escogida_NuclProt,self.nombre_ruta_Ref_Gydb,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada)
                            nombre_ruta_archivo_clasificacion_tsv='%s/%s.%s.clasificacion.tsv' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra4)
                            identificacion_archivo_clasificaion_tsv=open(nombre_ruta_archivo_clasificacion_tsv,'w')
                            listas_elementos_linea_clasificados=clasificar_segundo_paso(nombre_ruta_archivo_gff2,base_referencia_escogida,identificacion_archivo_clasificaion_tsv,self.nombre_ruta_Ref_Gydb)
                            clasificar_tercer_paso(listas_elementos_linea_clasificados,identificacion_archivo_clasificaion_tsv,nombre_ruta_archivo_seq2,nombre_ruta_archivo_clasificacion_tsv,nombre_ruta_archivo_gff2,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,base_referencia_escogida,self.nombre_direccion_archivo)
                            messagebox.showinfo(title="Finalización del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR HA FINALIZADO EXITÓSAMENTE."+"\n"+"\n"+str("Los archivos de salida se ubican en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK."))
       
                    #--Si base de referencia escogida es Gydb
                    if base_referencia_escogida==self.opcionesmenu_basereferencia_dominios[1]:
                        print("La base de referencia escogida es: "+str(self.opcionesmenu_basereferencia_dominios[1]))
                        comprobar_version_hmmer_bdRef(self.nombre_ruta_Ref_Gydb,nombre_programa='hmmscan') #LLamo a la funcion para comprobar que versiones del programa hmmscan y de base ref Rexdb tengan igual version
                        compresion_indexamiento_bd_ref(self.nombre_ruta_Ref_Gydb)
                        longitud_max_secuencia4=verificar_tamanio_secuencia(self.nombre_direccion_archivo)
                        if longitud_max_secuencia4<1e6 and tipo_secuencia_escogida=="Genomicas":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia4)+str(". Son secuencias de elementos transponibles y  no Genoma. Por favor escoja opcion TE"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de elementos transponibles y no de Genoma. Por favor escoja la opcion TE.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            #respuesta1=messagebox.askokcancel(title="Inicio del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR INICIARÁ."+"\n"+"\n"+str("El archivo de secuencias que se evaluará es:")+"\n"+str(nombre_archivo_entrada)+"\n"+"\n"+str("Los archivos de salida se ubicarán en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK para dar inicio y este proceso tardará unos minutos"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif longitud_max_secuencia4>1e6 and tipo_secuencia_escogida=="TE":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia4)+str(". Es un genoma y no secuencia de elementos transponibles. Por favor escoja opcion Genomicas"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de un genoma y no de elementos transponibles. Por favor escoja la opcion Genomicas.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif (longitud_max_secuencia4>1e6 and tipo_secuencia_escogida=="Genomicas") or (longitud_max_secuencia4<1e6 and tipo_secuencia_escogida=="TE"): 
                            #---Escanear con perfiles de modelos ocultos de markov
                            lista_archivos_traducidos2=escanear_hmm(self.nombre_direccion_archivo,ruta_nuevacarpeta_salida,self.nombre_ruta_Ref_Gydb,base_referencia_escogida,self.nombre_archivo_entrada,tipo_secuencia_escogida_NuclProt)
                            #----Proceso de clasificacion
                            palabra4=base_referencia_escogida
                            nombre_ruta_archivo_unif_tabladom2='%s/%s.%s.tabladom' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra4)
                            nombre_ruta_archivo_gff2,nombre_ruta_archivo_seq2=clasificar(lista_archivos_traducidos2,nombre_ruta_archivo_unif_tabladom2,tipo_secuencia_escogida,base_referencia_escogida,self.nombre_direccion_archivo,tipo_secuencia_escogida_NuclProt,self.nombre_ruta_Ref_Gydb,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada)
                            nombre_ruta_archivo_clasificacion_tsv='%s/%s.%s.clasificacion.tsv' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra4)
                            identificacion_archivo_clasificaion_tsv=open(nombre_ruta_archivo_clasificacion_tsv,'w')
                            listas_elementos_linea_clasificados=clasificar_segundo_paso(nombre_ruta_archivo_gff2,base_referencia_escogida,identificacion_archivo_clasificaion_tsv,self.nombre_ruta_Ref_Gydb)
                            clasificar_tercer_paso(listas_elementos_linea_clasificados,identificacion_archivo_clasificaion_tsv,nombre_ruta_archivo_seq2,nombre_ruta_archivo_clasificacion_tsv,nombre_ruta_archivo_gff2,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,base_referencia_escogida,self.nombre_direccion_archivo)
                            messagebox.showinfo(title="Finalización del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR HA FINALIZADO EXITÓSAMENTE."+"\n"+"\n"+str("Los archivos de salida se ubican en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK."))
            
                #---Si se va analizar a partir de Secuencias de Genoma
                if tipo_secuencia_escogida==self.opcionesmenu_tipo_secuencia[1]:
                    print("El tipo de secuencias es: "+str(self.opcionesmenu_tipo_secuencia[1]))
                    #---Si la base de referencia escogida es REXdb
                    if base_referencia_escogida==self.opcionesmenu_basereferencia_dominios[0]: 
                        print("La base de referencia es: "+str(self.opcionesmenu_basereferencia_dominios[0]))
                        comprobar_version_hmmer_bdRef(self.nombre_ruta_Ref_REXdb,nombre_programa='hmmscan') #LLamo a la funcion para comprobar que versiones del programa hmmscan y de base ref Rexdb tengan igual version
                        compresion_indexamiento_bd_ref(self.nombre_ruta_Ref_REXdb)
                        longitud_max_secuencia3=verificar_tamanio_secuencia(self.nombre_direccion_archivo)
                        if longitud_max_secuencia3<1e6 and tipo_secuencia_escogida=="Genomicas":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia3)+str(". Son secuencias de elementos transponibles y  no Genoma. Por favor escoja opcion TE"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de elementos transponibles y no de Genoma. Por favor escoja la opcion TE.")+"\n"+"\n"+str("De click en OK para cancelar el proceso"))
                            #respuesta1=messagebox.askokcancel(title="Inicio del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR INICIARÁ."+"\n"+"\n"+str("El archivo de secuencias que se evaluará es:")+"\n"+str(nombre_archivo_entrada)+"\n"+"\n"+str("Los archivos de salida se ubicarán en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK para dar inicio y este proceso tardará unos minutos"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif longitud_max_secuencia3>1e6 and tipo_secuencia_escogida=="TE":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia3)+str(". Es un genoma y no secuencia de elementos transponibles. Por favor escoja opcion Genomicas"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de un genoma y no de elementos transponibles. Por favor escoja la opcion Genomicas.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif (longitud_max_secuencia3>1e6 and tipo_secuencia_escogida=="Genomicas") or (longitud_max_secuencia3<1e6 and tipo_secuencia_escogida=="TE"):
                            palabra5=base_referencia_escogida
                            nombre_ruta_secuencias_recortadas='%s/secuencias_recortadas.%s.%s' % (ruta_nuevacarpeta_salida,palabra5,self.nombre_archivo_entrada)
                            #nombre_ruta_secuencias_recortads='{}/cu.{'.format(prueba44,'fasta')
                            print("Dentro se secuencias genomicas, esto es nombre ruta secuencias recortadas: "+str(nombre_ruta_secuencias_recortadas))
                            with open(nombre_ruta_secuencias_recortadas,'w') as identificacion_archivo_secuen_recortadas:
                                recortar_secuencias_genomicas(self.nombre_direccion_archivo,identificacion_archivo_secuen_recortadas)
                            #---Escanear con perfiles de modelos ocultos de markov
                            lista_archivos_traducidos2=escanear_hmm(nombre_ruta_secuencias_recortadas,ruta_nuevacarpeta_salida,self.nombre_ruta_Ref_REXdb,base_referencia_escogida,self.nombre_archivo_entrada,tipo_secuencia_escogida_NuclProt)
                            print("despues de escanear hmm")
                            #----Proceso de clasificacion
                            palabra6=base_referencia_escogida
                            nombre_ruta_archivo_unif_tabladom2='%s/%s.%s.tabladom' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra6)
                            nombre_ruta_archivo_gff2,nombre_ruta_archivo_seq2=clasificar(lista_archivos_traducidos2,nombre_ruta_archivo_unif_tabladom2,tipo_secuencia_escogida,base_referencia_escogida,nombre_ruta_secuencias_recortadas,tipo_secuencia_escogida_NuclProt,self.nombre_ruta_Ref_Gydb,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada)
                            clasificar_segundo_paso_genoma(nombre_ruta_archivo_gff2,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,base_referencia_escogida)
                            messagebox.showinfo(title="Finalización del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR HA FINALIZADO EXITÓSAMENTE."+"\n"+"\n"+str("Los archivos de salida se ubican en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK."))
       

                    #--Si base de referencia escogida es Gydb
                    if base_referencia_escogida==self.opcionesmenu_basereferencia_dominios[1]:
                        print("La base de referencia escogida es: "+str(self.opcionesmenu_basereferencia_dominios[1]))
                        comprobar_version_hmmer_bdRef(self.nombre_ruta_Ref_Gydb,nombre_programa='hmmscan') #LLamo a la funcion para comprobar que versiones del programa hmmscan y de base ref Rexdb tengan igual version
                        compresion_indexamiento_bd_ref(self.nombre_ruta_Ref_Gydb)
                        longitud_max_secuencia5=verificar_tamanio_secuencia(self.nombre_direccion_archivo)
                        if longitud_max_secuencia5<1e6 and tipo_secuencia_escogida=="Genomicas":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia5)+str(". Son secuencias de elementos transponibles y  no Genoma. Por favor escoja opcion TE"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de elementos transponibles y no de Genoma. Por favor escoja la opcion TE.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            #respuesta1=messagebox.askokcancel(title="Inicio del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR INICIARÁ."+"\n"+"\n"+str("El archivo de secuencias que se evaluará es:")+"\n"+str(nombre_archivo_entrada)+"\n"+"\n"+str("Los archivos de salida se ubicarán en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK para dar inicio y este proceso tardará unos minutos"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif longitud_max_secuencia5>1e6 and tipo_secuencia_escogida=="TE":
                            print("La longitud  de la secuencia mas larga es: "+str(longitud_max_secuencia5)+str(". Es un genoma y no secuencia de elementos transponibles. Por favor escoja opcion Genomicas"))
                            messagebox.showerror(title="ERROR DE TIPO DE SECUENCIA ESCOGIDA",message="EXISTE UN ERROR EN EL TIPO DE SECUENCIA ESCOGIDA."+"\n"+"\n"+str("Las secuencias son de un genoma y no de elementos transponibles. Por favor escoja la opcion Genomicas.")+"\n"+"\n"+str("De click en Ok para cancelar el proceso"))
                            print("PROCESO CANCELADO POR ERROR EN TIPO DE SECUENCIA ESCOGIDO")
                        elif (longitud_max_secuencia5>1e6 and tipo_secuencia_escogida=="Genomicas") or (longitud_max_secuencia5<1e6 and tipo_secuencia_escogida=="TE"):
                            palabra5=base_referencia_escogida
                            nombre_ruta_secuencias_recortadas='%s/secuencias_recortadas.%s.%s' % (ruta_nuevacarpeta_salida,palabra5,self.nombre_archivo_entrada)
                            #nombre_ruta_secuencias_recortads='{}/cu.{'.format(prueba44,'fasta')
                            print("Dentro se secuencias genomicas, esto es nombre ruta secuencias recortadas: "+str(nombre_ruta_secuencias_recortadas))
                            with open(nombre_ruta_secuencias_recortadas,'w') as identificacion_archivo_secuen_recortadas:
                                recortar_secuencias_genomicas(self.nombre_direccion_archivo,identificacion_archivo_secuen_recortadas)
                            #---Escanear con perfiles de modelos ocultos de markov
                            lista_archivos_traducidos2=escanear_hmm(nombre_ruta_secuencias_recortadas,ruta_nuevacarpeta_salida,self.nombre_ruta_Ref_Gydb,base_referencia_escogida,self.nombre_archivo_entrada,tipo_secuencia_escogida_NuclProt)
                            print("Despues de escanear hmm")
                            #----Proceso de clasificacion
                            palabra6=base_referencia_escogida
                            nombre_ruta_archivo_unif_tabladom2='%s/%s.%s.tabladom' % (ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,palabra6)
                            nombre_ruta_archivo_gff2,nombre_ruta_archivo_seq2=clasificar(lista_archivos_traducidos2,nombre_ruta_archivo_unif_tabladom2,tipo_secuencia_escogida,base_referencia_escogida,nombre_ruta_secuencias_recortadas,tipo_secuencia_escogida_NuclProt,self.nombre_ruta_Ref_Gydb,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada)
                            clasificar_segundo_paso_genoma(nombre_ruta_archivo_gff2,ruta_nuevacarpeta_salida,self.nombre_archivo_entrada,base_referencia_escogida)
                            messagebox.showinfo(title="Finalización del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR HA FINALIZADO EXITÓSAMENTE."+"\n"+"\n"+str("Los archivos de salida se ubican en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK."))
       

                #messagebox.showinfo(title="Finalización del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES LTR HA FINALIZADO EXITÓSAMENTE."+"\n"+"\n"+str("Los archivos de salida se ubican en la ruta: ")+"\n"+str(ruta_nuevacarpeta_salida)+"\n"+"\n"+str("De click en OK."))
            else:
                messagebox.showinfo(title="Cancelación del Proceso",message="EL PROCESO DE BÚSQUEDA DE RETROTRANSPOSONES SE CANCELÓ."+"\n"+"\n"+str("Para realizar una nueva búsqueda vuelva a escoger los campos en la ventana principal"))

        self.colocar_parametros_vacios()

#---Funcion para Ejecutar instrucciones como subproceso
def ejecutar(a):
    ejecutar_comando=subprocess.Popen(a,stdout=subprocess.PIPE,\
        stderr=subprocess.PIPE,shell=True)
    salida1=ejecutar_comando.communicate()
    #print("Dentro de ejecutar este es salida: "+str(salida1))
    estado1=ejecutar_comando.poll()
    #print("Dentro de ejecutar, este es estado: "+str(estado1))
    return salida1 + (estado1,)

##--Verificar programa hmmscan que se instalo con conda install hmmer y base de datos de referencia REXdb o Gydb tengan la misma version 3.3
def comprobar_version_hmmer_bdRef(nombre_ruta_bd_Ref,nombre_programa):
    primera_linea_bd_Ref=open(nombre_ruta_bd_Ref).readline()   #Metodo readline para leer la 1ra linea del archivo REXdb_protein..
    patron_version_bd_Ref=re.compile(r'HMMER\S+ \[([\w\.]+)')     #Es el patron que se va a buscar en la 1ra lineea del archivo REXdb
    result_patron_bd_Ref=patron_version_bd_Ref.search(primera_linea_bd_Ref)
    version_bd_Ref=result_patron_bd_Ref.groups()[0]
    version_reducida_bd_Ref=version_bd_Ref[:3]
    print("Esta es la version de base de referencia de dominios proteicos: "+str(version_reducida_bd_Ref)) #La version de base de referencia Rexdb en formato hmmer es 3.3 
    nombre_programa_hmmscan='which {}'.format(nombre_programa) #comprobar que exista programa hmmscan
    salida,error,estado=ejecutar(nombre_programa_hmmscan)
    print("Este es el estado de comprobacion de version:"+str(estado))
    if estado==0:  #si existe el programa hmmscan entonces busco su version, estado=0 indica que existe el programa hmmscan
        informacion_programa_hmmscan='{} -h'.format(nombre_programa)
        salida2,error2,estado2=ejecutar(informacion_programa_hmmscan)
        version_programa_hmmscan=re.compile(r'HMMER (\S+)').search(salida2.decode('utf-8')).groups()[0]
        version_reducida_programa_hmmscan=version_programa_hmmscan[:3]
        print("Esta es la version de programa hmmscan:"+str(version_reducida_programa_hmmscan)) #La version de hmmscan es 3.3
        if version_reducida_programa_hmmscan>=version_reducida_bd_Ref:  #Las versiones de hmmscan y de base de ref son correctas de 3.3
            print("Las versiones son correctas")
        elif version_reducida_programa_hmmscan<version_reducida_bd_Ref:
            print("La version del programa hmmscan es muy baja, por favor   vuelva a cargarlo con conda install")
        else:
            print("La version no es adecuada, asegurese que la base de referencia este en formato hmm")
    else:
        print("Las versiones no son correctas o el  programa hmmscan no esta instalado con conda install")

#------Compresion binaria e indexamiento del archivo de base de datos de referencia de dominios proteicos RExdb o Gydb
def compresion_indexamiento_bd_ref(nombre_ruta_bd_ref_2):  #Creacion de 4 archivos de hmmer, primero comprime a la base de referencia Rexdb o Gydb, etiqueta y crea 4 archivos hmmmr
    comando_hmmpress="hmmpress -f "+nombre_ruta_bd_ref_2      #Hmmpress permite crear 4 archivos de hmmer
    ejecuta_comando_hmmpress=subprocess.Popen(comando_hmmpress,shell=True,stdout=subprocess.PIPE)    #Subproceso para ejecutar comando hmmpress  en la terminal
    salida_indexam,error_indexam=ejecuta_comando_hmmpress.communicate()
    print("Esto es salida de indexamiento despues de  comando hmmpress:"+str(salida_indexam))
    print("Esto es error de indexamiento despues de comando hmmpress: "+str(error_indexam))

#---Comprobacion de tamano de secuencia
def verificar_tamanio_secuencia(nombre_ruta_Seq_In_1):
    #----Verificar que se trate de una secuencia de elementos transponibles y no de un genoma
    lista_longitud_secuencias=[]
    for secuencia1 in SeqIO.parse(open(nombre_ruta_Seq_In_1),'fasta'): #El metodo SeqIO.parse para leer secuencia y se crea como objeto SeqRecord iterable
        longitud_de_secuencia=len(secuencia1.seq)                  #Calculo longitud de secuencias, metodo seq selecciono solo las secuencias
        #print("Tipo de lista de longitud: "+str(type(longitud_de_secuencia)))
        lista_longitud_secuencias.append(longitud_de_secuencia)          #Creo una lista con la longitud de secuencias
    #print("Esta es la lista de longitud de secuencias: "+str(lista_longitud_secuencias))
    numero_secuencias=len(lista_longitud_secuencias)               #El metodo len determina el numero de elementos de la lista 
    #print("Este es el numero de secuencias a analizar: "+str(numero_secuencias))
    longitud_maxima_de_secuencia=max(lista_longitud_secuencias)           #Calculo la longitud maxima de la secuencia
    print("La longitud maxima es: "+str(longitud_maxima_de_secuencia))
    return longitud_maxima_de_secuencia

#-------Escanear secuencias con hmmscan de hmmer con base de referencia de dominios proteicos
def escanear_hmm(nombre_ruta_Seq_In_2,ruta_salida1,nombre_ruta_bd_ref_dominios,nombre_ref_dominios,nombre_Seq_In_2,tipo_secuencia_escogida_NuclProt2):
    #---Creando archivo de las secuencias a analizar
    nucleos=4
    #nombre_ruta_Seq_In_2=open(nombre_ruta_Seq_In_2)
    #Comprobacion de que secuencia es un tipo de base de datos
    comprobacion_de_secuencia=isinstance(nombre_ruta_Seq_In_2,IOBase)  #Compruebo que secuencia es una base de datos
    if not comprobacion_de_secuencia:  #Si el archivo no es de tipo IOBase entonces se abre el archivo y ya es de tipo IOBase base de datos
        nombre_ruta_Seq_In_2=open(nombre_ruta_Seq_In_2)
        print("Comprobando secuencia base de datos:"+str(comprobacion_de_secuencia))
    #---Creo y abro 8 archivos fasta, nucleos son 4 por 2 son 8 
    lista_nombres_archivos_fasta=[]  #Lista con nombres de los 8 archivos
    numero_archivos=nucleos*2
    for b in range(numero_archivos):        #Creo 8 archivos
        b+=1
        nombre_archivo_fasta='%s/archivo%s.fasta' % (ruta_salida1,b)
        lista_nombres_archivos_fasta+=[nombre_archivo_fasta]
        identificacion_archivo_fasta='archivo%s' % (b,)
        comando1='%s=open("%s","w")' % (identificacion_archivo_fasta,nombre_archivo_fasta)
        exec(comando1) #Ejecuta el comando
    #print("Lista de nombres de archivos fasta: "+str(lista_nombres_archivos_fasta))
    #---Se escribe en cada uno de los 8 archivos fasta las secuencias fasta del archivo de entrada
    y1=0
    objeto_iterable_todas_secuencias_fasta=identificar_agrupar_por_secuencia_fasta(nombre_ruta_Seq_In_2) #Guarda cada secuencia fasta recibida de la funcion idetificar-agrupar
    for secuencia2 in objeto_iterable_todas_secuencias_fasta:   #Itero cada secuencia del objeto iterable que tiene todas las secuencias fasta del archivo de entrada 
        numero_identificacion_archivo=(y1%8)+1 #Operacion modulo agrupar secuencias entre los 8 archivos creados, enviar cada secuencia a cada uno de los 8 creados anteriormente, rango de 0 a (N-1), de 0 a 7 por eso sumo 1 para que sea de 1 a 8 archivos   
        identificacion_archivo_unasecuencia_fasta='archivo%s' % (numero_identificacion_archivo)     #Creo la identificacion de 1 de los 8 archivos fasta en donde se escribira``
        comando2='%s.write(secuencia2)' % (identificacion_archivo_unasecuencia_fasta, ) #Escribo una secuencia dentro de uno de los 8 archivos fasta creados, guardo cada secuencia dentro de 1 de los 8 archivos fasta
        exec(comando2)  #Ejecuta el comando para ecribir una secuencia dentro de uno de los 8 archivos fasta 
        y1+=1
    #----Cerrar los 8 archivos fasta
    for y2 in range(numero_archivos):
        y2+=1
        comando3='archivo%s.close()' % (y2, )
        exec(comando3)
    #------Creacion de lista de archivos revisados llenos
    lista_nombres_archivos_revisados_fasta=[]  #Lista con nombres de los 8 archivos
    for nombre_archivo_fasta1 in lista_nombres_archivos_fasta:
        if os.path.getsize(nombre_archivo_fasta1)>0:   #Reviso archivos fasta que no estan vacios
            lista_nombres_archivos_revisados_fasta+=[nombre_archivo_fasta1] #Si archivo fasta no esta vacio creo lista de nombres de archivos fasta revisados
    #print("Esta lista de nombre de archivos revisados Fasta: "+str(lista_nombres_archivos_revisados_fasta))
    iterable_nombres_archivos_fasta=(nombre_archivo_fasta2 for nombre_archivo_fasta2 in lista_nombres_archivos_revisados_fasta) #creo tupla para iterar
    #print("Este es iterable nombres archivos fasta: "+str(iterable_nombres_archivos_fasta))
    #--Proceso de traduccion
    if tipo_secuencia_escogida_NuclProt2=='Nucleotidos': #Realiza la traduccion si la secuencia es nucleotidica
        lista_nombres_archivos_traducidos=list(multiprocesar(traducir,lista_nombres_archivos_revisados_fasta))  #lista de archivos traducidos: [1.fasta.aa,2.fasta.aa...], envio funcion traduccion e iterable de lista d nombres archivos fasta
        print("Dentro de escaner hmm, dentro de nucleotidos, Lista archivos de salida despues de traduccion: "+str(lista_nombres_archivos_traducidos))
    #----Si son secuencias de aminoaacidos  no se hace el proceso de traduccion
    if tipo_secuencia_escogida_NuclProt2=='Proteinas':
        lista_nombres_archivos_traducidos=lista_nombres_archivos_revisados_fasta
        print("En escanear hmm, en secuencia aminoacidos, esto es lista nombres archivos traducidos: "+str(lista_nombres_archivos_traducidos))

    #-----Creacion de lista de archivos de tabla de dominios proteicos, en base a los nombres de archivos traducidos a aminoacidos,  en los que se guardara el resultado de busqueda de secuencias contra base ref de dominios
    lista_nombres_archivos_tabla_dominios=[]  #Lista con nombres de archivos de tabla de dominios proteicos
    for nombre_archivo_traducido1 in lista_nombres_archivos_traducidos:
        lista_nombres_archivos_tabla_dominios+=[nombre_archivo_traducido1+'.tabladom']  #Aqui se almacenaran los nombres de los archivos de tabla de dominios proteicos
    print("Estos son nombres de archivos de tabla de dominios: "+str(lista_nombres_archivos_tabla_dominios))
    #-----Generacion de archivos de salida de tabla de dominios ,con el comando hmmscan para que compare secuencias contra la base de referencia con secuencias query 
    nombre_base_referencia_dominios=nombre_ruta_bd_ref_dominios
    lista_comandos_hmmscan=[
        'hmmscan --nobias --notextw --noali --domtblout {} {} {} > /dev/null'.format(
        nombre_archivo_tabla_dominios, nombre_base_referencia_dominios,nombre_archivo_traducido2)\
        for nombre_archivo_traducido2, nombre_archivo_tabla_dominios in zip(lista_nombres_archivos_traducidos,lista_nombres_archivos_tabla_dominios)]  #lista de archivos
    #print("esto es lista comandos hmmscan:"+str(lista_comandos_hmmscan))
    iterable_comandos=lista_comandos_hmmscan  #Lista con nombres de archivos de tabla de dominios proteicos
    #print("estoy dentro esto es iterble: "+str(iterable_comandos))
    #lista_resultados2=multiprocesar(ejecutar,iterable_comandos)
    try: [resultado2 for resultado2 in multiprocesar(ejecutar,iterable_comandos)] #Se itera porquela funcion multiprocesar devuelve una lista
    except:pass #Pass para pasar si hay errores
    #----Creo 1 unico archivo de salida de tablas de dominios que unificara los archivos antriormente creados y escribo en este archivo el contenido de los archivos de tabla de dominio
    palabra2=nombre_ref_dominios
    nombre_ruta_archivo_unif_tabladom='%s/%s.%s.tabladom' % (ruta_salida1,nombre_Seq_In_2,palabra2)
    #---Escribo en el archivo  unificado de salida de tablas de dominios, cada linea de los archivos de tabla de dominios copio y pongo en el archhivo unificado
    with open(nombre_ruta_archivo_unif_tabladom,'w') as identificacion_archivo_unif_tabladom:
        for archivo_tabla_dominio2 in lista_nombres_archivos_tabla_dominios:
            for linea2 in open(archivo_tabla_dominio2):
                identificacion_archivo_unif_tabladom.write(linea2)
    return lista_nombres_archivos_traducidos


#-----Resolver sobrelapes
def corregir_sobrelapes(lista_lineas):
    ultima_linea=None
    sobrelape_maximo=20
    lista_descartados=[]
    y3=0
    z3=0
    for linea4 in sorted(lista_lineas,key=lambda j:j[3]):
        descartado=None
        if ultima_linea:
            if linea4==ultima_linea:
                y3+=1
                pareja_lineas=[ultima_linea,linea4]
            else:
                sobrelape=max(0,min(linea4[4],ultima_linea[4])-max(linea4[3],ultima_linea[3]))
                sobrelape2=100*sobrelape/(min((linea4[4]-linea4[3]+1),(ultima_linea[4]-ultima_linea[3]+1)))
                if sobrelape2>sobrelape_maximo:
                    z3+=1
                    if linea4[5]>ultima_linea[5]:
                        pareja_lineas=[linea4,ultima_linea]
                    else:
                        pareja_lineas=[ultima_linea,linea4]
                else:
                    ultima_linea=linea4
                    continue
            mantenido,descartado=pareja_lineas
            lista_descartados+=[descartado]
        if not ultima_linea or descartado!=linea4:
            ultima_linea=linea4
    #print("Dentro de corregir sobrelapes, descartados:"+str(y3))
    print("Dentro de corregir sobrelapes, sobrelapes:"+str(z3))
    #print("Dentro de corregir sobrelapes, total:"+str(y3+z3))
    return sorted(set(lista_lineas)-set(lista_descartados),key=lambda k:k[3])

#-----Convierte variables a valores enteros o decimales
def convertir_valor_numerico(valor):
    try:return int(valor)
    except:
        try:return float(valor)
        except: return valor

#-----Obtengo una lista  de  los componentes de cada linea del archivo de tabla de dominio unificado
def listar_componentes_cadalinea_tabladominio(nombre_ruta_archivo_unif_tabladom4):
    print("En listar componentes cadalinea tabladominio, esto es archivo unificado salido tabla de dominios: "+str(nombre_ruta_archivo_unif_tabladom4))
    formato1='tabladom'
    #----Itero cada linea del archivo unificado de tabla de dominios
    for linea4 in open(nombre_ruta_archivo_unif_tabladom4):  #Abro el archivo paraProbarfasta.Rexdb.domtbl
        if linea4.startswith('#'):      #Si la linea inicia con # entonces salto y busco en la siguiente linea
            continue
        if formato1=='tabladom':      #si la linea no inicia con # corresponde a un elemento transponible identificado
            #print("En listar componentes cadalinea tabladominio, esto es linea completa:"+str(linea4))
            linea_sin_espacios=linea4.strip()     #Metodo strip remueve espacios al inicio y fin de la cadena que contiene info de la linea
            #----Defino nombres de compontente de la linea
            lista_componentes_linea=linea_sin_espacios.split()    # Metodo split separa una cadena en sus partes y le convierte en una lista, en la que cada palabra de la cadena es un item de la lista
            nombre_target,acesion_target,longitud_target,nombre_query,acesion_query,longitud_query, \
            valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia, \
            dominio1,dominio2,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio, \
            inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env,precision \
            =lista_componentes_linea[:22]   #Defino nombres de componentes de la linea del archivo unificado de tabla de dominios, que contiene los resultados de hmmer
            #print("esto es lita de componentes: "+str(lista_componentes_linea))
            #------Convertir valores en enteros
            longitud_target,longitud_query,dominio1,dominio2, \
            inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env= \
                list(map(convertir_valor_numerico, [longitud_target,longitud_query,dominio1,dominio2, \
                inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env])) #Map retorna un objeto iterable despues de aplicar una funcion a cada item iterable
            #-------Convertir valores en decimales
            valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio,precision= \
                list(map(convertir_valor_numerico, [valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio,precision]))
            #-------Se almacena el ultimo componente de la linea en una sola cadena y los separa con un espacio
            cadena_componentes_linea=' '.join(lista_componentes_linea[22:])  # Solo see almacena el ultimo componente que es descripcion  que tiene una -,almacena como cadena
            #print("estoy dentro de listar esto es cadena componentes linea:"+str(cadena_componentes_linea))
            #-----Calcular la cobertura 
            cobertura=round(1e2*(fin_hmm_query-inicio_hmm_query + 1) /  longitud_target, 1)
            yield lista_componentes_linea

#-----Analizar clado y dominio
def analizar_clado_dominio(nombre_target5,nombre_ref_dominios5):
    if nombre_ref_dominios5=='Gydb':
        lista_componentes_nombre_target=nombre_target5.split('_')  #Metodo split separa partes del string por ":" y convierte en una lista los items
        #print("estoy dentro de analizar_clado_dominio gydb")
        dominio=lista_componentes_nombre_target[0]
        #print("estoy dentro de analizar_clado_dominio gydb esto es dominio:"+str(dominio))
        clado='_'.join(lista_componentes_nombre_target[1:])
        #print("estoy dentro de analizar_clado_dominio gydb esto es clado:"+str(clado))
    elif nombre_ref_dominios5=='REXdb':
        #print("estoy dentro de analizar_clado_dominio rexdb") 
        lista_componentes_nombre_target=nombre_target5.split(':')  #Metodo split separa partes del string por ":" y convierte en una lista los items
        dominio=lista_componentes_nombre_target[1]
        #print("estoy dentro de analizar_clado_dominio rexdb esto es dominio:"+str(dominio))
        clado=lista_componentes_nombre_target[0].split('/')[-1]   
        #print("estoy dentro de analizar_clado_dominio Rexdb esto es clado:"+str(clado))
    else:
        #print("estody dentro de analizar_clado_dominioen else")
        clado,dominio=nombre_target5,nombre_target5
    return clado,dominio 

#----Analiza el mejor nombre de dominio y clado
def analizar_mejor_nombre_dominio_clado(mejor_nombre_tarjet,nombre_ref_dominios5):
    if nombre_ref_dominios5=='Gydb':
        #print("estoy dentro de analizar mejor nombre dominio clado dentro de gydb")
        lista_componentes_mejor_nombre_target=mejor_nombre_tarjet.split('_')
        dominio1,clado1=lista_componentes_mejor_nombre_target[0], '_'.join(lista_componentes_mejor_nombre_target[1:])
    elif nombre_ref_dominios5=='REXdb':
        #print("estoy dentro de analizar mejor nombre dominio clado en rexdb")
        dominio1=mejor_nombre_tarjet.split(':')[1]  
        #print("estoy dentro de analizar mejor nombre dominio clado esto es mejor nombre target :"+str(mejor_nombre_tarjet))
        #print("estoy dentro de analizar mejor nombre dominio clado esto es  dominio: "+str(dominio1))
        clado1=mejor_nombre_tarjet.split(':')[0].split('/')[-1]  
        #print("estoy dentro de analizar mejor nombre dominio esto es clado:"+str(clado1))
    else:
        #print("estody dentro de analizar mejor nombre dominio clado en else")
        dominio1,clado1=mejor_nombre_tarjet,mejor_nombre_tarjet
    return dominio1,clado1    

#-----Obtiene el mejor acierto del analisis de dominios proteicos con modelos HMM, en base a la comparacion de puntajes de los resultados obtenidos del analisis de dominios proteicos contra base de ref dominios  hmmer
def escoger_mejor_puntaje(nombre_ruta_archivo_unif_tabladom4,nombre_ref_dominios4,tipo_secuencia4,tipo_secuencia_escogida_NuclProt4):
    print("esto es archivo unificado salido tabla de dominios, dentro de escoger mejor puntaje: "+str(nombre_ruta_archivo_unif_tabladom4))
    diccionario_mejor_puntaje_mMarkov={} #Creo un diccionario 
    iterable_listas_componentes_linea=listar_componentes_cadalinea_tabladominio(nombre_ruta_archivo_unif_tabladom4)
    #----Itero cada lista que contiene los #componentes de cada linea del archivo unificado de tabla de dominios
    for lista_componentes_linea2 in iterable_listas_componentes_linea:
        #----Defino nombres de compontente de la linea
        nombre_target,acesion_target,longitud_target,nombre_query,acesion_query,longitud_query, \
        valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia, \
        dominio1,dominio2,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio, \
        inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env,precision \
        =lista_componentes_linea2[:22]   #Defino nombres de componentes de la linea del archivo unificado de tabla de dominios, que contiene los resultados de hmmer
        #print("esto es lita de componentes: "+str(lista_componentes_linea2))
        #-----Obtener el nombre sufijo del nombre del query
        palabra_sufijo2_query=nombre_query.split('|')[-1]  #Obtengo la palabra posterior del nombre del query, puede ser aa o rev_aa
        #print("esto es palabra sufijo: "+str(palabra_sufijo2_query))
        #------Obtener identificacion del query
        if tipo_secuencia_escogida_NuclProt4=='Nucleotidos' and (palabra_sufijo2_query.startswith('aa') or palabra_sufijo2_query.startswith('rev_aa')):
            #print("Dentro de escoger mejor puntaje y dentro de if secuencia nucleotidos")
            palabra_identificacin_query=nombre_query.split('|')[:-1]
            identificacion_query='|'.join(palabra_identificacin_query) #Obtengo la identificacion del query puede ser: Osat_TEdenovoGr-B-G1014-Map3
        else:
            identificacion_query=nombre_query
            #print("Dentro de escoger mejor puntaje y dentro de else")
        #print("Dentro de escoger mejor puntaje, esto es identificacion del query: "+str(identificacion_query))
        #----Obtener clado y dominio de la linea
        clado,dominio=analizar_clado_dominio(nombre_target,nombre_ref_dominios4)
        #---Creacion de diccionario
        clave=(identificacion_query,)
        if tipo_secuencia4=="Genomicas":
            clave+=(inicio_env,fin_env)
        #-----Normalizar puntaje
        #print("Dentro de escoger mejor puntaje esto es anterior puntaje secuencia: "+str(puntaje_secuencia))
        #print("estoy dentro de  escogr mejor puntaje esto es valor_e secuencia antes:"+str(valor_e_secuencia))
        puntaje_secuencia=round(float(puntaje_dominio)/ int(longitud_target),2)  
        valor_e_secuencia=i_valor_e_dominio  
        #print("Dentro de escoger mejor puntaje esto es puntaje dominio; "+str(puntaje_dominio))
        #print("dentro de escoger mejor puntaje esto  es  longitud target:"+str(longitud_target))
        #print("Dentro de escoger mejor puntaje esto es nuevo puntaje secuenia: "+str(puntaje_secuencia))
        #print("estoy dentro de escoger mejor putnaje esto  es nuevo valor_e_secuencia:"+str(valor_e_secuencia))
        lista_componentes_linea2=[nombre_target,acesion_target,longitud_target,nombre_query,acesion_query,longitud_query, \
        valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia, \
        dominio1,dominio2,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio, \
        inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env,precision \
        ] 
        #print("estoy dentro de escoger mejor puntaje esto  es nueva lista de componentes linea:"+str(lista_componentes_linea2))
        #-----
        #print("DENTRO DE ESCOGER MEJOR PUNTAJE ESTO ES DOMINIO: "+str(dominio))
        if nombre_ref_dominios4=='REXdb':
            dominio3=dominio.split('-')[1]
            #print("estoy dentro de escoger mejor puntaje esto es dominio3: "+str(dominio3))
            if dominio3=='aRH' and tipo_secuencia4=="TE":
                dominio3='RH'
                #print("DENTRO DE ESCOGER MEJOR PUNTAJE ENTRE por Tpase: esto es dominio3:  "+str(dominio3))
            if dominio3=='TPase' and tipo_secuencia4=="TE":
                dominio3='INT'
                #print("DENTRO DE ESCOGER MEJOR PUNTAJE ENTRE por arh: esto es dominio3:  "+str(dominio3))
            clave+=(dominio3,)
            if clave in diccionario_mejor_puntaje_mMarkov:
                lista_componentes_linea_mejor1=diccionario_mejor_puntaje_mMarkov[clave]
                #----Defino nombres de compontentes de la linea mejor
                mejor_nombre_target,mejor_acesion_target,mejor_longitud_target,mejor_nombre_query,mejor_acesion_query,mejor_longitud_query, \
                mejor_valor_e_secuencia,mejor_puntaje_secuencia,mejor_sesgo_correccion_secuencia, \
                mejor_dominio1,mejor_dominio2,mejor_c_valor_e_dominio,mejor_i_valor_e_dominio,mejor_puntaje_dominio,mejor_sesgo_correcion_dominio, \
                mejor_inicio_hmm_query,mejor_fin_hmm_query,mejor_inicio_alineamiento_target,mejor_fin_alineamineto_target,mejor_inicio_env,mejor_fin_env,mejor_precision \
                =lista_componentes_linea_mejor1[:22]   #Defino nombres de componentes de la mejor puntuacion de la linea del archivo unificado de tabla de dominios, que contiene los resultados de hmmer
                #print("esto es lista compoentes mejor1: "+str(lista_componentes_linea_mejor1))
                #print("esto es mejor puntaje secuencia: "+str(mejor_puntaje_secuencia))
                #print("esto es puntaje de seuencia: "+str(puntaje_secuencia))
                #print("esto es mejor nombre del tarjet "+str(mejor_nombre_target))
                #print("esto es valor e secuencia: "+str(valor_e_secuencia))
                #print("esto es mejor valor e secuencia: "+str(mejor_valor_e_secuencia))
                #----Normalizacion de mejor secuencia del puntaje de secuencia y del valor-e de secuencia
                mejor_puntaje_dominio=float(mejor_puntaje_dominio)
                mejor_longitud_target=int(mejor_longitud_target)
                mejor_puntaje_secuencia=round(mejor_puntaje_dominio/ mejor_longitud_target,2)  
                mejor_i_valor_e_dominio=float(mejor_i_valor_e_dominio)
                mejor_valor_e_secuencia=mejor_i_valor_e_dominio
                #print("eesto es mejor puntaje dominio: "+str(mejor_puntaje_dominio))
                #print("esto es mejor longitud del target: "+str(mejor_longitud_target))
                #print("esto es nuevo mejor puntaje de secuencia: "+str(mejor_puntaje_secuencia))
                #print("esto es mejor_i_valor_e_dominio: "+str(mejor_i_valor_e_dominio))
                #print("esto es nuevo mejor valor e secuencia: "+str(mejor_valor_e_secuencia))
                lista_componentes_linea_mejor1=[mejor_nombre_target,mejor_acesion_target,mejor_longitud_target,mejor_nombre_query,mejor_acesion_query,mejor_longitud_query, \
                mejor_valor_e_secuencia,mejor_puntaje_secuencia,mejor_sesgo_correccion_secuencia, \
                mejor_dominio1,mejor_dominio2,mejor_c_valor_e_dominio,mejor_i_valor_e_dominio,mejor_puntaje_dominio,mejor_sesgo_correcion_dominio, \
                mejor_inicio_hmm_query,mejor_fin_hmm_query,mejor_inicio_alineamiento_target,mejor_fin_alineamineto_target,mejor_inicio_env,mejor_fin_env,mejor_precision \
                ] 
                #print("estoy dentro de escoger mejor puntaje esto  es nueva lista de componentes mejor linea:"+str(lista_componentes_linea_mejor1))
                if puntaje_secuencia >mejor_puntaje_secuencia: 
                    mejor_dominio, _=analizar_mejor_nombre_dominio_clado(mejor_nombre_target,nombre_ref_dominios4)
                    if dominio==mejor_dominio:
                        diccionario_mejor_puntaje_mMarkov[clave]=lista_componentes_linea2
                    elif inicio_env <= mejor_fin_env and fin_env>=mejor_inicio_env:
                        diccionario_mejor_puntaje_mMarkov[clave]=lista_componentes_linea2
            else:
                diccionario_mejor_puntaje_mMarkov[clave]=lista_componentes_linea2
        else: #si es gydb
            #print("estoy dentro de escoger mejor puntaje, entre  por else gydb")
            clave+=(dominio,)
            #------Crear diccionario con mejores puntajes
            if clave in diccionario_mejor_puntaje_mMarkov: #Si la clave se encuentra en el diccionario le analizo el valor de puntaje, la clave es el nombre del query con el dominio como Ty-Gag
                lista_componentes_linea_mejor1=diccionario_mejor_puntaje_mMarkov[clave]
                #----Defino nombres de compontentes de la linea mejor
                mejor_nombre_target,mejor_acesion_target,mejor_longitud_target,mejor_nombre_query,mejor_acesion_query,mejor_longitud_query, \
                mejor_valor_e_secuencia,mejor_puntaje_secuencia,mejor_sesgo_correccion_secuencia, \
                mejor_dominio1,mejor_dominio2,mejor_c_valor_e_dominio,mejor_i_valor_e_dominio,mejor_puntaje_dominio,mejor_sesgo_correcion_dominio, \
                mejor_inicio_hmm_query,mejor_fin_hmm_query,mejor_inicio_alineamiento_target,mejor_fin_alineamineto_target,mejor_inicio_env,mejor_fin_env,mejor_precision \
                =lista_componentes_linea_mejor1[:22]   #Defino nombres de componentes de la mejor puntuacion de la linea del archivo unificado de tabla de dominios, que contiene los resultados de hmmer
                #print("esto es lista compoentes mejor1: "+str(lista_componentes_linea_mejor1))
                #print("esto es mejor puntaje secuencia: "+str(mejor_puntaje_secuencia))
                #print("esto es puntaje de seuencia: "+str(puntaje_secuencia))
                #print("esto es mejor nombre del tarjet "+str(mejor_nombre_target))
                #print("esto es valor e secuencia: "+str(valor_e_secuencia))
                #print("esto es mejor valor e secuencia: "+str(mejor_valor_e_secuencia))
                #----Normalizacion de mejor secuencia del puntaje de secuencia y del valor-e de secuencia
                mejor_puntaje_dominio=float(mejor_puntaje_dominio)
                mejor_longitud_target=int(mejor_longitud_target)
                mejor_puntaje_secuencia=round(mejor_puntaje_dominio/ mejor_longitud_target,2)  
                mejor_i_valor_e_dominio=float(mejor_i_valor_e_dominio)
                mejor_valor_e_secuencia=mejor_i_valor_e_dominio
                #print("eesto es mejor puntaje dominio: "+str(mejor_puntaje_dominio))
                #print("esto es mejor longitud del target: "+str(mejor_longitud_target))
                #print("esto es nuevo mejor puntaje de secuencia: "+str(mejor_puntaje_secuencia))
                #print("esto es mejor_i_valor_e_dominio: "+str(mejor_i_valor_e_dominio))
                #print("esto es nuevo mejor valor e secuencia: "+str(mejor_valor_e_secuencia))
                #lista_componentes_linea_mejor1=[nombre_target,acesion_target,longitud_target,nombre_query,acesion_query,longitud_query, \
                #valor_e_secuencia,puntaje_secuencia,sesgo_correccion_secuencia, \
                #dominio1,dominio2,c_valor_e_dominio,i_valor_e_dominio,puntaje_dominio,sesgo_correcion_dominio, \
                #inicio_hmm_query,fin_hmm_query,inicio_alineamiento_target,fin_alineamineto_target,inicio_env,fin_env,precision \
                #] 
                #print("estoy dentro de escoger mejor puntaje esto  es nueva lista de componentes mejor linea:"+str(lista_componentes_linea_mejor1)) 
                if puntaje_secuencia >mejor_puntaje_secuencia:   #Si la nueva secuencia tiene un puntaje mayor al anterior entonces le agrego al diccionario, si es menor no le agrego
                    diccionario_mejor_puntaje_mMarkov[clave]=lista_componentes_linea2
            else:  #Si la clave no se encuentra en el diccionario, le agrego al diccionario
                diccionario_mejor_puntaje_mMarkov[clave]=lista_componentes_linea2
    return diccionario_mejor_puntaje_mMarkov

#----Determinar orden y superamilia con Rexdb
def deter_orden_superfamilia_rexdb(clave_clado_mayor_frecuencia2):  
    if clave_clado_mayor_frecuencia2.startswith('Class_I/LTR/Ty1_copia'):
        orden,superfamilia='LTR','Copia'
        #print("Dentro de determinar orden superfamilia Rexdb, en if Class_I copia")
    elif clave_clado_mayor_frecuencia2.startswith('Class_I/LTR/Ty3_gypsy'):
        orden,superfamilia='LTR','Gypsy'
        #print("Dentro de determinar orden superfamilia Rexdb, en elif Class_I gypsy")
    elif clave_clado_mayor_frecuencia2.startswith('Class_I/LTR/'):
        orden,superfamilia=clave_clado_mayor_frecuencia2.split('/')[1:3]
        #print("Dentro de determinar orden superfamilia Rexdb, en elif Class_I LTR")
    elif clave_clado_mayor_frecuencia2.startswith('Class_I/'):
        #print("Dentro de determinar orden superfamilia Rexdb, en elif Class_I")
        try: orden,superfamilia=clave_clado_mayor_frecuencia2.split('/')[1:3]
        except ValueError: orden,superfamilia=clave_clado_mayor_frecuencia2.split('/')[1],'Desconocido'
    elif clave_clado_mayor_frecuencia2.startswith('Class_II/'):
        #print("Dentro de determinar orden superfamilia Rexdb, en elif Class_II")
        try: orden,superfamilia=clave_clado_mayor_frecuencia2.split('/')[2:4]
        except ValueError: orden,superfamilia=clave_clado_mayor_frecuencia2.split('/')[2],'Desconocido'
    elif clave_clado_mayor_frecuencia2.startswith('NA'):
        #print("Dentro de determinar orden superfamilia Rexdb, en elif NA")
        orden,superfamilia='LTR','Retrovirus'
    else:
        print("Clado desconocido")          
    return orden,superfamilia

#-----Determina el orden, superfamilia y clado
def determinar_orden_superfamilia_clado_Rexdb(lista_dominios6,lista_clados6):
    estructura_correcta={
        ('LTR','Copia'):['GAG','PROT','INT','RT','RH'],
        ('LTR','Gypsy'):['GAG','PROT','RT','RH','INT'],
        ('LTR','Bel-Pao'):['GAG','PROT','RT','RH','INT'],
        }
    diccionario_frecuencia_clado=Counter(lista_clados6)
    #print("Dentro de determinar orden superfamilia Clado Rexdb esto es lista de genes o dominios: "+str(lista_dominios6))
    #print("Dentro de determinar orden superfamilia Clado Rexdb esto es lista de clados: "+str(lista_dominios6))
    #print("Dentro de determinar orden superfamilia Clado Rexdb, esto es diccionario de frecuencia de clado: "+str(diccionario_frecuencia_clado))
    clado_mayor_frecuencia=max(diccionario_frecuencia_clado,key=lambda i:diccionario_frecuencia_clado[i])
    #---Determinar orden y superfamilia
    orden,superfamilia=deter_orden_superfamilia_rexdb(clado_mayor_frecuencia)
    #print("Dentro determinar orden superfamilia  Clado  RExcb, esto es clave de clado de mayor frecuencia: "+str(clado_mayor_frecuencia))
    #print("Dentro de determinar orden superfamilia Clado  Rexdb  esto es orden: "+str(orden))
    #print("Dentro de determinar orden superfamilia  Clado Rexdb eesto es superfamilia: "+str(superfamilia))
    if len(diccionario_frecuencia_clado)==1 or diccionario_frecuencia_clado[clado_mayor_frecuencia] >1:
        clado_mayor_frecuencia=clado_mayor_frecuencia.split('/')[-1]
        #print("Dentro de determinar orden superfamilia Clado RExdb esto es clado de mayor frecuencia: "+str(clado_mayor_frecuencia))
    elif len(diccionario_frecuencia_clado)>1:
        clado_mayor_frecuencia='mezcla'
        lista_superfamlias=[deter_orden_superfamilia_rexdb(clado2)[1] for clado2 in lista_clados6]
        if len(Counter(lista_superfamlias))>1:
            superfamilia='mezcla'
            lista_ordenes=[deter_orden_superfamilia_rexdb(clado3)[0] for clado3 in lista_clados6]
            if len(Counter(lista_ordenes))>1:
                orden='mezcla'
    try:
        genes_ordenados=estructura_correcta[(orden,superfamilia)]
        #print("Dentro de determinar orden superfamilia Clado Rexdb, dentro de try esto es genes ordenados: "+str(genes_ordenados))
        genes_revisar=[gen3 for gen3 in lista_dominios6 if gen3 in set(genes_ordenados)]
        #print("Dentro de determinar orden superfamilia Clado Rexdb, dentro de try esto es genes revisar: "+str(genes_revisar))
        if genes_ordenados==genes_revisar:
            codificacion='Si'
            #print("Dentro de determinar orden superfamilia Clado Rexdb, en try, en if, genes ordenados y codificacion si")
        else:
            codificacion='No'
            #print("Dentro de determinar orden superfamilia Clado Rexdb, en try, en else coddicacion no")
    except KeyError:
        codificacion='Desconocida'
    #print("dentro de determinar orden superfamilia Clado Rexdb, esto es clado de mayor frecuencia: "+str(clado_mayor_frecuencia))
    if superfamilia not in {'Copia','Gypsy'}:
        clado_mayor_frecuencia='Desconocido'
        #print("Dentro de determinar orden superfamilia Clado Rexdb, en if superamilia no esta en copia gypsy")
        #print("Dentro de determinar orden superfamilia Clado Rexdb, en if superamilia no esta en copia gypsy, esto es clado mayor frecuencia: "+str(clado_mayor_frecuencia))
    if clado_mayor_frecuencia.startswith('Ty'):
        clado_mayor_frecuencia='Desconocido'
        #print("Dentro de determinar orden superfamilia Clado Rexdb,  en   if clado mayour frecuencia inicia con start Ty")
        #print("Dentro de determinar orden superfamilia  Clado Rexdb, en if clado mayor frecuencia inicia con starta Ty, esto es clado mayor frecuencia: "+str(clado_mayor_frecuencia))
    return orden,superfamilia,clado_mayor_frecuencia,codificacion

#----Determinar orden, clado, superfamilias,  de inforamcion de Gydb
def determinar_referencia_orden_superfamilia_clado_Gydb(nombre_ruta_Ref_Gydb7):
    nombre_ruta_Ref_Gydb_infor='%s.information' % (nombre_ruta_Ref_Gydb7)
    #print("En determinar  referencia orden superfamilia clado Gydb, esto es nombre ruta gydb informacion: "+str(nombre_ruta_Ref_Gydb_infor))
    diccionario_clado={
        'Ty_(Pseudovirus)':'pseudovirus',
        'Cer2-3':'cer2-3',
        '412/Mdg1':'412_mdg1',
        'TF1-2':'TF',
        'Micropia/Mdg3':'micropia_mdg3',
        'CoDi-I':'codi_I',
        'CoDi-II':'codi_II',
        '17.6':'17_6',
        }
    y7=0
    for linea7 in open(nombre_ruta_Ref_Gydb_infor):
        y7+=1
        lista_elementos_linea7=linea7.strip().split('\t')
        if y7==1:
            lista_enunciados_claves_linea7=lista_elementos_linea7
            continue
        lista_tupla_enunciados_elementos_linea7=list(zip(lista_enunciados_claves_linea7,lista_elementos_linea7))  #----Tupla
        diccionario_enunciados_elementos_linea7=dict(lista_tupla_enunciados_elementos_linea7)
        #print("En determinar referencia  orden superfamilia clado Gydb, esto ees diccionario enunciados elementos linea7: "+str(diccionario_enunciados_elementos_linea7))
        if diccionario_enunciados_elementos_linea7['Clade']=='NA':
            clado7=diccionario_enunciados_elementos_linea7['Cluster_or_genus']
        else:
            clado7=diccionario_enunciados_elementos_linea7['Clade']
        familia7=diccionario_enunciados_elementos_linea7['Family']
        superfamilia7=familia7.split('/')[-1]
        if superfamilia7=='Retroviridae':
            clado7=diccionario_enunciados_elementos_linea7['Cluster_or_genus'].replace('virus','viridae')
        if superfamilia7=='Retrovirus':
            superfamilia7='Retroviridae'
        if diccionario_enunciados_elementos_linea7['System'] in {'LTR_retroelements','LTR_Retroelements','LTR_retroid_elements'}:
            orden7='LTR'
        else:  
            orden7=diccionario_enunciados_elementos_linea7['System']
        yield clado7,orden7,superfamilia7
        if clado7 in diccionario_clado:
            clado7=diccionario_clado[clado7]
            yield clado7,orden7,superfamilia7
            if clado7=='412_mdg1':
                clado7='412-mdg1'
                yield clado7,orden7,superfamilia7

        clado7=clado7.replace('-','_')
        yield clado7,orden7,superfamilia7
        clado7=clado7.lower()
        yield clado7,orden7,superfamilia7
    orden7,superfamilia7,clado7,diccionario_enunciados_elementos_linea7=['LTR','Copia','ty1/copia',{}]
    yield clado7,orden7,superfamilia7
    diccionario_otros_clados={
        'retroelement':'LTR',
        'retroviridae':'LTR',
        'B-type_betaretroviridae':'LTR',
        'D-type_betaretroviridae':'LTR',
        'caulimoviruses':'LTR',
        'caulimoviridae_dom2':'LTR',
        'errantiviridae':'LTR',
        'retropepsins':'LTR',
        'VPX_retroviridae':'LTR',
        'cog5550':'Desconocido',
        'ddi':'Desconocido',
        'dtg_ilg_template':'Desconocido',
        'saspase':'Desconocido',
        'GIN1':'Desconocido',
        'shadow':'Desconocido',
        'all':'Desconocido',
        'pepsins_A1a':'Desconocido',
        'pepsins_A1b':'Desconocido',
        }
    for clado8,orden8 in list(diccionario_otros_clados.items()):
        orden7,superfamilia7,clado7,diccionario_enunciados_elementos_linea7=[orden8,'Desconocido',clado8,{}]
        yield clado7,orden7,superfamilia7

#------Determinar orden, superfamilia, clado para Gydb
def determinar_orden_superfamilia_clado_Gydb(lista_dominios6,lista_clados6,nombre_ruta_Ref_Gydb6):
    estructura_correcta={
        ('LTR','Copia'):['GAG','AP','INT','RT','RNaseH'],
        ('LTR','Gypsy'):['GAG','AP','RT','RNaseH','INT'],
        ('LTR','Pao'):['GAG','AP','RT','RNaseH','INT'],
        ('LTR','Retroviridae'):['GAG','AP','RT','RNaseH','INT','ENV'],
        ('LTR','Caulimoviridae'):['GAG','AP','RT','RNaseH',],
        }
    #----Obtener informacion de clado
    #Crear diccionario de referencia de Gydb 
    diccionario_clado={}
    lista_valores=determinar_referencia_orden_superfamilia_clado_Gydb(nombre_ruta_Ref_Gydb6)
    for elemento9 in lista_valores:
        #print("En determinar orden superfamiilia clado Gydb esto es elemento: "+str(elemento9))
        clado9=elemento9[0]
        orden9=elemento9[1]
        superfamilia9=elemento9[2]
        #diccionario_clado={clado9:(orden9,superfamilia9)}
        diccionario_clado[clado9]=(orden9,superfamilia9)
   #print("En determinar orden superfamilia clado Gydb esto es diccionario clado: "+str(diccionario_clado))

    diccionario_frecuencia_clado=Counter(lista_clados6)

    #print("Dentro de determinar orden superfamilia cladoGydb, esto es diccionario frecuencia clado: "+str(diccionario_frecuencia_clado))
    clave_clado_mayor_frecuencia=max(diccionario_frecuencia_clado,key=lambda i:diccionario_frecuencia_clado[i])
    #print("Dentro de determinar orden superfamilia cladoGydb, esto es clave clado mayor frecuencia: "+str(clave_clado_mayor_frecuencia))
    try: (orden,superfamilia)=diccionario_clado[clave_clado_mayor_frecuencia]
    except KeyError:
        (orden,superfamilia)=('Desconocido','Desconocido')
        #print("Dentro de determinar orden superfamiia clado Gydb, clado desconocido")
    if len(diccionario_frecuencia_clado)==1 or diccionario_frecuencia_clado[clave_clado_mayor_frecuencia]>1:
        clave_clado_mayor_frecuencia=clave_clado_mayor_frecuencia
    elif len(diccionario_frecuencia_clado)>1:
        clave_clado_mayor_frecuencia='mezcla'
        lista_superfamilias=[diccionario_clado.get(clado3,[None,None])[1] for clado3 in lista_clados6]
        if len(Counter(lista_superfamilias))>1:
            superfamilia='mezcla'
            lista_ordenes=[diccionario_clado.get(clado4,[None,None])[0] for clado4 in lista_clados6]
            if len(Counter(lista_ordenes))>1:
                orden='mezcla'
    try:
        genes_ordenados=estructura_correcta[(orden,superfamilia)]
        genes_revisar=[gen3 for gen3 in lista_dominios6 if gen3 in set(genes_ordenados)]
        if genes_ordenados==genes_revisar:
            codificacion='Si' #Estructura completa de genes y en el mismo orden
        else:
            codificacion='No'
    except KeyError:
        codificacion='Desconocido'
    return orden,superfamilia,clave_clado_mayor_frecuencia,codificacion

#-----Realiza clasificacion de elementos transponibles
def clasificar(lista_nombres_archivos_traducidos3,nombre_ruta_archivo_unif_tabladom3,tipo_secuencia3,nombre_ref_dominios3,nombre_direccion_archivo_entrada3,tipo_secuencia_escogida_NuclProt3,nombre_ruta_Ref_Gydb3,ruta_nuevacarpeta_salida3,nombre_archivo_entrada3):
    cobertura_minima=20
    maximo_valor_e=1e-3
    minima_precision=0.6
    lista_lineas_gff=[]
    diccionario_mejor_puntaje=escoger_mejor_puntaje(nombre_ruta_archivo_unif_tabladom3,nombre_ref_dominios3,tipo_secuencia3,tipo_secuencia_escogida_NuclProt3) #obtiene el mejor puntaje del resultado del analisis de dominios proteicos contra base de ref de dominios
    #print("dentro de clasificar esto es diccionario mejor puntaje: "+str(diccionario_mejor_puntaje))
    diccionario_componentes_secuencias_traducidas={}
    for nombre_archivo_traducido in lista_nombres_archivos_traducidos3:
        for componente_secuencia_traducida in SeqIO.parse(nombre_archivo_traducido,'fasta'):
            diccionario_componentes_secuencias_traducidas[componente_secuencia_traducida.id]=componente_secuencia_traducida
    #----
    for clave2,valor2 in list(diccionario_mejor_puntaje.items()): #Metodo items obtiene clave y valor del diccionario y le pongometodo list para obtener en una lista el resultado
        if tipo_secuencia3=="Genomicas":
            #print("Dentro de clasificar dentro de secuencia genomica")
            identificacion_query2,inicio_env2,fin_env2,dominio5=clave2
        else:
            identificacion_query2,dominio5=clave2
         #----Defino nombres de compontente de la linea
        nombre_target2,acesion_target2,longitud_target2,nombre_query2,acesion_query2,longitud_query2, \
        valor_e_secuencia2,puntaje_secuencia2,sesgo_correccion_secuencia2, \
        dominio1_2,dominio2_2,c_valor_e_dominio2,i_valor_e_dominio2,puntaje_dominio2,sesgo_correcion_dominio2, \
        inicio_hmm_query2,fin_hmm_query2,inicio_alineamiento_target2,fin_alineamineto_target2,inicio_env2,fin_env2,precision2 \
        =valor2[:22]   #Defino nombres de componentes de la linea del archivo unificado de tabla de dominios, que contiene los resultados de hmmer
        #print("Dentro de clasificar esto es lita de componentes de la linea: "+str(valor2))
        #-----Calcular la cobertura 
        cobertura2=round(1e2*(int(fin_hmm_query2)-int(inicio_hmm_query2) + 1) /  int(longitud_target2), 1)
        #print("Dentro de clasificar esto es cobertura;  "+str(cobertura2))
        if cobertura2<cobertura_minima or float(valor_e_secuencia2)>maximo_valor_e or float(precision2)<minima_precision:
            #print("LA cobertura, evalue y probabilidad es muy baja")
            continue
        solo_identificacion_query=identificacion_query2
        gen5,clado5=analizar_mejor_nombre_dominio_clado(nombre_target2,nombre_ref_dominios3)
        #print("Estoy dentro de clasificar despues de analizar  dominio y clado esto es gen: "+str(gen5))
        #print("Estoy dentro de clasificar despues de analizar  dominio y clado esto es clado: "+str(clado5))
        #-----Obtener el valor del dominio
        if nombre_ref_dominios3=='REXdb':
            dominio5=gen5.split('-')[1]
        #---Colocar el nombre de identificacion para archivo gff 
        #print("Dentro de clasificar esto es identificacion del query: "+str(identificacion_query2))
        patron2=re.compile(r'[;=\|]')  #Se busca el patron
        resultado_verif_patron2=patron2.sub("_",identificacion_query2)  #Se busca el patron en el nombre de identificacion del query
        #print("Dentro de clasificar esto es patron2: "+str(patron2))
        #print("Dento de clasificar esto es resultado_verif_patron2: "+str(resultado_verif_patron2))
        identificacion_queryTarget_gff='{}|{}'.format(resultado_verif_patron2,nombre_target2)  #Coloco la identificacion para el archivo 
        #print("Dentro de clasificar esto es identificacion del archivo gff: "+str(identificacion_queryTarget_gff))
        #-----Obtener secuencia de aminoacidos del archivo fasta traducido 
        #valor2=diccionario_componentes_secuencias_traducidas[nombre_query2]
        #secuencia2=diccionario_componentes_secuencias_traducidas[nombre_query2].seq[inicio_env2-1:fin_env2]
        try: secuencia_archivo_gff=diccionario_componentes_secuencias_traducidas[nombre_query2].seq[int(inicio_env2)-1:int(fin_env2)]  #Obtengo la secuencia fasta.aa en esas posiciones
        except KeyError as e:
            raise KeyError('{}\nEl archivo dombtl Hmm no es consistente con el archivo de secuencia')
        secuencia_archivo_gff=str(secuencia_archivo_gff)  #Se convierte la cadena en formato string
        #print("Dentro de clasificar esto es secuencia en formato string: "+str(secuencia_archivo_gff))
        #-----Analisis de la hebra y el marco
        if tipo_secuencia_escogida_NuclProt3=='Nucleotidos':  #Si la secuencia de entrada es nucleotidica se realiza traduccion para obtener aminoacidos, pero si las secuencias de entrada ya son aminoacidos o proteinas no se realiza la traduccion
            #---Crear una diccionario con identificacion de secuencia y longitud de secuencia
            #print("Dentro de clasificar dentro de if tipo de secuencia nucleotidos")
            diccionario_id_longitud_secuencia3={}
            for grupo_secuencias3 in SeqIO.parse(open(nombre_direccion_archivo_entrada3),'fasta'): #El metodo SeqIO.parse para leer secuencia y se crea como objeto SeqRecord iterable
                id_secuencia3=grupo_secuencias3.id                  #Obtengo identificacion de secuencias, metodo id selecciono solo la identificacion
                #print("La id de la secuencia es: "+str(id_secuencia3))
                longitud_secuencia3=len(grupo_secuencias3.seq)                  #Calculo longitud de secuencias, metodo seq selecciono solo las secuencias
                #print("Longitud se secuencia es: "+str(longitud_secuencia3))
                diccionario_id_longitud_secuencia3[id_secuencia3]=longitud_secuencia3   
            #print("Dentro de clasificar el diccionario de id y longitud de secuencias es"+str(diccionario_id_longitud_secuencia3))  #Creo un diccionario con claves la identificacion, y valor la longitud de secuencia 
            terminacion_nombre_query=nombre_query2.split('|')[-1] #Obtengo la terminacion del nombre del query puede ser rev de reveersa o aa
            #print("Dentro de clasificar esto es  terminacion nombre query: "+str(terminacion_nombre_query))
            if terminacion_nombre_query.startswith('rev'):
                hebra='-'
                #print("Dentro de clasificar dentro de nucleotidos en if reversa es hebra negativa, esto es  hebra: "+str(hebra))
            elif terminacion_nombre_query.startswith('aa'):
                hebra='+'
                #print("Dentro de clasificar dentro de  nucleotidos en if aa es positiva,  esto es hebra: "+str(hebra))
            else:
                hebra='.'
                marco_lectura='.'
                #print("Dentro de clasificar dentro de  aminoacidos no nucleotidos en else,  esto es hebra: "+str(hebra))
                #print("Dentro de clasificar dentro de  aminoacidos no nucleotidos en else,  esto es marco: "+str(marco_lectura))
            marco_lectura=int(terminacion_nombre_query[-1]) -1  # De la terminacion de aa y rev le acompana un numero que correspodne con el marco de lectura
            #print("Dentro de clasificar esto  es  marco lectura: "+str(marco_lectura))
            if  hebra=='+':
                inicio_nucleotido=(int(inicio_env2)-1)*3+marco_lectura+1
                fin_nucleotido=int(fin_env2)*3 +marco_lectura
                #print("Dentro de clasificar dentro de if strand + esto es inicio env: "+str(inicio_env2))
                #print("Dentro de clasificar dentro de if strand + esto es fin emv: "+str(fin_env2))
                #print("Dentro de clasificar dentro de if strand + esto es inicio nucleotido: "+str(inicio_nucleotido))
                #print("Dentro de clasificar dentro de if strand + esto es fin nucleotido: "+str(fin_nucleotido))
            elif hebra=='-':
                longitud_secuencia_nucleotido=diccionario_id_longitud_secuencia3[identificacion_query2]
                inicio_nucleotido=longitud_secuencia_nucleotido - (int(fin_env2) *3+marco_lectura)+1
                fin_nucleotido=longitud_secuencia_nucleotido -((int(inicio_env2)-1)*3 +marco_lectura)
                #print("Dentro de clasificar dentro de elf strand - esto es inicio env: "+str(inicio_env2))
                #print("Dentro de clasificar dentro de elif strand - esto es fin evn: "+str(fin_env2))
                #print("Dentro de clasificar dentro de elif strand - esto es longitud secuencia nucleotido: "+str(longitud_secuencia_nucleotido))
                #print("Dentro de clasificar dentro de elif esto es inicio nucleotido: "+str(inicio_nucleotido))
                #print("Estoy dentro de hmm2best dentro de elif esto es nuc_end: "+str(fin_nucleotido))
            else: # si no se realizo la traduccion
                inicio_nucleotido=inicio_env2
                fin_nucleotido=fin_env2
        elif tipo_secuencia_escogida_NuclProt3=='Proteinas': #Si la secuencia de entrada es de proteinas
            #print("Dentro de clasificar dentro de elif la secuencia de entrada son proteinas")
            diccionario_id_longitud_secuencia3=None  #Se coloca el diccionario como ninguno
            hebra='+'
            marco_lectura='.'
            inicio_nucleotido=inicio_env2
            fin_nucleotido=fin_env2
        #---Busco patron en identificacion del query
        patron3=re.compile(r'(\S+?):(\d+)[\.\-]+(\d+)')
        resultado_verif_patron3=patron3.match(identificacion_query2)  #Busco patron en identificacion del query
        #print("Dentro de clasificar esto es identificacion query: "+str(identificacion_query2))
        #print("Dentro de clasificar esto ees el patron: "+str(patron3))
        #print("Dentro de clasificar esto es resultado verif patron: "+str(resultado_verif_patron3))
        if resultado_verif_patron3: #Si existe coincidencia con el patron
            identificacion_query2,inicio_linea,fin_ltr=resultado_verif_patron3.groups()
            inicio_linea=int(inicio_linea)
            inicio_nucleotido=inicio_linea+inicio_nucleotido-1
            fin_nucleotido=inicio_linea+fin_nucleotido-1
            #print("Dentro de clasificar, dentro de coincidencia        clasificar_segundo_paso_genoma(nombre_ruta_archivo_gff2)de patron")
        espacio1=''
        linea1_archivo_gff=(identificacion_query2,'ARTHUR_LTRanalizer','CDS',inicio_nucleotido,fin_nucleotido,puntaje_secuencia2,hebra,marco_lectura, )
        #print("Dentro de clasificar esto es linea archivo gff:"+str(linea1_archivo_gff))
        if tipo_secuencia3=="Genomicas":
            identificacion_queryTarget_gff='{}:{}-{}|{}'.format(identificacion_query2,inicio_nucleotido,fin_nucleotido,nombre_target2)  #Coloco la identificacion para el archivo gff
            #print("Dentro de clasificar dentro de secuencias genomicas esto es identificacion del archivo gff: "+str(identificacion_queryTarget_gff))
            #----Obtener diccionario con  elementos  de la linea 2 para el archivo gff
            elementos_adicionales_linea_archivo_gff=({'ID':identificacion_queryTarget_gff,'Dominio':dominio5,'Clado':clado5},) 
            linea_completa_archivo_gff=linea1_archivo_gff+({'ID':identificacion_queryTarget_gff,'Dominio':dominio5,'Clado':clado5},) 
            #print("Dentro de cladsificar, en secuencias genomicas, esto es linea completa archivo gff: "+str(linea_completa_archivo_gff)) 
            identificacion_ltr='|'.join(identificacion_queryTarget_gff.split('|')[:-1])
            nombre5=identificacion_queryTarget_gff.split('|')[-1].split(':')[0]
            #print("Dentro  de clasificar, dentro de secuencias genomicas,  esto ees identificacion  ltr:  "+str(identificacion_ltr))
            #print("Dentro  de clasificar, en secuencias genomicas, esto es nombre5:  "+str(nombre5))  
            #----Determinar orden, superfamilia,clado
            if nombre_ref_dominios3=='REXdb':
                #print("Dentro de clasificar, dentro de secuencias genomicas, en Rexdb")
                orden3,superfamilia3,clado_mayor_frecuencia3,codificacion3=determinar_orden_superfamilia_clado_Rexdb([dominio5],[nombre5])
            elif nombre_ref_dominios3=='Gydb':
                orden3,superfamilia3,clado_mayor_frecuencia3,codificacion3=determinar_orden_superfamilia_clado_Gydb([dominio5],[clado5],nombre_ruta_Ref_Gydb3)
            else:
                orden3,superfamilia3,clado_mayor_frecuencia3,codificacion3='Desconocido','Desconocido','Desconocido','Desconocido'
            if orden3=='Desconocido':
                #print("Dentro de clasificar, en secuencias genomicas, elemento deconocido es exlido: "+str(identificacion_queryTarget_gff))
                continue
            #----Formato para archivo de salida cls, de etiqueta
            lista_nombres_para_cls=[orden3,superfamilia3,clado_mayor_frecuencia3]
            lista_nombres=[]
            for  nombre4 in lista_nombres_para_cls:
                if nombre4=='Desconocido' or nombre4 in  set(lista_nombres):
                    continue
                lista_nombres+=[nombre4]
            lista_nombres_paraCls_final='/'.join(lista_nombres)
            #print("Dentro de clasificar, en secuencias genomicas, eto es listanombres para archivo Cls: "+str(lista_nombres_paraCls_final))
            numero_parada=list(secuencia_archivo_gff).count('*')
            nombre_target_inicio_fin='{} {} {}'.format(nombre_target2,inicio_hmm_query2,fin_hmm_query2)
            espacio1='Clasificacion={};Target={};Nparada={};'.format(lista_nombres_paraCls_final,nombre_target_inicio_fin,numero_parada)
        nombre6='{}-{}'.format(clado5,dominio5)
        #print("Dentro de clasificar esto es nombre: "+str(nombre6))
        elementos_adicionales_linea_gff='ID={};Nombre={};{}Dominio={};Clado={};Cobertura={};Valor_e={};Probabilidad={}'.format(
            identificacion_queryTarget_gff,nombre6,espacio1,dominio5,clado5,cobertura2,valor_e_secuencia2,precision2)
        #print("Dentro de clasificar esto es elementos adcionales lineea gff: "+str(elementos_adicionales_linea_gff))
        linea1_archivo_gff=linea1_archivo_gff+(elementos_adicionales_linea_gff,valor_e_secuencia2,cobertura2,precision2,solo_identificacion_query,identificacion_queryTarget_gff,secuencia_archivo_gff)
        #print("Dentro de clasificar esto es linea1 archivo gff: "+str(linea1_archivo_gff))
        lista_lineas_gff.append(linea1_archivo_gff)
    #-----Creo los nombres con rutas para los archivos de salida: gff3, amin que tiene las secuencias aminoacidos, tsv archivo separado por tabs
    palabra5=nombre_ref_dominios3
    nombre_ruta_archivo_gff='%s/%s.%s.dominios.gff3' % (ruta_nuevacarpeta_salida3,nombre_archivo_entrada3,palabra5)
    nombre_ruta_archivo_seq='%s/%s.%s.dominios.amin' % (ruta_nuevacarpeta_salida3,nombre_archivo_entrada3,palabra5)
    nombre_ruta_archivo_tsv='%s/%s.%s.dominios.tsv' % (ruta_nuevacarpeta_salida3,nombre_archivo_entrada3,palabra5)
    print("Dentro de clasificar, esto es nombre ruta archivo gff: "+str(nombre_ruta_archivo_gff))
    print("Dentro de clasificar, esto es nombre ruta archivo de secuencias: "+str(nombre_ruta_archivo_seq))
    print("Dentro de clasificar, esto es nombre ruta archivo tsv: "+str(nombre_ruta_archivo_tsv))
    #print("Dentro de clasiicar, en secuencias esto es lista lineas gff: "+str(lista_lineas_gff))
    #----Ordeno las secuencias en base a la identificacion del query
    if tipo_secuencia3=="Genomicas":
        lista_lineas_gff=sorted(lista_lineas_gff,key=lambda i:i[0]) #Ordena elementos en base a la identificacion del query
        #print("Dentro de clasificar, en secuencias genomicas esto es lista lineas gff ordenadas: "+str(lista_lineas_gff))
        #----Resolucion de sobrelapes
        lista_lineas_sin_sobrelapes=[]
        for cromosoma,elementos in itertools.groupby(lista_lineas_gff,key=lambda i:i[0]):
            #print("Dentro de clasificar, en secuencias genomicas, resolviendo sobrelapes")
            lista_lineas_sin_sobrelapes+=corregir_sobrelapes(list(elementos))
        lista_lineas_gff=lista_lineas_sin_sobrelapes
    else:
        lista_lineas_gff=sorted(lista_lineas_gff,key=lambda j:(j[0],j[-3],j[3]))
        #print("Dentro de clasificar, en else esto es lista lineas gff ordenadas: "+str(lista_lineas_gff))
    #----Creo los archivos y los abro, con nombres de identificacion
    identificacion_archivo_gff=open(nombre_ruta_archivo_gff,'w')
    identificacion_archivo_seq=open(nombre_ruta_archivo_seq,'w')
    identificacion_archivo_tsv=open(nombre_ruta_archivo_tsv,'w')
    #------Coloco la informacion de titulos dentro del archivo tsv, separado por tabs
    print('\t'.join(['#Id','Longitud','Valor-e','Cobertura','Probabilidad','Puntaje']),file=identificacion_archivo_tsv)
    for linea5 in  lista_lineas_gff:
        linea1_archivo_gff=linea5[:9]
        #print("Dentro de clasificar, esto es linea para publickar en archivo gff: "+str(linea1_archivo_gff))
        linea1_archivo_gff=list(map(str,linea1_archivo_gff))
        #---Coloco la informacion dentro del archivo gff3
        print('\t'.join(linea1_archivo_gff),file=identificacion_archivo_gff)
        identificacion_queryTarget_gff,secuencia_archivo_gff=linea5[-2:]
        elementos_adicionales_archivo_seq=linea5[8]
        #print("Dentro de clasificar eesto es identificacion query target para archivo  secuencias: "+str(identificacion_queryTarget_gff))
        #print("Dentro de clasificar, esto es secuencia para archivos secuencias; "+str(secuencia_archivo_gff))
        #print("Dentro de clasificar, esto es elementos adicionales para publicar en archivo secuencias: "+str(elementos_adicionales_archivo_seq))
        #----Coloco la informacion dentro del archivo de secuencias
        print('>{} {}\n{}'.format(identificacion_queryTarget_gff,elementos_adicionales_archivo_seq,secuencia_archivo_gff),file=identificacion_archivo_seq)
        valor_e_secuencia2,cobertura2,precision2=linea5[-6:-3]
        puntaje_secuencia2=linea1_archivo_gff[5]
        linea5=[identificacion_queryTarget_gff,len(secuencia_archivo_gff),valor_e_secuencia2,cobertura2,precision2,puntaje_secuencia2]
        #print("Dentro de clasificar esto es linea 5 para colcoar en archivo tsv: "+str(linea5))
        #------Coloco la informacion dentro del archivo tsv
        print('\t'.join(map(str,linea5)),file=identificacion_archivo_tsv)
    #----Cerrar archivos gff y seq de secuencias
    identificacion_archivo_gff.close()
    identificacion_archivo_seq.close()
    print("Estoy dentro de clasificar al final")
    return nombre_ruta_archivo_gff,nombre_ruta_archivo_seq

#-----Realiza el segundo paso de clasificacion
def clasificar_segundo_paso(nombre_ruta_archivo_gff2,nombre_ref_dominios4,identificacion_archivo_clasificacion_tsv4,nombre_ruta_Ref_Gydb4):
    if nombre_ref_dominios4=='REXdb':
        dominios={'GAG','PROT','INT','RT','RH'}
    elif nombre_ref_dominios4=='Gydb':
        dominios={'GAG','AP','INT','RT','RNaseH','ENV'}
    #-----Coloco la primera linea de enunciado en archivo de clasificacion tsv separado por tabs
    linea1_archivo_clasif_tsv=['#TE','Orden','Superfamilia','Clado','Completo','Hebra','Dominio']
    print("Dentro de clasificar segundo paso, esto es identficacion achivo clasif tsv: "+str(identificacion_archivo_clasificacion_tsv4))
    print('\t'.join(linea1_archivo_clasif_tsv),file=identificacion_archivo_clasificacion_tsv4)
    #------Creo una lista que agrupa a los elementos en base a su identificacion
    lista_agupacion_lineas_por_id2=agrupar_elementos(nombre_ruta_archivo_gff2)
    #print("Dentro de clasificar segundo paso, esto es lista agrupacion lineas por id: "+str(lista_agupacion_lineas_por_id2))
    #------Determinar hebra, lista dominios
    for agrupacion_lineas_por_id in lista_agupacion_lineas_por_id2:
        agrupacion_lineas_por_id2=agrupacion_lineas_por_id
        #print("Dentro de clasificar segundo paso, dentro de for, esto es agrupacion lineas por id: "+str(agrupacion_lineas_por_id))
        #-----Creo una lista de hebras
        lista_hebras=[]
        for linea6 in agrupacion_lineas_por_id2:
            if isinstance(linea6,str):
                lista_linea6=linea6.strip().split('\t')
                #print("Dentro de clasificar segundo paso, en for, en if es cadena, eesto es lista_lineea6: "+str(lista_linea6))
            lista_hebras+=lista_linea6[6]
        #print("Dentro de clasificar segundo paso, dentro de for, esto es lista hebras: "+str(lista_hebras))
        numero_hebras=len(set(lista_hebras))
        #print("Dentro de clasificar segundo paso, dentro de for, esto es numero de hebras: "+str(numero_hebras))
        if numero_hebras>1:
            hebra4='?'
            #print("Dentro de clasificar segundo paso, en for, en if hebra >1 esto es hebra: "+str(hebra4))
        elif numero_hebras==1:
            hebra4=lista_hebras[0]
            #print("Dentro de clasificar segundo paso, en for, en elif dehebra =1 esto es hebra: "+str(hebra4))
        else:
            continue
        if hebra4=='-':
            agrupacion_lineas_por_id2.reverse()
            #print("Dentro de clasificar segundo paso, en hebra -, esto es agrupacion lineas por id2 reversa: "+str(agrupacion_lineas_por_id2))
        #---Obtener identificacion de la primera linea de la agrupacion
        primera_linea_deagrupacion=agrupacion_lineas_por_id2[0]
        if isinstance(primera_linea_deagrupacion,str):
                lista_primera_linea_deagrupacion=primera_linea_deagrupacion.strip().split('\t')
                #print("Dentro de clasificar segundo paso, en if es cadena, esto es lista_primera linea de agrupacion: "+str(lista_primera_linea_deagrupacion))
        elementos_adicionales4=lista_primera_linea_deagrupacion[8]
        #print("Dentro de clasificar segundo paso, dentro de for, esto es elementos adicionales de la primera linea: "+str(elementos_adicionales4))
        diccionario_elementos_adicionales_linea=dict(elemento2.split('=',1) for elemento2 in elementos_adicionales4.split(';') if elemento2)
        #print("Dentro de clasificar segundo paso, esto es diccionario con elementos adicionales linea: "+str(diccionario_elementos_adicionales_linea))
        identificacion_ltr7='|'.join(diccionario_elementos_adicionales_linea['ID'].split('|')[:-1])
        #print("Dentro  de clasificar segundo paso, dentro de segundo for, esto es identificacion ltr7: "+str(identificacion_ltr7))
        #-------Obtener lista de dominios
        lista_dominios7=[]
        lista_clados7=[]
        lista_nombres7=[]
        lista_dominios_clados=[]
        for linea7 in agrupacion_lineas_por_id:
            if isinstance(linea7,str):
                lista_linea7=linea7.strip().split('\t')
                #print("Dentro de clasificar segundo paso, en for dominios, en if es cadena, esto es lista de linea7: "+str(lista_linea7))
            elementos_adicionales5=lista_linea7[8]
            #print("Dentro de clasificar segundo paso, dentro de for dominios, esto es elementos adicionales: "+str(elementos_adicionales5))
            diccionario_elementos_adicionales_linea3=dict(elemento3.split('=',1) for elemento3 in elementos_adicionales5.split(';') if elemento3)
            #print("Dentro de clasificar segundo paso, en for dominios, esto es diccionario con elementos adicionales linea: "+str(diccionario_elementos_adicionales_linea3))
            dominio6=diccionario_elementos_adicionales_linea3['Dominio']
            clado6=diccionario_elementos_adicionales_linea3['Clado']
            nombre6=diccionario_elementos_adicionales_linea3['ID'].split('|')[-1].split(':')[0]
            #a2=''.join(['{}|{}'.format(dominio6,clado6)])
            a2=['{}|{}'.format(dominio6,clado6)]
            #print("Dentro de clasificar segundo paso, esto es a2: "+str(a2))
            lista_dominios_clados+=a2
            #print("Dentro de clasificar, en for dominios, esto es lista de dominios: "+str(lista_dominios_clados))
            #----Obtener lista dominios
            lista_dominios7+=[dominio6]
            lista_clados7+=[clado6]
            lista_nombres7+=[nombre6]
        lista_dominios_clados=' '.join(lista_dominios_clados)
        #print("Dentro de clasificar sgundo paso, despues for dominios, esto es lista de dominios_clados: "+str(lista_dominios_clados))
        #print("Dentro de clasificar segundo paso, despues for dominios, esto es lista de dominios: "+str(lista_dominios7))
        #print("Dentro de clasificar segundo paso, despues for dominios, esto es lista de cladoos: "+str(lista_clados7))     
        #print("Dentro de clasificar segundo paso, despues for dominios, esto es lista de nombres: "+str(lista_nombres7)) 
        if nombre_ref_dominios4=='REXdb':
            orden2,superfamilia2,clado_maximo2,codificacion2=determinar_orden_superfamilia_clado_Rexdb(lista_dominios7,lista_nombres7)
        elif nombre_ref_dominios4=='Gydb':
            orden2,superfamilia2,clado_maximo2,codificacion2=determinar_orden_superfamilia_clado_Gydb(lista_dominios7,lista_clados7,nombre_ruta_Ref_Gydb4)
        else:
            orden2,superfamilia2,clado_maximo2,codificacion2='Desconocido','Desconocido','Desconocido','Desconocido'
        linea8=[identificacion_ltr7,orden2,superfamilia2,clado_maximo2,codificacion2,hebra4,lista_dominios_clados]  
        #----Coloco lineas en archivo de clasificacion tsv separado por tabs
        #print("Dentro de clasificar segundo paso, en for esto es orden2: "+str(orden2))
        #print("Dentro de clasificar segundo paso, en for esto es linea8: "+str(linea8))
        print('\t'.join(linea8),file=identificacion_archivo_clasificacion_tsv4)
        yield linea8

#----Clasificar continuacion
def clasificar_tercer_paso(listas_elementos_linea_clasificados6,identificacion_archivo_clasificacion_tsv6,nombre_ruta_archivo_seq6,nombre_ruta_archivo_clasificacion_tsv6,nombre_ruta_archivo_gff6,ruta_nueva_carpeta_salida6,nombre_archivo_entrada6,nombre_ref_dominios6,nombre_ruta_archivo_entrada6):
    diccionario_identif_elementos_clasificados=OrderedDict() #Creo un diccionario vacio
    for  lista_elementos_linea_clasif in listas_elementos_linea_clasificados6:
        diccionario_identif_elementos_clasificados[lista_elementos_linea_clasif[0]]=lista_elementos_linea_clasif
    print("En clasificar tercer paso, esto es diccionario identif elementos clasificados: "+str(diccionario_identif_elementos_clasificados))
    identificacion_archivo_clasificacion_tsv6.close()
    numero_elementos_clasificados=len(diccionario_identif_elementos_clasificados)
    #-----Escribiendo  libreria para RepetMasker en archivo clasificacion libreria
    palabra5=nombre_ref_dominios6
    nombre_ruta_archivo_clasif_libreria='%s/%s.%s.clasificacion.libr' % (ruta_nueva_carpeta_salida6,nombre_archivo_entrada6,palabra5)
    identificacion_archivo_clasif_libreria=open(nombre_ruta_archivo_clasif_libreria,'w')
    for componentes_secuencia in SeqIO.parse(open(nombre_ruta_archivo_entrada6),'fasta'):
        identificacion_secuencia=componentes_secuencia.id 
        #print("En Clasificar Tercer paso, esto es el id: "+str(identificacion_secuencia))
        if componentes_secuencia.id in diccionario_identif_elementos_clasificados:
            #print("En Clasificar Tercer paso, dentro de if, este id esta en diccionario: "+str(identificacion_secuencia))
            elementos=diccionario_identif_elementos_clasificados[componentes_secuencia.id]
            #print("En clasificar tercer paso, dentro de if id esta en diccionario, esto es elementos: "+str(elementos))
            hebra6=elementos[5]
            #print("En clasificar tercer paso, dentro de if id esta en diccionario, esto es hebra: "+str(hebra6))
            lista_orden_superfamilia_clado=[elementos[1],elementos[2],elementos[3]]
            lista_elementos=[]
            for elemento in lista_orden_superfamilia_clado:
                if elemento=='Desconocido' or elemento=='Desconocida' or elemento in set(lista_elementos):
                    continue
                lista_elementos+=[elemento]
            elementos='/'.join(lista_elementos)
            #print("En clasificar tercer paso, dentro de if id esta en diccionario,   esto es elementos: "+str(elementos))
            if hebra6=='-':
                componentes_secuencia.seq=componentes_secuencia.seq.reverse_complement()
                #print("En clasificar Tercer paso, en if hebra -, esto es secuencia:"+str(secuencia))
        else:
            elementos='Desconocido'
        componentes_secuencia.id=componentes_secuencia.id.split('#')[0] + '#' + elementos
        #print("En clasifica tercer paso, esto es identificacion secuencia: "+str(componentes_secuencia.id))
        SeqIO.write(componentes_secuencia,identificacion_archivo_clasif_libreria,'fasta')
    identificacion_archivo_clasif_libreria.close()
    #---Escribir en archivo clasificacion las secuencias
    nombre_ruta_archivo_clasif_amin='%s/%s.%s.clasificacion.amin' % (ruta_nueva_carpeta_salida6,nombre_archivo_entrada6,palabra5)
    identificacion_archivo_clasif_amin=open(nombre_ruta_archivo_clasif_amin,'w')
    for componente_secuencia2 in SeqIO.parse(nombre_ruta_archivo_seq6,'fasta'):
        identificacion2='|'.join(componente_secuencia2.id.split('|')[:-1])
        elementos2=diccionario_identif_elementos_clasificados[identificacion2]
        lista_orden_superfamilia_clado2=[elementos2[1],elementos2[2],elementos2[3]]
        lista_elementos2=[]
        for elemento2 in lista_orden_superfamilia_clado2:
            if elemento2=='Desconocido' or elemento2=='Desconocida' or elemento2 in set(lista_elementos2):
                continue
            lista_elementos2+=[elemento2]
            elementos2='/'.join(lista_elementos2)
            #print("En clasificar tercer paso, dentro de if id esta en diccionario,   esto es elementos: "+str(elementos2))
        diccionario_descripcion=dict([pareja_elementos.split('=',1)for pareja_elementos in componente_secuencia2.description.split()[-1].split(';')])
        dominio2=diccionario_descripcion['Dominio']
        clado2=diccionario_descripcion['Clado']
        identificacion_modificada='{}#{}#{}|{}'.format(identificacion2.split('#')[0],elementos2,dominio2,clado2)
        componente_secuencia2.id=identificacion_modificada
        SeqIO.write(componente_secuencia2,identificacion_archivo_clasif_amin,'fasta')
    identificacion_archivo_clasif_amin.close()
    #-----Resumen del proces de clasificacion
    diccionario_resumen={}
    for  identificacion_resumen,elementos_resumen  in diccionario_identif_elementos_clasificados.items():
        #print("En clasificar tercer paso, en for resumen, esto es elementos resumen: "+str(elementos_resumen))
        clave_resumen=(elementos_resumen[1],elementos_resumen[2])
        diccionario_resumen[clave_resumen]=[0,0,[],0] #inicializo el diccionario
    for  identificacion_resumen,elementos_resumen  in diccionario_identif_elementos_clasificados.items():
        clave_resumen=(elementos_resumen[1],elementos_resumen[2])
        diccionario_resumen[clave_resumen][0]+=1
        if elementos_resumen[3] not in {'Desconocido','Mezcla'}:
            diccionario_resumen[clave_resumen][1]+=1
            diccionario_resumen[clave_resumen][2]+=[elementos_resumen[3]]
        if elementos_resumen[4]=='Si':
            diccionario_resumen[clave_resumen][3]+=1
    #print("Dentro de clasificar tercer paso, despues de for, esto es diccionario resumen: "+str(diccionario_resumen))
    formato2='{:^17}{:^17}{:^30}{:^17}{:^30}{:^14}{:^26}'.format('Orden','Superfamilia','Clado','# De  Secuencias','# De Secuencias con Clados','# De Clados','# De Dominios Completos') # Dentro de llaves se coloca el numero de caracteres y alineacion, mayor es alienacion a izq, menor a der
    #print("En clasificar 3, esto es formato2: "+str(formato2))
    nombre_ruta_archivo_resumen='%s/%s.%s.resumen' % (ruta_nueva_carpeta_salida6,nombre_archivo_entrada6,palabra5)
    identificacion_archivo_resumen=open(nombre_ruta_archivo_resumen,'w')
    print(formato2,file=identificacion_archivo_resumen)  
    lista_ordenes=['LTR','pararetrovirus','DIRS','Penelope','LINE','SINE','TIR','Helitron','Maverick','mezcla','Desconocido','Total']
    items_ordenados=sorted(list(diccionario_resumen.items()), key=lambda k: (lista_ordenes.index(k[0][0]), k[0][1])) # ordena en base a orden y superfamilia y despues en base al numero de elementos identificados
    #print("En clasificar tercer paso, esto es items ordenadods:"+str(items_ordenados))
    for (orden3,superfamilia3), resumen3 in items_ordenados:
        numero_clados=len(set(resumen3[2]))
        linea_resumen2=[orden3,superfamilia3,resumen3[2],resumen3[0],resumen3[1],numero_clados,resumen3[3]]  
        linea_resumen2=list(map(str,linea_resumen2))
        linea_resumen2='{:^17}{:^17}{:^30}{:^17}{:^30}{:^14}{:^26}'.format(*linea_resumen2)
        print(linea_resumen2,file=identificacion_archivo_resumen)
    print("En clasificar tercer paso, el numero de elementos clasiicados es: "+str(numero_elementos_clasificados))
    print("En clasificar tercer paso, las secuencias de dominios proteicos se ubican en : "+str(nombre_ruta_archivo_seq6))
    print("En clasificar tercer paso, el archivo de anotacion gff se ubica en : "+str(nombre_ruta_archivo_gff6))
    print("En clasificar tercer paso, las secuencias clasificadas se ubican en : "+str(nombre_ruta_archivo_clasificacion_tsv6))
    print("En clasificar tercer paso, la libreria para RepeatMasker se ubica en: "+str(nombre_ruta_archivo_clasif_libreria))
    print("En clasificar tercer paso, las secuencias clasificadas con nueva identificacion se ubican en: "+str(nombre_ruta_archivo_clasif_amin))
    print("En clasificar tercer paso, el resumen de la clasificacion se ubica en: "+str(nombre_ruta_archivo_resumen))
    print("EN CLASIFICAR3, EL PROCESO HA TERMINADO")

#-----Generar resumen
def clasificar_segundo_paso_genoma(nombre_ruta_archivo_gff4,ruta_nueva_carpeta_salida7,nombre_archivo_entrada7,nombre_ref_dominios7):
    valor_final=0
    dicccionario_elementos_resumen={}
    for linea11 in open(nombre_ruta_archivo_gff4):
        #print("Dentro  clasificar segundo paso genomoa, esto es linea11: "+str(linea11))
        if isinstance(linea11,str): # Si la  linea es una cadena de caracteres
            lista_elementos_linea_gff4=linea11.strip().split('\t') #Se convierte la cadena de caracteres separada por tabs en una lista
            #print("Dentro de clasificar segundo paso genoma, en linea es un caracter, esto es la lista de elementos de linea gff: "+str(lista_elementos_linea_gff4))
        else:
            lista_elementos_linea_gff4=linea11
            #print("Dentro de clasificar segundo paso genoma, en linea no es un caracter, esto es la lista de elementos de linea gff: "+str(lista_elementos_linea_gff4))
        if len(lista_elementos_linea_gff4)==8:
            lista_elementos_linea_gff4=list(lista_elementos_linea_gff4)+['']
            #print("Dentro de clasificar segundo paso genoma, dentro de if longitud =8")
        #Obtengo elementos de cada linea del archivo gff
        identificacion_query4,fuente_analisis4,tipo_analisis4,inicio_nucleotido4,fin_nucleotido4,puntaje_secuencia4,hebra4,marco_lectura4,elementos_adicionales4=lista_elementos_linea_gff4
        diccionario_elementos_adicionales_linea_gff4=dict(elemento4.split('=',1) for elemento4 in elementos_adicionales4.split(';') if elemento4)
        #print("Dentro de clasificar segundo paso genoma, esto es diccionario con elementos adicionales linea gff4: "+str(diccionario_elementos_adicionales_linea_gff4))
        #dominio6=diccionario_elementos_adicionales_linea_gff['Dominio']
        #clado6=diccionario_elementos_adicionales_linea_gff['Clado']
        clasificacion4=diccionario_elementos_adicionales_linea_gff4['Clasificacion']
        #print("En clasificar segundo paso genoma, esto es clasificacion4 antes:"+str(clasificacion4))
        clasificacion4=tuple(clasificacion4.split('/'))
        #print("En clasificar segundo paso genoma, esto es clasificacion4 despues:"+str(clasificacion4))
        #print("En clasificar segundo paso genoma,  esto es elemento 0 de clasificaion4; "+str(clasificacion4[0]))
        ordenes4=['LTR','Maverick','pararetrovirus','DIRS','LINE','TIR','Penelope','SINE','Helitron','Total','mezcla','Desconocido']
        assert len(clasificacion4)<=4  #3
        assert clasificacion4[0] in set(ordenes4),'Orden  Desconocido: {}'.format(clasificacion4)
        if int(fin_nucleotido4)<valor_final:
            continue
        elif int(inicio_nucleotido4)<valor_final:
            inicio4=valor_final+1
        else:
            inicio4=int(inicio_nucleotido4)
        longitud4=int(fin_nucleotido4)-inicio4+1
        lista_clasificacion4=[]
        if len(clasificacion4)==3:
            #print("En clasificar segundo paso en longitud 3 de ty/copia, no hace nada, esto es clasificacion4: "+str(clasificacion4))
            lista_clasificacion4+=[clasificacion4[:1],clasificacion4[:2]]
        elif len(clasificacion4)==2:
            #print("En clasificar segundo paso en longitud 2 de ty/copia, no hace nada, esto es clasificacion4: "+str(clasificacion4))
            lista_clasificacion4+=[clasificacion4[:1]]
        elif len(clasificacion4)==4:
            print("En clasificar segundo paso en longitud 4 de ty/copia, no hace nada, esto es clasificacion4: "+str(clasificacion4))
        if len(clasificacion4)<4:
            #print("En clasificar segundo paso en <4 de ty/copia, no hace nada, esto es clasificacion4: "+str(clasificacion4))
            lista_clasificacion4+=[clasificacion4]+[('Total', )]
        #print("En clasificar segundo paso genoma, esto es lista clasificacion4: "+str(lista_clasificacion4))
        for clasif  in lista_clasificacion4:
            #print("En clasificar segundo paso genoma, entre por el for")
            try:dicccionario_elementos_resumen[clasif]+=[longitud4]
            except KeyError:dicccionario_elementos_resumen[clasif]=[longitud4]
    
    linea_elementos4=['Orden','Superfamilia','Clado','Numero','Longitud_Total','Longitud_media']
    linea_elementos4='{:^15}{:^22}{:^32}{:^16}{:^24}{:^24}'.format(*linea_elementos4)
    palabra7=nombre_ref_dominios7
    nombre_ruta_archivo_resumen4='%s/%s.%s.resumen' % (ruta_nueva_carpeta_salida7,nombre_archivo_entrada7,palabra7)
    identificacion_archivo_resumen4=open(nombre_ruta_archivo_resumen4,'w')
    print(linea_elementos4,file=identificacion_archivo_resumen4)
    #identificacion_archivo_resumen4.write('\t'.join(linea_elementos4)+'\n')
    for clasificacion5,longitud5 in sorted(dicccionario_elementos_resumen.items(), key=lambda x:(ordenes4.index(x[0][0]),x[0])):
        clasificacion5=['']*(len(clasificacion5)-1)+[clasificacion5[-1]]+['']*(3-len(clasificacion5))
        #print("En clasificar segundo paso genoma, esto es clasificacion5: "+str(clasificacion5))
        numero5,total5=len(longitud5),sum(longitud5)
        media5=round(total5/numero5,1)
        linea_elementos4=clasificacion5+[numero5,total5,media5]
        #print("En  clasificacion segundo  paso  genoma esto es linea  elementos4: "+str(linea_elementos4))
        linea_elementos4=map(str,linea_elementos4)
        #print("En  clasificacion segundo  paso  genoma esto es linea  elementos4: "+str(linea_elementos4))
        #identificacion_archivo_resumen4.write('\t'.join(linea_elementos4)+'\n')
        linea_elementos4='{:^15}{:^22}{:^32}{:^16}{:^24}{:^24}'.format(*linea_elementos4)
        print(linea_elementos4,file=identificacion_archivo_resumen4)
    print("PROCESO FINALIZADO")


#-----Creo una lista que agrupa a los elementos en base a su identificacion
def agrupar_elementos(nombre_ruta_archivo_gff2_2):
    lista_agupacion_lineas_por_id=[]
    caracter_final=''
    for linea5 in open(nombre_ruta_archivo_gff2_2):
        #print("Dentro de agrupar elementos, esto es linea: "+str(linea5))
        if linea5.startswith('#'):
                continue
        if isinstance(linea5,str): # Si la  linea es una cadena de caracteres
            lista_elementos_linea_gff=linea5.strip().split('\t') #Se convierte la cadena de caracteres separada por tabs en una lista
            #print("Dentro de agrupar, en linea es un caracter, esto es la lista de elementos de linea gff: "+str(lista_elementos_linea_gff))
        else:
            lista_elementos_linea_gff=linea5
            #print("Dentro de agrupar, en linea no es un caracter, esto es la lista de elementos de linea gff: "+str(lista_elementos_linea_gff))
        if len(lista_elementos_linea_gff)==8:
            lista_elementos_linea_gff=list(lista_elementos_linea_gff)+['']
            #print("Dentro de agrupar, dentro de if longitud =8")
        #Obtengo elementos de cada linea del archivo gff
        identificacion_query3,fuente_analisis3,tipo_analisis3,inicio_nucleotido3,fin_nucleotido3,puntaje_secuencia3,hebra3,marco_lectura3,elementos_adicionales3=lista_elementos_linea_gff
        #inicio_nucleotido3=int(inicio_nucleotido3)
        #fin_nucleotido3=int(fin_nucleotido3)
        #try: puntaje_secuencia3=float(puntaje_secuencia3)
        #except: pass
        #try: marco_lectura3=int(marco_lectura3)
        #except: pass
        diccionario_elementos_adicionales_linea_gff=dict(elemento.split('=',1) for elemento in elementos_adicionales3.split(';') if elemento)
        #print("Dentro de agrupar elementos, esto es diccionario con elementos adicionales linea gff: "+str(diccionario_elementos_adicionales_linea_gff))
        #dominio6=diccionario_elementos_adicionales_linea_gff['Dominio']
        #clado6=diccionario_elementos_adicionales_linea_gff['Clado']
        identificacion_ltr6='|'.join(diccionario_elementos_adicionales_linea_gff['ID'].split('|')[:-1])
        #nombre6=diccionario_elementos_adicionales_linea_gff['ID'].split('|')[-1].split(':')[0]
        if lista_agupacion_lineas_por_id and not identificacion_ltr6==caracter_final:
            yield lista_agupacion_lineas_por_id
            lista_agupacion_lineas_por_id=[]
        lista_agupacion_lineas_por_id.append(linea5)
        caracter_final=identificacion_ltr6
    #print("Dentro de agrupar, despues de for, esto es lista agrupacion lineeas por id: "+str(lista_agupacion_lineas_por_id))
    #print("Dentro de agrupar, despues de for, esto es diccionario con elementos adicionales linea gff: "+str(diccionario_elementos_adicionales_linea_gff))
    yield lista_agupacion_lineas_por_id

#------Multiprocesamiento, se paraleliza la ejecucion de una funcion que tiene multiples valores de entrada
def multiprocesar(funcion,iterable,**kargs):
    procesos=8
    multiproceso=multiprocessing.Pool(procesos)     #Defino el numero de procesosscomando Pool del paquete multiprocessing y le indico que trabjara el multiprocesamiento con 8 procesadores
    #print("versi estoy dentr de multiproceso_traduccion")
    objeto_iterable_resultados_multiproceso=multiproceso.imap(funcion,iterable,**kargs) #Guardo los resultados del multiproceso de realizar la funcion especificada
    for valor_de_retorno in objeto_iterable_resultados_multiproceso: #Itero cada valor de retorno dentro del objeto iterable resultados
        #print("Esto es valor de retorno de pool imAP:"+str(valor_de_retorno))
        yield valor_de_retorno
    multiproceso.close()     #Terminar el proceso de pool
    multiproceso.join()      #Para esperar que todo el proceso se termine


#------Idenfificar si es una secuencia en formato fasta or fastq, lee cada linea de archivo de entrada y separa por secuencias, el contenido de 1 secuencia agrupa en una lista
def identificar_agrupar_por_secuencia_fasta(nombre_ruta_Seq_In2): 
    j=0
    lista_contenido_una_secuencia=[]      #Lista en la que se almacenara cada secuencia
    for linea in nombre_ruta_Seq_In2:
        j+=1
        if j==1:
            if linea[0]=='>':     #Analiza que archivo de secuencias sean fasta porque inicia con >
                formato_de_secuencia='fasta'
                #print("estoy en el primer if y es fasta")
            elif linea[0]=='@':
                formato_de_secuencia='fastq'
            else:
                raise ValueError("Secuencia desconocida, no es un archivo fasta")
            lista_contenido_una_secuencia.append(linea)         #Aniado la primera linea a lista
            #print("Primer append",str(lista_contenido_una_secuencia))
            continue  #salta de nuevo al for para chequear la siguiente linea
        if (formato_de_secuencia=='fasta' and linea[0]=='>') or (formato_de_secuencia=='fastq' and j% 4==1):  #Si encuentra el inicio de una nueva secuencia, retorna toda la secuencia, en f fastq cada 4 lineas inicia una nueva secuencia 
            contenido_una_secuencia=''.join(lista_contenido_una_secuencia) #Join une todos los items de una lista y entrega en formato string, que una con un 1 espacio es:''
            yield contenido_una_secuencia  #Yield  retorna una secuencia completa, un grupo grande de valores 
            #print("estoy despues de yield"+str(lista_contenido_una_secuencia))
            lista_contenido_una_secuencia=[]               #Reinicio la lista para iniciar con una nueva secuencia
        lista_contenido_una_secuencia.append(linea)        #Se va aniadiendo las demas lineas a la lista 
        #print("soy externo de if"+str(lista_contenido_una_secuencia))
    contenido_una_secuencia=''.join(lista_contenido_una_secuencia) #Join une todos los items de una lista y entrega en formato string, que una con un 1 espacio es:''      
    #print("Esto es contenido de toda la secuencia: "+str(contenido_una_secuencia))
    yield contenido_una_secuencia                          #Retorna una secuencia, es para la ultima secuencia, yield retorna un grupo grande de valores

def traducir(nombre_ruta_archivo_fasta2):
#def traduccion(nombre_archivo_fasta2):
    #print("Dentro de traducir, esto es nombre ruta archivo fasta2:"+str(nombre_ruta_archivo_fasta2))
    #nombre_archivo_fasta_secuencia=nombre_archivo_fasta2   #El nombre de uno de los 8 archivos creados como:rchivo1.fasta
    #sobreescribir=True
    formato_secuencia2='fasta'
    tabla_traduccion=1
    #prefix=nombre_archivo_fasta_secuencia
    nombre_ruta_archivo_traducido=nombre_ruta_archivo_fasta2+'.amin'   #El nombre del archivo de salida traducido: archivo1.fasta.amin
    #sobreescribir2=not(os.path.exists(nombre_archivo_traducido) and os.path.getsize(nombre_archivo_traducido) >0) or sobreescribir
    #print("dentro de traduccion valor de sobreescribir ospath: "+str(sobreescribir2))
    #if not sobreescribir2:  #Si existe el archivo 
        #print("dentro de traduccion y dentro de ifnot sobreescribir")
        #return nombre_archivo_traducido
    #Creo y abro un archivo fasta traducido que contendra las secuencias fasta traducidas
    with open(nombre_ruta_archivo_traducido,'w') as identificacion_archivo_traducido_amin: #fp contiene informacion del archivo fasta.aa recien creado
        #print("dentro de traducir, dentro withopen este es nombre ruta archivo traducido"+str(identificacion_archivo_traducido_amin))
        #six_frame_translate(nombre_archivo_fasta_secuencia,informacion_archivo_traducido)
        #--Traduccion de secuencias con 6 marcos, se crean 6 archivos para cada secuencia, 3 archivos aa +1,+2,+3 y 3 archivos reversa -1,-2,-3
        longitud_secuencia_nucleotidos= OrderedDict()  #Creo un diccionario ordenado
        for secuencia3 in SeqIO.parse(open(nombre_ruta_archivo_fasta2),formato_secuencia2):   #Abro un archivo archivo1.fasta le pongo como objto SeqIO y Recorro cada linea
            for solo_secuencia,nombre1 in zip([secuencia3.seq,secuencia3.seq.reverse_complement()],['aa','rev_aa']):  #Metodo Zip crea una tupla entre secuencia y tipo, (secuencia normal,aa), (secuencia reversa,rev_aa)
                #print("este es el valor de seq de zip: "+str(solo_secuencia))
                #print("eeste es el valor de suffix: "+str(nombre1))
                for i in range(0,3):  #lazo que se repite 3 veces
                    #print("este es el numero de marco que hace: "+str(i))
                    secuencia_nucleotidos=solo_secuencia[i:]
                    #print("valor de secuencia de nucleotidos dentro de traduccion 6 marcos: "+str(secuencia_nucleotidos))
                    secuencia_aminoacidos=secuencia_nucleotidos.translate(table=tabla_traduccion) 
                    try: secuencia_aminoacidos=secuencia_nucleotidos.translate(table=tabla_traduccion) 
                    #try: secuencia_aminoacidos=translate_seq(secuencia_nucleotidos,table=tabla_traduccion)
                    except CodonTable.TranslationError:continue
                    nombre_secuencia_aminoacidos='|{}{}'.format(nombre1,i+1)
                    print('>{}{}\n{}'.format(secuencia3.id,nombre_secuencia_aminoacidos,secuencia_aminoacidos),file=identificacion_archivo_traducido_amin)
            longitud_secuencia_nucleotidos[secuencia3.id]=len(secuencia3.seq)

    return nombre_ruta_archivo_traducido

#------Se recortan secuencias genomicas
def recortar_secuencias_genomicas(nombre_direccion_archivo_entrada,identificacion_archivo_secuen_recortadas):
    tamanio_ventana=int(270000)  #270000
    sobrelape_ventana=int(30000) #30000
    #print("En funcion Recortar_secuencias_genomicas")
    #----Coloco  secuencias recortadas en un nuevo  archivo:secuencias recortadas.fasta
    for secuencia5 in SeqIO.parse(open(nombre_direccion_archivo_entrada),'fasta'):
        longitud_secuencia5=len(secuencia5.seq)
        for longitud_pedazo_secuencia in range(0,longitud_secuencia5+1,tamanio_ventana):
            valor1=longitud_pedazo_secuencia+tamanio_ventana+sobrelape_ventana
            if valor1>longitud_secuencia5:
                valor1=longitud_secuencia5
                valor2=valor1-longitud_pedazo_secuencia
                if valor2<tamanio_ventana/10:
                    longitud_pedazo_secuencia=max(0,longitud_pedazo_secuencia-tamanio_ventana/10)
            secuencia_recortada5=secuencia5.seq[longitud_pedazo_secuencia:valor1]
            identifiicacion_secuencia_recortada5='%s:%s-%s'%(secuencia5.id,longitud_pedazo_secuencia+1,valor1)
            identifiicacion_secuencia_recortada5+=' longitud=%s' % len(secuencia_recortada5)
            print('>%s\n%s' % (identifiicacion_secuencia_recortada5,secuencia_recortada5), file=identificacion_archivo_secuen_recortadas)
            if  valor1==longitud_secuencia5:
                continue

##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Escoger archivo de secuencias
#boton_escoger_archivo_entrada=Button(subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="ESCOGER ARCHIVO",padx=5,pady=2,width=17,command=escoger_archivo_entrada)
#boton_escoger_archivo_entrada.grid(row=0,column=2,sticky="w")
##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger tipo de secuencia
#tipo_secuencia_menu=OptionMenu(subsubmarco_inferior_left_superior,varmenu_tipo_secuencia,*opcionesmenu_tipo_secuencia) 
#tipo_secuencia_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
#tipo_secuencia_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
#tipo_secuencia_menu.grid(row=1,column=2,padx=0,pady=1,sticky=W)
##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger tipo de secuencia Nucleotidos o Proteinas
#tipo_secuencia_NuclProt_menu=OptionMenu(subsubmarco_inferior_left_superior,varmenu_tipo_secuencia_NuclProt,*opcionesmenu_tipo_secuencia_NuclProt) 
#tipo_secuencia_NuclProt_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
#tipo_secuencia_NuclProt_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
#tipo_secuencia_NuclProt_menu.grid(row=2,column=2,padx=0,pady=1,sticky=W)
##-----Creo Opcion desplegable en la ventana del marco subsubmarco_inferior_left, para Escoger archivo de referencia si es REXdb o Gydb
#baseref_dominios_menu=OptionMenu(subsubmarco_inferior_left_superior,varmenu_basereferencia_dominios,*opcionesmenu_basereferencia_dominios) 
#baseref_dominios_menu.config(width=11,font=("arial",14,"bold"),bg="white",fg="black",bd=4,relief=GROOVE,borderwidth=4,activebackground="lightcyan",activeforeground="blue",anchor=W,direction=RIGHT)
#baseref_dominios_menu['menu'].config(font=("arial",14)) #aumenta el tamanio de letra del menu desplegable
#baseref_dominios_menu.grid(row=3,column=2,padx=0,pady=1,sticky=W)
##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Inciar el proceso
#boton_inicio=Button(subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="INICIAR BUSQUEDA",padx=5,pady=2,width=17,command=iniciar_proceso)
#boton_inicio.grid(row=5,column=2,sticky="w")
##-----Creo Boton en la ventana del marco subsubmarco_inferior_left para Escoger ruta de salida en donde se guardaran loa archivos de salida
#boton_escoger_archivo_salida=Button(subsubmarco_inferior_left_superior,bd=4,font=("arial",11,"bold"),bg="darkgreen",fg="white",text="ESCOGER RUTA",padx=5,pady=2,width=17,command=escoger_ruta_salida)
#boton_escoger_archivo_salida.grid(row=4,column=2,sticky="w")

#----Configuracion para poder cambiar tamano de pantalla usando grid
#subsubmarco_inferior_left_superior.columnconfigure(0, weight=1)
#subsubmarco_inferior_left_superior.columnconfigure(1, weight=1)
#subsubmarco_inferior_left_superior.columnconfigure(2, weight=1)
#subsubmarco_inferior_left_superior.columnconfigure(3, weight=1)
#subsubmarco_inferior_left_superior.rowconfigure(0, weight=1)
#subsubmarco_inferior_left_superior.rowconfigure(1, weight=1)
#subsubmarco_inferior_left_superior.rowconfigure(2, weight=1)
#subsubmarco_inferior_left_superior.rowconfigure(3, weight=1)
#subsubmarco_inferior_left_superior.rowconfigure(4, weight=1)
#submarco_left_superior.rowconfigure(0, weight=1)
#submarco_left_superior.columnconfigure(0, weight=1)
#marco_inferior.rowconfigure(0, weight=1)
#marco_inferior.rowconfigure(1, weight=1)
#marco_inferior.rowconfigure(2, weight=1)
#marco_inferior.columnconfigure(0, weight=1)
#marco_botones.rowconfigure(0, weight=1)
#marco_botones.columnconfigure(0, weight=1)
#marco_superior.rowconfigure(0, weight=1)
#marco_superior.columnconfigure(0, weight=1)
#ventana.rowconfigure(0, weight=1)

#------Para poner 4 puntos en la parte inferior derecha de la pantalla para poder modificar tamano de pantalla
#puntos_extender_esquina = ttk.Sizegrip(ventana)   #Coloco puntos para  arrastrar diagonalmente y modificar tamano de pantalla
#puntos_extender_esquina.pack(expand = True, fill = BOTH, anchor = SE)

def main():
    ventana=Tk()
    interfaz_grafica(ventana)
    ventana.mainloop()   #Lazo loop para que la ventana de la interfaz grafica funcione continuamente
    print("en main al final antes de return 0")
    return 0

if __name__ == '__main__':
    #ventana=Tk()
    #ventana.mainloop()   #Lazo loop para que la ventana de la interfaz grafica funcione continuamente
    #main()
    print("en if name main")
    sys.exit(main()) 
#ArchivoSecuencias=input("Ingrese archivo que tiene secuencias")
#print(ArchivoSecuencias)
