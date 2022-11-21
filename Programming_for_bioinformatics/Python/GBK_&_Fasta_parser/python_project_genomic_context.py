#!/usr/bin/env python3

#Importamos los módulos necesarios para este script:

from Bio import SeqIO
import sys
import subprocess
import pandas as pd
import os
import glob

#Definimos la función de ayuda para que nos devuelva un
#mensaje en caso de que se de un número incorrecto de 
#argumentos o el frame esté fuera de rango

def ayuda():
	print ("""
Ejercicio optativo python v1.0
Gonzalo Cardenal Antolin, 2021

This script will run a blast analysis against a protein query
and, if there is at least one hit, it will print the genomic context,
with the locus_tag and product of each neighboring gene.


Usage:  ./script_name.py file.gff protein_query.fa [2-4]

file.gff         : a gff file containing the information of the living organism of interest
protein_query.fa : a fasta file containing protein query sequences for blast
[2-4]            : number of neighboring genes you would like to obtain.
	""")

#Definimos la función eliminar para borrar todos aquellos
#archivos creados por el script, tanto si tenemos hit
#como si no en el análisis de blast

def eliminar():
	if os.path.exists("blast_result.txt"):
		os.remove("blast_result.txt")
	if os.path.exists("result_filtered.txt"):
		os.remove("result_filtered.txt")

#Con el módulo de glob obtenemos una lista con todos los nombres
#que hagan match con un patrón. Nos ayuda a borrar los diferentes
#archivos que ha creado makeblastdb que empiezan por output.f(nuestro patrón).

	fileList = glob.glob("output.f*", recursive=True)

	for file in fileList:
		if os.path.exists(file):
			os.remove(file)


if sys.argv[1] == "-h":
        ayuda()
        sys.exit(0)

#Comprobamos que tenemos los argumentos necesarios (3 + 1(nombre))

if len(sys.argv) != 4:
	print ("Incorrect number of arguments\n")
	ayuda()
	sys.exit(0)

input_file = sys.argv[1]
query = sys.argv[2]
frame1 = sys.argv[3]

#Definimos el tercer argumento como un número, ya que 
#el argumento lo pasa como una str, creándonos problemas 
#posteriormente.

frame = int(frame1)

#Definimos la función is_fasta para comprobar que el archivo
#entregado como segundo argumento es un archivo con
#formato fasta

def is_fasta(query):
	with open(query, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")

		return any(fasta) #devuelve falso si no hay ninguna sec fasta

#Comprobamos que el archivo con la secuencia de la proteína existe y,
#de ser cierto, comprobamos que es fasta, si no, abortamos el script

if os.path.exists(query):
	f = is_fasta(query)
	if f == False :
		print ("File: " + query + " has no fasta format. Try again.")
		ayuda()
		sys.exit()
else:
	print ("File: " + query + " does not exist or it is not in the current directory")
	ayuda()
	sys.exit()

#Comprobamos que el archivo con formato gff al menos
#existe antes de correr el script. En caso contrario,
#que devuelva el mensaje de ayuda y se termina el script.

if not os.path.exists(input_file):
	print ("File: " + input_file + " does not exist or it is not in the current directory")
	ayuda()
	sys.exit()

#Comprobamos que el frame se encuentra en el 
#rango deseado [2-4]

if frame < 2 or  frame > 4:
	print ("Frame out of range\n")
	ayuda()
	sys.exit(0)

#Definimos la función para extraer del archivo gff todos
#los locus y crear un diccionario que une el locus
#con el producto y la cadena

def make_list_and_dict(input_file):

	list_of_genes = []
	dict_of_genes = {}

	with open(input_file, "r") as input_handle:
		for record in SeqIO.parse(input_handle, "genbank"):
			for feature in record.features:
				if feature.type == 'CDS':
					locus = feature.qualifiers['locus_tag'][0] 
					list_of_genes.append(locus)

					strand = feature.location.strand
					try:
						product = feature.qualifiers['product'][0] 
					except:
						product = "NA"
					dict_of_genes[locus]=[strand,product]


	return  dict_of_genes,list_of_genes

#Definimos la función extractprotein para extraer las secuencias de aquellos
#genes con tipo CDS que tengan traducción. Crearemos un archivo que contenga
#todas las secuencias de proteínas para poder hacer a partir del mismo
#la base de datos necesaria para poder hacer el análisis de Blast.


def extractprotein(input_file):

        list_of_genes = []
        dict_of_genes = {}

        output = open("output.fa","w")
        with open(input_file, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                        for feature in record.features:
                                if feature.type == 'CDS':
                                        locus = feature.qualifiers['locus_tag'][0]
                                        try:
                                                translation = feature.qualifiers['translation'][0]
                                        except:
                                                translation = "NA"
                                        output.write(">" + locus + "\n" + translation + "\n")
        output.close()

#Definimos la función para hacer el blast desde Python.
#Debido a que no se puede utilizar directamente las líneas de blast
#desde un script de python (solo desde terminal), podemos crear
#subprocesos para que se lleva a cabo el análisis de blast.


def blast():
	print ("\nCreating the database to run the local blast")

#Primero, creamos una base de datos con las secuencias del archivo creado
#anteriormente que contiene todas las sec de proteínas del archivo gff:

	p = subprocess.Popen(['makeblastdb', '-in', 'output.fa', '-dbtype', 'prot'])

#Debido a que los subprocesos se hacen en otro shell, es necesario esperar a que
#acaben para continuar. Para ello usamos la función poll(), que devuelve 'None'
#mientras no haya acabado

	while True:
		if p.poll() is not None:
			print ("\nRunning BLAST analysis and result filtering...")
			b = subprocess.Popen(['blastp', '-query', query ,'-db', 'output.fa', '-out', 'blast_result.txt', '-outfmt', '6 qseqid sseqid pident qcovs'])
			break
		else:
			continue
#Cuando acabe el subproceso del blast, empezamos con la filtración. Para ello,
#utilizaremos awk con otro subproceso
	while True:
		if b.poll() is not None:
			file = open('result_filtered.txt', 'w')
			cat = ['cat', 'blast_result.txt']
			awk = ['awk', 'BEGIN {OFS = "\t"} {if ($3 > 40 && $4 >= 51) {print $0}}']
			cmd_1 = subprocess.Popen(cat,stdout=subprocess.PIPE)
			cmd_2 =subprocess.Popen(awk, stdin = cmd_1.stdout, stdout=file)
			file.close()
			break
		else:
			continue

#Finalmente, leeremos el archivo filtrado y cogeremos el locus correspondiente
#a la primera línea, si que existe, ya que corresponde con el mejor de los resultados.
#Guardaremos la variable de locus para posteriormente printear el contexto genómico
#En caso de no tener ningún hit se informará y se terminará el script.

	while True:
		if cmd_2.poll() is not None:
			file = open("result_filtered.txt", "r")
			nonempty_lines = [line.strip("\n") for line in file if line != "\n"]
			line_count = len(nonempty_lines)
			file.close()
			if line_count != 0:
				f1 = open("result_filtered.txt", "r")
				line1=f1.readline()
				locus=(line1.split("\t")[1])
				f1.close()
			else:
				print ("\nUps...no coincidence in blast result. Try again.")
				eliminar()
				sys.exit ()

			break
		else:
			continue
	return locus, line1

#Definimos la función para obtener el locus de los genes vecinos.
#Para ello primeramente obtenemos la posición del locus de interés
#y posteriormente, y dependiendo del frame escogido [2-4], añadimos
#a la lista el locus de los genes vecinos al locus de interés


def get_genomic_neighbourhood_list(locus_tag, list_of_genes, frame):


	genomic_neighbourhood_list = []
	position = list_of_genes.index(locus_tag)
	for p in range(-frame, frame+1, 1):
		p_n = position + p
		genomic_neighbourhood_list.append(list_of_genes[p_n])

	return genomic_neighbourhood_list


#Definimos la función para que nos printee los resultados.
#El método escogido es por dataframe del módulo de pandas.
#En la primera columna se muestran las posiciones relativas
#de los locus vecinos, siendo el locus correspodiente a la posición 0
#el locus de interés (aquel con mejor resultado en el Blast).
#En la segunda columna se muestran los locus.En la tercera
#el producto de los mismos. En la cuarta, la cadena (+/-)


def print_result(genomic_neighbourhood_list,dict_of_genes,frame, line1):

	position_list = []
	product_list = []
	cadena_list = []
	for n in range(-frame, frame+1, 1):
		position_list.append(n)
	for p in genomic_neighbourhood_list:
		product = dict_of_genes[p][1]
		product_list.append(product)
	for c in genomic_neighbourhood_list:
		cadena = dict_of_genes[c][0]
		if cadena == 1:
			cadena_list.append("Forward")
		elif cadena == -1:
			cadena_list.append("Reverse")
	data = {'Position' : position_list, 'Locus': genomic_neighbourhood_list,'Product': product_list, 'Strand' : cadena_list}

	df = pd.DataFrame(data)
	line1 = line1.strip("\n") #Eliminamos el \n de la línea
	#Utilizamos los [] para evitar pasar valores escalares y tener que poner un index. Así, lo pasamos en forma de lista.
	hit = {'qseqid' : [line1.split("\t")[0]],  'sseqid' : [line1.split("\t")[1]],  'pident' : [line1.split("\t")[2]],  'qcovs' : [line1.split("\t")[3]]}
	best = pd.DataFrame(hit)
	print ("""
				***RESULTS***
""")
	print ("Best hit in local blast\n" + "\n" +  best.to_string(index=False) + "\n")

	print ("				Genomic context")
	print (df.to_string(index=False))


extractprotein(input_file)
locus, line1 = blast()
dict_of_genes,list_of_genes= make_list_and_dict(input_file)
genomic_neighbourhood_list = get_genomic_neighbourhood_list(locus,list_of_genes,frame)
print_result(genomic_neighbourhood_list,dict_of_genes,frame, line1)
eliminar()
