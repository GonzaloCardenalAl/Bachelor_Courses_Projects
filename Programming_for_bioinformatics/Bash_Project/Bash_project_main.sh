#!/bin/bash
#Ejercicio Evaluable bash v1.0
#Gonzalo Cardenal Antolín 2021

# This script allow to performs blastp analysis using a urls_file to download the subject proteomes
# this file have to contain two columns, first column contains the species (identifier) for that proteome
# and the url. The script will generate a output project folder which contains two subfolders, data (here
# downloaded proteomes are stored and input fasta files) and results (here blast result, in addition a folder is created
# for every query protein to store blast result, aligment and trees for each protein)

#Para que ejecute el script en inglés
export LC_ALL=C

#Function to print help in terminal
ayuda () {
   echo 'Ejercicio Evaluable bash v1.0'
   echo 'Gonzalo Cardenal Antolin 2021'
   echo -e '\nusage: ejercicio.sh <query_sequences.fa> <ncbi_urls_file> <output_folder> <blast-identity> <blast-coverage>\n'
   echo -e 'query_sequences.fa : a fasta/multifasta file containing protein query sequences for blast'
   echo -e 'ncbi_urls_file     : a text plain file containing species name and url to download fasta protein file'
   echo -e 'output_folder	   : folder in which data and results will be stored'
   echo -e 'blast-identity     : sequence identity cut off value 0-100'
   echo -e 'blast-coverage     : sequence coverage cut off value 0-100\n'
   exit
}

#Arguments assgination
query=$1
ncbi_urls_file=$2
project_name=$3
iden=$4 #70
cov=$5 #40

#Controls help message

case $1 in 
	-help | -h | HELP | -H | -Help ) ayuda;;
esac

#Control of arguments number 

if [[ $# -lt 5 ]]
then
	echo -e "error:too few aguments\n"
	ayuda 
	exit
fi

#Control del formato fasta del query y del rango de coverage e identity
	#De esta forma sabemos cuántos query en formato fasta contiene nuestro archivo a analizar
	#Aunque no lo pida el guión, se pararía la ejecución de programa si el query no está en formato fasta
nfasta=$(grep '^>' $query | wc -l )
if [ $nfasta -ne 0 ]
then
	echo -e  "\tPerfect! File $1 contains $nfasta query sequences in Fasta format\n"
else
echo -e "$1 does not contain fasta format sequences\n"
	exit
fi

if [ $4 -gt 100 ]
then
	echo "error: blast-identity value bigger than 100\n"
	ayuda
	exit
fi

if [ $5 -gt 100 ]
then
	echo "error: blast-coverage value bigger than 100\n"
	ayuda
	exit
fi

#Create project directories and output_files
	
	#Detecta si existe ya un directorio con el $project_name si es así, escribe la pregunta y permite escribir por medio de read
	#Answer que es la variable escrita por read entra en un case y en función de la respuesta ejecuta la acción correcta
if [ -d $project_name ]
then
	echo -e "Warning!!! Directory $project_name already exists, would you like to continue ?\t"
	echo -e "If you continue, current $project_name directory and all the files inside will be deleted (y/n)\n"
	read answer
	case $answer in
		#La respuesta "y/Y/yes/YES" continuará con el script y reescribirá el directorio
		y | Y | Yes | yes | YES) rm -r $project_name ;;
		n | N | No | no | NO) echo -e "Relauch the programe introducing a new project_name\n"; exit ;;
		*)echo -e "Incorrect answer. Please,relauch the programe"; exit ;;
	esac
fi
	
	#Creamos los subdirectorios data y results
mkdir -p $project_name/data/
mkdir -p $project_name/results/

#parsing $ncbi_urls_file

echo -e "$0 is Running...\n"
echo "Fasta files will be stored in /$project_name/data/ "

	#Descagarmos los archivos y los descomprimimos
	#La línea con el pipe introduce cada url (columna 2) en la variable $line y uno a uno los va descargando y descomprimiendo
awk '{print $2}' $ncbi_urls_file | while read line
do 
	wget --directory-prefix ./$project_name/data/ $line 2> ./$project_name/log
	gzip -d ./$project_name/data/"${line##*/}" #los ## eliminan los urls dejando solo el nombre del archivo
done
################
	
################
cp $ncbi_urls_file ./$project_name/data 
cp $query ./$project_name/data		
	
	#Creamos el archivo proteome.fa y Identifiers.tsv
cd ./$project_name/data/
for fastafile in *.faa #recorre todos los .faa del directorio
do
	cat $fastafile >> proteome.fa
	#coge la especie correspondiente al archivo y la almacena en la variable specie
	specie=$(grep $fastafile $ncbi_urls_file | awk '{print $1}')
	echo -e "$specie\t$project_name/data/$fastafile" >> Identifiers.tsv
done
rm $ncbi_urls_file


	#Generate a big file that contains for every species all the fasta headers, necesary for species assgination in blast hits
echo "Generating species_fasta_ID.tsv file.. "

	#Aquí cogemos todos los IDs de cada proteoma gracias a ^>. y $1, los almacenamos en IDs.
	#con un grep del nombre del archivo en Identifiers.tsv obtenemos a que especie se corresponde cada archivo 
	#finalmente copiamos cada ID en el archivo species_fasta_ID.tsv con su respectiva especie almacenada en $specie y seperados por un tab.
for i in *.faa
do
   awk '/^>./ { print $1 }' $i > IDs 
   specie=$(grep $i Identifiers.tsv | awk '{print $1}')
   cat IDs | while read ids; do echo -e "$specie\t$ids" >> species_fasta_ID.tsv ; done
done
rm IDs

#### BLAST SECTION ####

#Blast analysis and result filtering
echo "Running Blast analsys and result filtering... "
	
	#Volvemos al directorio de trabajo
cd .. 
cd ..
	
	#Corremos el blastp
blastp -query query.fa -subject ./$project_name/data/proteome.fa -outfmt "6 qseqid sseqid pident qcovs sseq" 2>./log 1>./$project_name/results/blast_result.txt
	
	#Filtramos con el porcentaje de coverage e identity -v nos permite introducir variables del shell dentro de awk
awk -v iden=$iden -v cov=$cov '{if ($3 >= iden && $4 >= cov) {print $0}}' ./$project_name/results/blast_result.txt > ./$project_name/results/blast_result_filtered.tsv

	#Avanzamos al directorio del proyecto
cd ./$project_name/
	
	#Generamos el archivo blast_result_final.tsv
	#Utilizamos sed para quitar el ">" de fasta y así coincida con los IDs de blast_result_filtered
	#Con awk generamos el archivo final, la línea de código funciona de la siguiente forma. NR==FNR hace que como solo coinciden el número de columnas totales con el número
	#de columnas locales del primer archivo solo genere el array asociativo a partir de species_ID. Este array guarda para cada ID su correspondiente especie. 
sed 's/>//g' ./data/species_fasta_ID.tsv > species_ID
awk  'BEGIN{OFS="\t"} NR==FNR{specie[$2]=$1; next} ($2 in specie){print $1,specie[$2],$2,$3,$4,$5}'  species_ID ./results/blast_result_filtered.tsv > ./results/blast_result_final.tsv
rm species_ID

	#Creamos el archivo uniq_query_list.txt con los resultados de blast_result_final
	#donde con awk buscamos las proteínas que han dado match y las almacena de forma
	#única en un array que printeamos en el archivo
awk -F' ' '{protein[$1];}END{for (i in protein)print i;}' ./results/blast_result_final.tsv > ./results/uniq_query_list.txt
	
	#Creamos las carpetas para cada proteína con su archivo $proteina.fa con sus hits
	#para ello en la línea de awk la proteína actual del bucle debe coincidir con la de
	#blast_result_final y de esta forma se printean los hits de su respectiva proteína
cat ./results/uniq_query_list.txt | while read line
do
	mkdir ./results/$line
	awk -v protein=$line -F' ' '$1==protein {print ">"$2"\t"$3"\n"$6}' ./results/blast_result_final.tsv > ./results/$line/$line.fa
done

### MUSCLE SECTION ####

#Make MUSCLE Phylogenetic trees for each query protein

echo "Making MUSCLE alaignment and phylogenetic trees..."

	#Corremos el muscle de cada proteína para obtener el alineamiento con proteínas homólogas
	#y su árbol filogenético
cat ./results/uniq_query_list.txt | while read line
do
muscle -in ./results/$line/$line.fa -out ./results/$line/$line.aln 2> log
muscle -maketree -in ./results/$line/$line.aln -out ./results/$line/$line.aln.nw -cluster neighborjoining 2> log
done


### TERMINAL PRINTING ####

# Final stats to show in terminal 

echo -e "\n*** DONE ***"
echo "Results are available at /$project_name/results/"
	
	#NR acumula el número de filas, al ser un .tsv este número
	# se corresponderá con el número de hits ya que en cada fila hay 6 columnas siendo
	#última la secuencia
echo "Total Blast hits $(awk 'END{print NR}' ./results/blast_result_final.tsv)"
echo "hits were found for query proteins:"
	
	#Contamos los hits específicos para cada proteina, counhitsSPI funciona como un contador
	#donde cada vez que encuentra una proteína suma 1 al contador de esa específicamente.
	#realmente como los arrays en AWK son asociativos, la proteina sería el key
	#y el número el value que se va sumando 1 cada vez que encuentra su key.
awk 'BEGIN { ORS=" " };{counthitsSPI[$1]++}; END { for(SPIprotein in counthitsSPI) print counthitsSPI[SPIprotein], SPIprotein }' ./results/blast_result_final.tsv ; echo -e "\n"



