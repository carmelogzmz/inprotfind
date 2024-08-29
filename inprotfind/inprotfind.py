
# external libraries
import pandas as pd
import matplotlib.pyplot as plt
from colorama import Fore, Back, Style, init
from Bio import SeqIO, Phylo
import requests
from tqdm import tqdm
from ete3 import Tree, TreeStyle

# interanal libraries
import os
import subprocess
import shutil
import time
import importlib.resources as pkg_resources
import argparse
import tarfile

# resetting colorama
init(autoreset=True)

#############
# FUNCTIONS #
#############

'''
# verification function
#######################

This group of functions checks if the needed external software (mmseqs2, mafft
and fasttree) is installed in the working environment. If not, it will ask the
user to install it.
'''

# It checks if mmseqs2 is installed
def verifying_mmseqs2():
    try:
        subprocess.run(["mmseqs"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print("Error verifying mmseqs2:", e)
    except FileNotFoundError:
        raise EnvironmentError("MMseqs2 is not installed or is not in the PATH. Please, install it before using the 'inprotfind' library")

# It checks if mafft is installed
def verifying_mafft():
    try:
        subprocess.run(["mafft", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print("Error verifying MAFFT:", e)
    except FileNotFoundError:
        raise EnvironmentError("MAFFT is not installed or is not in the PATH. Please, install it before using 'inprotfind.build_tree' function")
        
# It checks if fasttree is installed
def verifying_fasttree():
    try:
        subprocess.run(["FastTree", "-help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print("Error verifying FastTree:", e)
    except FileNotFoundError:
        raise EnvironmentError("FastTree is not installed or is not in the PATH. Please, install it before using 'inprotfind.build_tree' function")
   
'''
# example results function
##########################
'''

def show_example_result(example):
    example = int(example)
    if example > 0 and example <= 10:
        result = pkg_resources.files("inprotfind").joinpath("query_examples/protein_names.txt")
        df = pd.read_csv(result, sep=';', header=0)
        if example > 0 and example < 10:
            fila_deseada = df[df['example'] == f'query_example0{example}']
        elif example == 10:
            fila_deseada = df[df['example'] == f'query_example{example}']
        
        # Imprimir la cabecera y la fila deseada de forma tabulada
        if not fila_deseada.empty:
            print(pd.concat([df.head(0), fila_deseada]).to_string(index=False))
        else:
            print(Fore.RED + Style.BRIGHT + "Examples are numered from 1 to 10. Please, select a number in this range")
    else:
        print(Fore.RED + Style.BRIGHT + "Examples are numered from 1 to 10. Please, select a number in this range")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return   
    
    
'''
# get_database
##############

This function download and install the database arthropods_OrthoDB in the 
location of the library inprotfind in the current environment. In case an
installed database is corrupted or malfunction, it is possible to reinstall
it using this function
'''

def get_database(fm_calling = False):
        
    mmseqs_targetdir = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB")
    mmseqs_target_metadata = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB_metadata.parquet")

    if fm_calling == False:
        if not os.path.exists(mmseqs_targetdir) and not os.path.exists(mmseqs_target_metadata):
            install = True
            if not os.path.exists(pkg_resources.files("inprotfind").joinpath("databases")):
                os.mkdir(str(pkg_resources.files("inprotfind").joinpath("databases")))
        else:
            value = input(Fore.RED + Style.BRIGHT + f"A version of the database is already installed in your computer in {mmseqs_targetdir}. Do you want to reinstall it? (yes/no): ")
                
            if value == "yes" or value == "y" or value == "YES" or value == "Y":
                mmseqs_target_metadata = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB_metadata.parquet")
                try:
                    if os.path.exists(mmseqs_targetdir):
                        shutil.rmtree(mmseqs_targetdir)
                    if os.path.exists(mmseqs_target_metadata):
                        os.remove(mmseqs_target_metadata)
                    print(Fore.GREEN + Style.BRIGHT + "Previous database removed")
                    install = True
                except:
                    print("Error removing previous database")
                    print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
                    install = False
                    return
            else:
                print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
                return
    else:
        try:
            if os.path.exists(mmseqs_targetdir):
                shutil.rmtree(mmseqs_targetdir)
            if os.path.exists(mmseqs_target_metadata):
                os.remove(mmseqs_target_metadata)
            print(Fore.GREEN + Style.BRIGHT + "Previous database removed")
            install = True
        except:
            print("Error removing previous database")
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            install = False
            return
        
    if install == True:
        
        if not os.path.exists(pkg_resources.files("inprotfind").joinpath("databases")):
            os.mkdir(str(pkg_resources.files("inprotfind").joinpath("databases")))

        save_path = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB.tar.gz")
        extract_path = pkg_resources.files("inprotfind").joinpath("databases")
        doi = "https://zenodo.org/records/13386908/files/arthropodsDB.tar.gz?download=1"
        
        response = requests.get(doi, stream=True)
        total_size = int(response.headers.get('content-length', 0))  # Tamaño total del archivo
        block_size = 8192  # Tamaño del bloque a leer cada vez
    
        # Configurar tqdm para mostrar la barra de progreso
        #progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)
        print(Fore.GREEN + Style.BRIGHT + "Downloading database...")
        with tqdm(total=total_size, unit='iB', unit_scale=True, leave=False) as progress_bar:
            with open(save_path, 'wb') as file:
                for data in response.iter_content(block_size):
                    file.write(data)
                    progress_bar.update(len(data))

        progress_bar.close()
    
        if total_size != 0 and progress_bar.n != total_size:
            print(Fore.RED + Style.BRIGHT + "Error downloading the database")
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return
        else:
            print(Fore.GREEN + Style.BRIGHT + "The database was downloaded successfuly")
                
        if not os.path.exists(extract_path):
            os.makedirs(extract_path)  # Crear el directorio si no existe
    
        try:
            print(Fore.GREEN + Style.BRIGHT + "Installing database...")
            with tarfile.open(save_path, "r:gz") as tar:
                # Obtener el nombre de la carpeta raíz dentro del tar.gz
                root_dir = os.path.commonpath([member.name for member in tar.getmembers()])
                
                for member in tar.getmembers():
                    # Quitar la carpeta raíz del nombre del miembro
                    member.name = os.path.relpath(member.name, root_dir)
                    tar.extract(member, path=extract_path)
            
            # Eliminar el archivo tar.gz después de descomprimir
            os.remove(save_path)
            print(Fore.GREEN + Style.BRIGHT + f"Database have been installed in the inprotfind library (Path: {save_path})")

        except Exception as e:
            print(Fore.RED + Style.BRIGHT + f"Error extracting the files: {e}")
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return

###############################################################################
###############################################################################

'''
# find_matches
##############

This function queries the database (arthropods_OrthoDB) to search for sequences 
similar to the query sequence (query_path). Once found, it creates a file named
best_matches.m8 containing similarity statistics with the query sequence, and
it also adds metadata about the sequences in the table (in addition to the 
protein ID, it adds species, genome ID, gene ID, and protein description). All 
created files will be saved in a folder with a chosen name (job_name). If this 
function is executed and the database is not installed yet, it will download it
and install it.
'''

def find_matches(job_name, query_path):
        
    start_time = time.time()
    verifying_mmseqs2()
            
    db_name = "arthropods_OrthoDB"
    
    mmseqs_targetdir = pkg_resources.files("inprotfind").joinpath(f"databases/{db_name}")
    metadata_targetdir = pkg_resources.files("inprotfind").joinpath(f"databases/{db_name}_metadata.parquet")

    if not os.path.exists(mmseqs_targetdir):
        print("The default database (arthropods_OrthoDB) is not yet installed or is corrupted. The database will be downloaded and installed now.")
        get_database(True)
    elif not os.path.exists(metadata_targetdir):
        print("The metadata file for the default database (arthropods_OrthoDB) is not yet installed or is corrupted. The database will be downloaded and installed now.")
        get_database(True)
          
    # Crear carpetas temporales para MMseqs2
    mmseqs_workdir = job_name
    mmseqs_querydir = mmseqs_workdir + "/queryDB"
    mmseqs_resultdir = mmseqs_workdir + "/resultDB"
    mmseqs_tmp = mmseqs_workdir + "/tmp"
    
    if(query_path == "example1"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example01.fa")
    elif(query_path == "example2"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example02.fa")
    elif(query_path == "example3"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example03.fa")
    elif(query_path == "example4"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example04.fa")
    elif(query_path == "example5"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example05.fa")
    elif(query_path == "example6"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example06.fa")
    elif(query_path == "example7"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example07.fa")
    elif(query_path == "example8"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example08.fa")
    elif(query_path == "example9"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example09.fa")
    elif(query_path == "example10"):
        query_path = pkg_resources.files("inprotfind").joinpath("query_examples/query_example10.fa")
    
    # Verificar si el directorio existe y es un directorio
    if os.path.isdir(mmseqs_workdir):
        answer = input(Fore.RED + Style.BRIGHT + f"There is already a job folder named '{mmseqs_workdir}' in this directory. Do you want to replace it? ALL THE CURRENT FILES in the job folder will be ERASED (yes/no): ")
        if answer == "yes" or answer == "y":
            shutil.rmtree(mmseqs_workdir)
            os.mkdir(mmseqs_workdir)     
            print(Fore.GREEN + Style.BRIGHT + "Previous job erased.")
        else:
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return
    
    if not os.path.isdir(mmseqs_workdir):
        os.mkdir(mmseqs_workdir)
    if not os.path.isdir(mmseqs_querydir):
        os.mkdir(mmseqs_querydir)
    if not os.path.isdir(mmseqs_targetdir):
        os.mkdir(mmseqs_targetdir)
    if not os.path.isdir(mmseqs_resultdir):
        os.mkdir(mmseqs_resultdir)
    if not os.path.isdir(mmseqs_tmp):
        os.mkdir(mmseqs_tmp)
    
    print(Fore.GREEN + Style.BRIGHT + "Converting query to mmseqs2 format...")
    # Crear la base de datos MMseqs2
    subprocess.run(f"mmseqs createdb {query_path} {mmseqs_querydir}/queryDB", shell=True)

    # Ejecutar la búsqueda MMseqs2
    subprocess.run(f"mmseqs search {mmseqs_querydir}/queryDB {mmseqs_targetdir}/{db_name}DB {mmseqs_resultdir}/resultDB {mmseqs_tmp} --max-seqs 100", shell=True)
    end_time = time.time()

    
    print(Fore.GREEN + Style.BRIGHT + "Passing results to table and adding metadata...")
    # Convertir los resultados al formato tabular
    subprocess.run(f"mmseqs convertalis {mmseqs_querydir}/queryDB {mmseqs_targetdir}/{db_name}DB {mmseqs_resultdir}/resultDB {mmseqs_tmp}/best_matches_tmp.m8", shell=True)

    meta_df = pd.read_parquet(str(metadata_targetdir))
    best_matches_df_all = pd.read_csv(f"{mmseqs_tmp}/best_matches_tmp.m8", sep="\t", header=None)
    
    # Crear diccionarios para cada columna de interés
    code_to_organism = meta_df.set_index('PubProtID')['Organism'].to_dict()
    code_to_genomeid = meta_df.set_index('PubProtID')['GenomeID'].to_dict()
    code_to_pubgeneid = meta_df.set_index('PubProtID')['PubGeneID'].to_dict()
    code_to_description = meta_df.set_index('PubProtID')['Description'].to_dict()

    # Mapear cada una de las nuevas columnas a best_matches_df_all
    best_matches_df_all['Organism'] = best_matches_df_all[1].map(code_to_organism)
    best_matches_df_all['GenomeID'] = best_matches_df_all[1].map(code_to_genomeid)
    best_matches_df_all['PubGeneID'] = best_matches_df_all[1].map(code_to_pubgeneid)
    best_matches_df_all['Description'] = best_matches_df_all[1].map(code_to_description)
    best_matches_df_all['Organism'] = best_matches_df_all['Organism'].str.replace(' ', '_')
    best_matches_df_all['Description'] = best_matches_df_all['Description'].str.replace(' ', '_')
    
    header = ["qseqid", "tseqid","pident", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "organism", "genomeid", "geneid", "description"]
    best_matches_df_all.columns = header
    
    best_matches_df_all.to_csv(f"{mmseqs_workdir}/best_matches_all.m8", sep="\t", index=False, header=True)
    best_matches_df = best_matches_df_all[:30]
    best_matches_df.to_csv(f"{mmseqs_workdir}/best_matches.m8", sep="\t", index=False, header=True)
    with open(f"{mmseqs_workdir}/db_name.txt", "w") as file:
        file.write(db_name)
    
    # Limpiar archivos temporales de MMseqs2
    if os.path.exists(mmseqs_tmp):
        shutil.rmtree(mmseqs_tmp)
        print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} removed.")
    else:
        print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} does not exist or has already been removed.")
    
    end_time = time.time()     
    
    # Abriendo el archivo en modo de adición ('a' de append)
    with open("processing_time.txt", 'a') as archivo:
        # Añadiendo las líneas adicionales
        archivo.writelines(f"Searching for {job_name} complete in {end_time - start_time:.2f} seconds.\n")

    print(Back.GREEN + Fore.BLACK + f"Done! You can find all the matches in '{mmseqs_workdir}/best_matches_all.m8', and just the first 30 in '{mmseqs_workdir}/best_matches.m8'")
    print(Fore.GREEN + Style.BRIGHT + f"Searching for {job_name} complete in {end_time - start_time:.2f} seconds")
        
###############################################################################
###############################################################################

'''
# align_sequences
#################

Esta función utiliza el archivo best_matches.m8 conseguido con la función
anterior para conseguir las secuencias de aminoácidos de las proteínas
incluidas en esta tabla y añade la secuencia de aminoácidos de la secuencia
problema. Lo guarda en un archivo llamado filtered_sequences.fasta. Después,
realiza el alineamiento de las secuencias mediante Mafft y lo guarda en un 
archivo llamado aligned_sequences.fasta. La función trabaja con los datos que
encuentre en la carpeta del trabajo elegido (job_name). Requiere que la base
de datos con la que se creó el trabajo siga existiendo.
'''

def align_sequences(job_name):
        
    verifying_mmseqs2()
    verifying_mafft()
                
    # Crear carpetas temporales para MMseqs2
    mmseqs_workdir = job_name
    
    with open(f"{mmseqs_workdir}/db_name.txt", "r") as file:
        db_name = file.read()
    
    if db_name == None or db_name == "arthropods_OrthoDB":
        db_name = "arthropods_OrthoDB"
        mmseqs_targetdir = pkg_resources.files("inprotfind").joinpath(f"databases/{db_name}")
        metadata_targetdir = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB_metadata.parquet")
    
        if not os.path.exists(mmseqs_targetdir):
            raise FileNotFoundError(f"The default database (arthropods_OrthoDB) is missing from {mmseqs_targetdir}")
    
        if not os.path.exists(metadata_targetdir):
            raise FileNotFoundError(f"The default metadata file (arthropods_OrthoDB_metadata.parquet) is missing from {metadata_targetdir}")
    else:
        mmseqs_targetdir = "databases/" + db_name
        
    mmseqs_querydir = mmseqs_workdir + "/queryDB"
    mmseqs_filtereddir = mmseqs_workdir + "/filteredDB"
    mmseqs_tmp = mmseqs_workdir + "/tmp"
    
    if not os.path.exists(f"{mmseqs_targetdir}/{db_name}DB"):
        print(Fore.RED + Style.BRIGHT + f"The database '{db_name}' was used for this job, but now does not exist or it has been renamed/moved/removed/corrupted. It will be necessary to fix it or to create it again. Yo may use the 'create_targetDB' function to create it. Make sure you use the same name ({db_name}) and the same fasta file used the first time.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
        
    # Verificar si el directorio existe y es un directorio
    if not os.path.isdir(mmseqs_workdir):
        print(Fore.RED + Style.BRIGHT + f"The job folder named {mmseqs_workdir} does not exist. Please, choose an existing job folder or run first the 'find_matches' function.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
    
    if not os.path.isdir(mmseqs_filtereddir):
        os.mkdir(mmseqs_filtereddir)
    if not os.path.isdir(mmseqs_tmp):
        os.mkdir(mmseqs_tmp)
       
    print(Fore.GREEN + Style.BRIGHT + "Collecting sequences of the best matches...")
    # Carga del archivo best_matches.m8 (producido con la función "find_matches")
    result_file = f"{mmseqs_workdir}/best_matches.m8"
    # Leer identificadores de secuencias desde el archivo MMseqs2
    df = pd.read_csv(result_file, sep='\t', skiprows=1, header=None)
    sequence_ids = set(df[1].str.strip())
    
    # Guarda los ids de las secuencias en best_matches.m8
    sequence_ids = pd.DataFrame(sequence_ids)
    sequence_ids.to_csv(f"{mmseqs_tmp}/sequence_ids.txt", sep="\t", index=False, header=False)
    
    # Filtrar la base de datos utilizando la lista de sequence_ids_temp.txt
    subprocess.run(f"mmseqs createsubdb {mmseqs_tmp}/sequence_ids.txt {mmseqs_targetdir}/{db_name}DB {mmseqs_filtereddir}/filteredDB --id-mode 1", shell=True)
    
    # Convertir la sub-base de datos filtrada a formato FASTA
    subprocess.run(f"mmseqs convert2fasta {mmseqs_filtereddir}/filteredDB {mmseqs_tmp}/filtered_sequences_tmp.fasta", shell=True)
    
    # Convertir la sub-base de datos filtrada a formato FASTA
    subprocess.run(f"mmseqs convert2fasta {mmseqs_querydir}/queryDB {mmseqs_tmp}/query.fasta", shell=True)
    
    # Leer la secuencia problema
    query_seq = list(SeqIO.parse(f"{mmseqs_tmp}/query.fasta", "fasta"))[0]

    # Escribir la secuencia problema junto con las secuencias similares
    with open(f"{mmseqs_tmp}/filtered_sequences_tmp.fasta", "r") as infile, open(f"{mmseqs_workdir}/filtered_sequences.fasta", "w") as outfile:
        SeqIO.write(query_seq, outfile, "fasta")  # Escribir primero la secuencia problema
        outfile.write(infile.read())  # Luego escribir el resto de las secuencias
        
    print(Fore.GREEN + Style.BRIGHT + "Sequences aligning...")
    start_time = time.time()
    subprocess.run(f"mafft --auto {mmseqs_workdir}/filtered_sequences.fasta > {mmseqs_workdir}/aligned_sequences.fasta", shell=True)
    print(Fore.GREEN + Style.BRIGHT + f"Aligned sequences saved in {mmseqs_workdir}/aligned_sequences.fasta")
    end_time = time.time()
    print(Fore.GREEN + Style.BRIGHT + f"Alignment with MAFFT completed in {end_time - start_time:.2f} seconds.")
    
    # Limpiar archivos temporales de MMseqs2
    if os.path.exists(mmseqs_tmp):
        shutil.rmtree(mmseqs_tmp)
        print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} removed.")
    else:
        print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} does not exist or has already been removed.")
    
    print(Back.GREEN + Fore.BLACK + f"Done! You can find the filtered sequences in '{mmseqs_workdir}/filtered_sequences.fasta', and the alignment results in '{mmseqs_workdir}/aligned_sequences.fasta'")
    

###############################################################################
###############################################################################

'''
# build_tree
############

Esta función utiliza el archivo align_sequences.fasta para crear el árbol
filogenético donde situa la secuencia problema junto con las 30 secuancias más
similares. Crea un archivo llamado tree.nwk y dibuja un árbol. Se puede decidir
si dibujar el árbol de forma sencilla (por defecto) o la forma interactiva
con la librería "ete3"
'''

def build_tree(job_name, tree_type = 'simple'):
    
    verifying_fasttree()
    
    # Crear carpetas temporales para MMseqs2
    mmseqs_workdir = job_name
    
    # Verificar si el directorio existe y es un directorio
    if not os.path.isdir(mmseqs_workdir):
        print(Fore.RED + Style.BRIGHT + f"The job folder named {mmseqs_workdir} does not exist. Please, choose an existing job folder or run first the 'find_matches' function.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
        
    subprocess.run(f"FastTree -nt {mmseqs_workdir}/aligned_sequences.fasta > {mmseqs_workdir}/tree.nwk", shell=True)
    print(Fore.GREEN + Style.BRIGHT + f"Tree saved in {mmseqs_workdir}/tree.nwk")
    
    if(tree_type == 'simple'):
        print(Fore.GREEN + Style.BRIGHT + "Default style tree drawn.")
        tree = Phylo.read(f"{mmseqs_workdir}/tree.nwk", "newick")
        fig = plt.figure(figsize=(12,12))
        print(Back.GREEN + Fore.BLACK + "Done! You can find the tree drawn in the Plots tab or in a pop-up window")
        Phylo.draw(tree, do_show=True, show_confidence=True, axes=plt.gca())
        ax = fig.gca()
        ax.set_axis_off()
    
    if(tree_type == "interactive"):
        print(Fore.GREEN + Style.BRIGHT + "Interactive tree drawn.")
        # ARBOL CON ETE3        
        # Cargar el árbol
        t = Tree(f"{mmseqs_workdir}/tree.nwk")
        
        # Definir el estilo del árbol
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "c"
        ts.arc_start = -180  # 0 degrees = 3 o'clock
        ts.arc_span = 180
        
        print(Back.GREEN + Fore.BLACK + "Done! You can find the interactive tree drawn in a pop-up window")
        # Mostrar el árbol
        t.show(tree_style=ts)
    
    if(tree_type == "plt"):
        
        # Leer el árbol
        tree = Phylo.read(f"{mmseqs_workdir}/tree.nwk", "newick")
        
        # Crear una figura y un eje
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Dibujar el árbol en el eje
        Phylo.draw(tree, do_show=False, axes=ax)
        
        # Personalizar las ramas, etiquetas, etc.
        for clade in tree.find_clades():
            if clade.is_terminal():
                # Personalizar las etiquetas de los terminales
                label = clade.name
                ax.text(clade.branch_length, 0, label, va='center', ha='right', fontsize=12)
        
        # Eliminar los ejes
        ax.set_axis_off()
        
        print(Back.GREEN + Fore.BLACK + "Done! You can find the tree drawn in the Plots tab or in a pop-up window")
        # Mostrar la figura
        plt.show()
    
        
###############################################################################
###############################################################################

'''
# Code to run the functions from terminal
############

To run the functions, the command inprotfind before the function name is needed

For example, to run the find_matches function:
    
    inprotfind find_matches --job_name test01 --query_path queries/sequences01.fa
'''

def main_function():
    parser = argparse.ArgumentParser(description='Procesamiento de secuencias proteicas')
    
    # Subparsers para los comandos
    subparsers = parser.add_subparsers(dest='command', help='Comandos disponibles')
    
    parser_get_database = subparsers.add_parser('get_database', help="To download and install the target database")
    parser_get_database.add_argument("--fm_calling", type=bool, default=False, help="Controls if the function is call it from find_matches or not")
    # Subparser para find_matches
    parser_find_matches = subparsers.add_parser('find_matches', help='To find coincidences in the database')
    parser_find_matches.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for find_matches")
    parser_find_matches.add_argument("--query_path", type=str, required=True, help="Path to query file for find_matches")

    # Subparser para align_sequences
    parser_align_sequences = subparsers.add_parser('align_sequences', help='To align sequences')
    parser_align_sequences.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for align_sequences")

    # Subparser para build_tree
    parser_build_tree = subparsers.add_parser('build_tree', help='Construir un árbol filogenético')
    parser_build_tree.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for build_tree")
    parser_build_tree.add_argument("--tree_type", type=str, default=None, help="Tree type for build_tree. It may be 'simple' (default) or 'interactive'")

    args = parser.parse_args()

    if args.command == "get_database":
        get_database(args.fm_calling)
    elif args.command == "find_matches":
        find_matches(args.job_name, args.query_path)
    elif args.command == "align_sequences":
        align_sequences(args.job_name)
    elif args.command == "build_tree":
        build_tree(args.job_name, args.tree_type)
    else:
        print("The command was not recognised.")

if __name__ == "__main__":
    main_function()
    