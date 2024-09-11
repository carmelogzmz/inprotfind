'''
#################
# inprotfind.py #
#################

- author: Carmelo Gómez-Martínez
- github: https://github.com/carmelogzmz
- email: carmelogzmz@gmail.com

This script manage all the functions available in the library. 
'''

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


##################
# MAIN FUNCTIONS #
##################

'''
# get_database
##############

This function download and install the database arthropods_OrthoDB in the 
location of the library inprotfind in the current environment. In case an
installed database is corrupted or malfunction, it is possible to reinstall
it using this function. The argument fm_calling is just for communication
between functions and does not need to be changed by the user.
'''

def get_database(fm_calling = False):
    
    # setting directories
    mmseqs_targetdir = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB")
    mmseqs_target_metadata = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB_metadata.parquet")

    # checking the status of the database in the environment
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
        
    # Installing the database
    if install == True:
        if not os.path.exists(pkg_resources.files("inprotfind").joinpath("databases")):
            os.mkdir(str(pkg_resources.files("inprotfind").joinpath("databases")))

        # finding the database in zenodo.org
        save_path = pkg_resources.files("inprotfind").joinpath("databases/arthropods_OrthoDB.tar.gz")
        extract_path = pkg_resources.files("inprotfind").joinpath("databases")
        doi = "https://zenodo.org/records/13622813/files/arthropodsDB.tar.gz?download=1"
        
        response = requests.get(doi, stream=True)
        total_size = int(response.headers.get('content-length', 0))  # Tamaño total del archivo
        block_size = 8192  # Tamaño del bloque a leer cada vez
    
        # downloading
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
             
        # installing
        if not os.path.exists(extract_path):
            os.makedirs(extract_path)  # Crear el directorio si no existe
        try:
            print(Fore.GREEN + Style.BRIGHT + "Installing database...")
            with tarfile.open(save_path, "r:gz") as tar:
                root_dir = os.path.commonpath([member.name for member in tar.getmembers()])
                
                for member in tar.getmembers():
                    member.name = os.path.relpath(member.name, root_dir)
                    tar.extract(member, path=extract_path)
            
            # removing downloaded file to save disk space
            os.remove(save_path)
            save_path = pkg_resources.files("inprotfind").joinpath("databases")
            print(Fore.GREEN + Style.BRIGHT + f"Database have been installed in the inprotfind library (Path: {save_path})")

        except Exception as e:
            print(Fore.RED + Style.BRIGHT + f"Error extracting the files: {e}")
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return


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

def find_matches(job_name, query_path, evalue = 0.0000000001, min_seq_id = 0.7):
        
    start_time = time.time()
    verifying_mmseqs2() # verifies if mmseqs2 is installed
    
    # managing the database
    db_name = "arthropods_OrthoDB"
    
    mmseqs_targetdir = pkg_resources.files("inprotfind").joinpath(f"databases/{db_name}")
    metadata_targetdir = pkg_resources.files("inprotfind").joinpath(f"databases/{db_name}_metadata.parquet")

    if not os.path.exists(mmseqs_targetdir):
        print("The default database (arthropods_OrthoDB) is not yet installed or is corrupted. The database will be downloaded and installed now.")
        get_database(True)
    elif not os.path.exists(metadata_targetdir):
        print("The metadata file for the default database (arthropods_OrthoDB) is not yet installed or is corrupted. The database will be downloaded and installed now.")
        get_database(True)
          
    # setting directory's names for file storage
    mmseqs_workdir = job_name
    mmseqs_querydir = mmseqs_workdir + "/queryDB"
    mmseqs_resultdir = mmseqs_workdir + "/resultDB"
    mmseqs_tmp = mmseqs_workdir + "/tmp"
    
    # checks if query_path contains any of the example's names and pass it the correct path
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
    
    # checks if the directory with the job_name exists or not to create it
    if os.path.isdir(mmseqs_workdir):
        answer = input(Fore.RED + Style.BRIGHT + f"There is already a job folder named '{mmseqs_workdir}' in this directory. Do you want to replace it? ALL THE CURRENT FILES in the job folder will be ERASED (yes/no): ")
        if answer == "yes" or answer == "y":
            shutil.rmtree(mmseqs_workdir)
            os.mkdir(mmseqs_workdir)     
            print(Fore.GREEN + Style.BRIGHT + "Previous job erased.")
        else:
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return
    
    # create the rest of directories
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
    
    # building the mmseqs2 query database
    print(Fore.GREEN + Style.BRIGHT + "Converting query to mmseqs2 format...")
    # Crear la base de datos MMseqs2
    subprocess.run(f"mmseqs createdb {query_path} {mmseqs_querydir}/queryDB", shell=True)

    # executing mmseqs2 search in the database
    subprocess.run(f"mmseqs search {mmseqs_querydir}/queryDB {mmseqs_targetdir}/{db_name}DB {mmseqs_resultdir}/resultDB {mmseqs_tmp} --max-seqs 100 -e {evalue} --min-seq-id {min_seq_id}", shell=True)
    
    print(Fore.GREEN + Style.BRIGHT + "Passing results to table and adding metadata...")
    # transforming results to tabular format
    subprocess.run(f"mmseqs convertalis {mmseqs_querydir}/queryDB {mmseqs_targetdir}/{db_name}DB {mmseqs_resultdir}/resultDB {mmseqs_tmp}/best_matches_tmp.m8", shell=True)

    if not os.path.exists(f"{mmseqs_tmp}/best_matches_tmp.m8"):
        no_matches = "None of the sequences in the database match with the query sequence"
        with open(f"{mmseqs_workdir}/no_matches.txt", 'w') as file:
            file.write(no_matches)
            print(Fore.RED + Style.BRIGHT + no_matches)
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return
    elif os.path.getsize(f"{mmseqs_tmp}/best_matches_tmp.m8") == 0:
        no_matches = "None of the sequences in the database match with the query sequence"
        with open(f"{mmseqs_workdir}/no_matches.txt", 'w') as file:
            file.write(no_matches)
            print(Fore.RED + Style.BRIGHT + no_matches)
            print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
            return
    else:
        # reading the metadata of the database and the results
        meta_df = pd.read_parquet(str(metadata_targetdir))
        best_matches_df_all = pd.read_csv(f"{mmseqs_tmp}/best_matches_tmp.m8", sep="\t", header=None)
        
        # adding the metadata to the result file
        code_to_organism = meta_df.set_index('ID')['Organism'].to_dict()
        code_to_genomeid = meta_df.set_index('ID')['GenomeID'].to_dict()
        code_to_pubprotid = meta_df.set_index('ID')['PubProtID'].to_dict()
        code_to_pubgeneid = meta_df.set_index('ID')['PubGeneID'].to_dict()
        code_to_description = meta_df.set_index('ID')['Description'].to_dict()

        best_matches_df_all['Organism'] = best_matches_df_all[1].map(code_to_organism)
        best_matches_df_all['GenomeID'] = best_matches_df_all[1].map(code_to_genomeid)
        best_matches_df_all['PubProtID'] = best_matches_df_all[1].map(code_to_pubprotid)
        best_matches_df_all['PubGeneID'] = best_matches_df_all[1].map(code_to_pubgeneid)
        best_matches_df_all['Description'] = best_matches_df_all[1].map(code_to_description)
        
        # adding header to result file
        header = ["qseqid", "tseqid","pident", "length", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "organism", "genomeid", "proteinid", "geneid", "description"]
        best_matches_df_all.columns = header
        
        # saving result file as best_matches_all.m8 and best_matches.m8
        best_matches_df_all.to_csv(f"{mmseqs_workdir}/best_matches_all.m8", sep="\t", index=False, header=True)
        best_matches_df = best_matches_df_all.groupby('qseqid').head(30).reset_index(drop=True)
        best_matches_df.to_csv(f"{mmseqs_workdir}/best_matches.m8", sep="\t", index=False, header=True)
        
        # saving database_name to a file (not longer necessary)
        with open(f"{mmseqs_workdir}/db_name.txt", "w") as file:
            file.write(db_name)
        
        # cleaning temporal files
        if os.path.exists(mmseqs_tmp):
            shutil.rmtree(mmseqs_tmp)
            print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} removed.")
        else:
            print(Fore.GREEN + Style.BRIGHT + f"Temporal folder {mmseqs_tmp} does not exist or has already been removed.")
        
        end_time = time.time()     
        
        # saving processing time in a file, just if want to have a log of it
        with open("processing_time.txt", 'a') as archivo:
            archivo.writelines(f"Searching for {job_name} complete in {end_time - start_time:.2f} seconds.\n")

        # DONE
        print(Back.GREEN + Fore.BLACK + f"Done! You can find all the matches in '{mmseqs_workdir}/best_matches_all.m8', and just the first 30 in '{mmseqs_workdir}/best_matches.m8'")
        print(Fore.GREEN + Style.BRIGHT + f"Searching for {job_name} complete in {end_time - start_time:.2f} seconds")


'''
# align_sequences
#################

This function uses the best_matches.m8 file obtained with the previous function
to retrieve the amino acid sequences of the proteins included in this table 
and adds the amino acid sequence of the query sequence. It saves it in a file 
called filtered_sequences.fasta. Then, it aligns the sequences using Mafft and
saves the result in a file called aligned_sequences.fasta. The function works
with the data found in the folder of the selected job (job_name). It requires 
that the database used to create the job still exists.
'''

def align_sequences(job_name, ids_to_align=None):
        
    verifying_mmseqs2()
    verifying_mafft()
                
    # managing directories
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
    
    # Verifying if the database exists
    if not os.path.exists(f"{mmseqs_targetdir}/{db_name}DB"):
        print(Fore.RED + Style.BRIGHT + f"The database '{db_name}' was used for this job, but now does not exist or it has been renamed/moved/removed/corrupted. It will be necessary to fix it or to create it again. Yo may use the 'create_targetDB' function to create it. Make sure you use the same name ({db_name}) and the same fasta file used the first time.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
        
    # Verifying if the directory exists
    if not os.path.isdir(mmseqs_workdir):
        print(Fore.RED + Style.BRIGHT + f"The job folder named {mmseqs_workdir} does not exist. Please, choose an existing job folder or run first the 'find_matches' function.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
    
    # creating directories
    if not os.path.isdir(mmseqs_filtereddir):
        os.mkdir(mmseqs_filtereddir)
    if not os.path.isdir(mmseqs_tmp):
        os.mkdir(mmseqs_tmp)
     
    # starting the aligning
    print(Fore.GREEN + Style.BRIGHT + "Collecting sequences of the best matches...")
    
    # reading best_matches.m8
    result_file = f"{mmseqs_workdir}/best_matches.m8"
    df = pd.read_csv(result_file, sep='\t', skiprows=1, header=None)
    
    # Crear una carpeta para almacenar los alineamientos
    os.makedirs(f"{mmseqs_workdir}/alignments", exist_ok=True)
    
    # Obtener el listado único de secuencias de consulta en queryDB
    query_sequences = df[0].unique()
    if ids_to_align is not None:
            query_sequences = [seq for seq in query_sequences if seq in ids_to_align]
    
    start_time = time.time()
    # Procesar cada secuencia de consulta de forma individual
    for query_id in query_sequences:
        # Comprobar si el archivo de alineamiento ya existe
        aligned_file = f"{mmseqs_workdir}/alignments/{query_id}_aligned.fasta"
        if os.path.exists(aligned_file):
            print(f"Alineamiento para {query_id} ya existe, saltando...")
            continue  # Saltar esta secuencia si ya está alineada
        
        print(f"Procesando secuencia de consulta: {query_id}")
        
        # Filtrar las secuencias homólogas para la secuencia de consulta actual
        homologs_df = df[df[0] == query_id]
        sequence_ids = set(homologs_df[1].str.strip())
        
        # Guardar los IDs de las secuencias homólogas en un archivo temporal
        sequence_ids_df = pd.DataFrame(sequence_ids)
        sequence_ids_df.to_csv(f"{mmseqs_tmp}/sequence_ids.txt", sep="\t", index=False, header=False)
        
        # Filtrar la base de datos targetDB usando los IDs de las secuencias homólogas
        subprocess.run(f"mmseqs createsubdb {mmseqs_tmp}/sequence_ids.txt {mmseqs_targetdir}/{db_name}DB {mmseqs_filtereddir}/filteredDB_{query_id} --id-mode 1", shell=True)
        
        # Convertir la base de datos filtrada a formato FASTA
        subprocess.run(f"mmseqs convert2fasta {mmseqs_filtereddir}/filteredDB_{query_id} {mmseqs_tmp}/filtered_sequences_tmp.fasta", shell=True)
        
        # Obtener los metadatos de las secuencias del archivo .m8 (pubprotid)
        metadatos = pd.Series(homologs_df[14].values, index=homologs_df[1]).to_dict()
        
        # Sustituir NaN por el identificador original
        for key, value in metadatos.items():
            if pd.isna(value):  # Si el valor es NaN
                metadatos[key] = key  # Reemplazar con el key (ID original)
            
        # Leer el archivo FASTA filtrado y reemplazar los IDs por los pubprotid
        with open(f'{mmseqs_tmp}/filtered_sequences_tmp.fasta', 'r') as input_fasta, open(f'{mmseqs_tmp}/filtered_sequences_pubprotid.fasta', 'w') as output_fasta:
            for record in SeqIO.parse(input_fasta, "fasta"):
                original_id = record.id
                if original_id in metadatos:
                    # Reemplazar el ID en la cabecera con el PubProtID
                    record.id = metadatos[original_id]
                    record.description = metadatos[original_id]  # Actualizar también la descripción
                SeqIO.write(record, output_fasta, "fasta")
        
        # Convertir la base de datos queryDB a FASTA
        subprocess.run(f"mmseqs convert2fasta {mmseqs_querydir}/queryDB {mmseqs_tmp}/query.fasta", shell=True)
        
        # Leer la secuencia de la consulta actual en queryDB
        query_seq = [seq for seq in SeqIO.parse(f"{mmseqs_tmp}/query.fasta", "fasta") if seq.id == query_id][0]
        
        # Agregar la secuencia de consulta a las secuencias homólogas filtradas
        with open(f"{mmseqs_tmp}/filtered_sequences_pubprotid.fasta", "r") as infile, open(f"{mmseqs_workdir}/alignments/{query_id}_sequences.fasta", "w") as outfile:
            SeqIO.write(query_seq, outfile, "fasta")  # Escribir primero la secuencia de consulta
            outfile.write(infile.read())  # Luego, escribir las secuencias homólogas
        
        # Alinear las secuencias con MAFFT
        print(f"Alineando secuencias para {query_id}...")
        subprocess.run(f"mafft --auto {mmseqs_workdir}/alignments/{query_id}_sequences.fasta > {mmseqs_workdir}/alignments/{query_id}_aligned.fasta", shell=True)
    
        print(f"Resultados guardados en {mmseqs_workdir}/alignments/{query_id}_aligned.fasta")
        
    end_time = time.time()
    print(f"Alineamientos completados para todas la secuencias de consulta en {end_time - start_time:.2f} segundos.")
    
'''
# build_tree
############

This function uses the align_sequences.fasta file to create a phylogenetic 
tree, placing the query sequence alongside the 30 most similar sequences. It 
creates a file called tree.nwk and draws a tree. You can choose to draw the 
tree in a simple way (default) or interactively using the "ete3" library.
'''

def build_tree(job_name, query_id, tree_type = 'simple'):
    
    verifying_fasttree()
    
    # Managing directories
    mmseqs_workdir = job_name
    qseqid = query_id
    mmseqs_qseqid = f"{qseqid}_aligned.fasta"
    
    # Crear una carpeta para almacenar los alineamientos
    os.makedirs(f"{mmseqs_workdir}/trees", exist_ok=True)
    
    # Verifying if the directory exists
    if not os.path.isdir(mmseqs_workdir):
        print(Fore.RED + Style.BRIGHT + f"The job folder named {mmseqs_workdir} does not exist. Please, choose an existing job folder or run first the 'find_matches' function.")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return
    
    # Creating tree.nwk with FastTree
    subprocess.run(f"FastTree {mmseqs_workdir}/alignments/{mmseqs_qseqid} > {mmseqs_workdir}/trees/{qseqid}_tree.nwk", shell=True)
    print(Fore.GREEN + Style.BRIGHT + f"Tree saved in {mmseqs_workdir}/trees/{qseqid}_tree.nwk")
    
    
    # Drawing simple tree
    if(tree_type == 'simple'):
        print(Fore.GREEN + Style.BRIGHT + "Default style tree drawn.")
        tree = Phylo.read(f"{mmseqs_workdir}/trees/{qseqid}_tree.nwk", "newick")
        fig = plt.figure(figsize=(12,12))
        print(Back.GREEN + Fore.BLACK + "Done! You can find the tree drawn in the Plots tab or in a pop-up window")
        Phylo.draw(tree, do_show=True, show_confidence=True, axes=plt.gca())
        ax = fig.gca()
        ax.set_axis_off()
        print(Back.GREEN + Fore.BLACK + "Done! You can find the tree drawn in the Plots tab or in a pop-up window")
    
    # Drawing interactive tree
    if(tree_type == "interactive"):
        print(Fore.GREEN + Style.BRIGHT + "Interactive tree drawn.")     
        # reading the tree
        t = Tree(f"{mmseqs_workdir}/trees/{qseqid}_tree.nwk")
        
        # defining tree style
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "r"
        ts.arc_start = -180  # 0 degrees = 3 o'clock
        ts.arc_span = 180
        
        print(Back.GREEN + Fore.BLACK + "Done! You can find the interactive tree drawn in a pop-up window")
        # plot tree
        t.show(tree_style=ts)


'''
# Plot results
########################

This function run the script ipf_report.py that build a report with streamlit and 
show it in the browser.
'''

def show_results(job_name, id_to_show="all"):

    # Ruta del archivo app.py que ejecuta Streamlit
    app_script = pkg_resources.files("inprotfind").joinpath("ipf_report.py")
    
    # Define los archivos basados en el job_name
     
    # Ejecutar el script de Streamlit mediante un subproceso
    subprocess.run(["streamlit", "run", app_script, "--", 
                    "--job_name", job_name, "--id_to_show", id_to_show])

############################
#  COMPLEMENTARY FUNCTIONS #
############################

'''
# verification functions
########################

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

This function allow to check the metadata of any of the query examples included
in the library. It just need the number of the example (1 to 10).
'''

def show_example_result(example):
    example = int(example)
    if example > 0 and example <= 10:
        result = pkg_resources.files("inprotfind").joinpath("query_examples/protein_names.txt")
        df = pd.read_csv(result, sep='\t', header=0)
        if example > 0 and example < 10:
            fila_deseada = df[df['Example'] == f'example0{example}']
        elif example == 10:
            fila_deseada = df[df['Example'] == f'example{example}']
        
        # print the header and the desired row
        if not fila_deseada.empty:
            print(pd.concat([df.head(0), fila_deseada]).to_string(index=False))
        else:
            print(Fore.RED + Style.BRIGHT + "Examples are numered from 1 to 10. Please, select a number in this range")
    else:
        print(Fore.RED + Style.BRIGHT + "Examples are numered from 1 to 10. Please, select a number in this range")
        print(Fore.GREEN + Style.BRIGHT + "Execution stopped. Returning to the prompt line.")
        return   
   

###########################################
# RUNNING THE FUNCTIONS FROM THE TERMINAL #
###########################################

'''
To run the functions, the command inprotfind before the function name is needed.
For example, to run the find_matches function:
    
  inprotfind find_matches --job_name test01 --query_path queries/sequences01.fa
'''

def main_function():
    parser = argparse.ArgumentParser(description='Procesamiento de secuencias proteicas')
    
    # Subparsers for commands
    subparsers = parser.add_subparsers(dest='command', help='Comandos disponibles')
    
    # Subparser for get_database
    parser_get_database = subparsers.add_parser('get_database', help="To download and install the target database")
    parser_get_database.add_argument("--fm_calling", type=bool, default=False, help="Controls if the function is call it from find_matches or not")
    
    # Subparser for find_matches
    parser_find_matches = subparsers.add_parser('find_matches', help='To find coincidences in the database')
    parser_find_matches.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for find_matches")
    parser_find_matches.add_argument("--query_path", type=str, required=True, help="Path to query file for find_matches")
    parser_find_matches.add_argument("--evalue", type=float, default=0.0000000001, help="e value treshold")
    parser_find_matches.add_argument("--min_seq_id", type=float, default=0.7, help="minimum sequence identity")
    
    # Subparser for align_sequences
    parser_align_sequences = subparsers.add_parser('align_sequences', help='To align sequences')
    parser_align_sequences.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for align_sequences")
    parser_align_sequences.add_argument("--ids_to_align", type=list, default=None, help="List of query ids to align")

    # Subparser for build_tree
    parser_build_tree = subparsers.add_parser('build_tree', help='To build the phylogenetic tree')
    parser_build_tree.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for build_tree")
    parser_build_tree.add_argument("--query_id", type=str, required=True, help="Name of the query to build its tree")
    parser_build_tree.add_argument("--tree_type", type=str, default=None, help="Tree type for build_tree. It may be 'simple' (default) or 'interactive'")

    # subparser for show_results
    parser_show_results = subparsers.add_parser('show_results', help='To show the results with Streamlit')
    parser_show_results.add_argument("--job_name", type=str, required=True, help="Name of the 'job' for show_results")
    parser_show_results.add_argument("--id_to_show", type=str, default="all", help="Name of the 'job' for show_results")
    
    # subparser for show_example_results
    parser_show_example_result = subparsers.add_parser('show_example_result', help='To show the results with Streamlit')
    parser_show_example_result.add_argument("--example", type=str, required=True, help="Name of the 'job' for show_results")
    args = parser.parse_args()

    if args.command == "get_database":
        get_database(args.fm_calling)
    elif args.command == "find_matches":
        find_matches(args.job_name, args.query_path, args.evalue, args.min_seq_id)
    elif args.command == "align_sequences":
        align_sequences(args.job_name, args.ids_to_align)
    elif args.command == "build_tree":
        build_tree(args.job_name, args.query_id, args.tree_type)
    elif args.command == "show_results":
        show_results(args.job_name, args.id_to_show)
    elif args.command == "show_example_result":
        show_example_result(args.example)
    else:
        print("The command was not recognised.")

if __name__ == "__main__":
    main_function()

'''
END OF THE SCRIPT
'''
