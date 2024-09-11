#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
import argparse
import os
from io import BytesIO

from ete3 import Tree, TreeStyle, NodeStyle, faces

# Función para dibujar y guardar el árbol de ete3
def draw_with_ete(tree_file, output_file, label_size=10, highlight_seq=None, vertical_margin=10):
    try:
        # Cargar el árbol usando el formato adecuado
        tree = Tree(tree_file, format=1)  # Ajusta 'format' si es necesario

        # Crear un estilo personalizado para el árbol
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.show_branch_length = False
        ts.show_branch_support = False
        ts.mode = "r"  # Modo circular (puedes cambiarlo a "r" para radial)

        # Ajustar la separación vertical entre las ramas
        ts.branch_vertical_margin = vertical_margin  # Aumentar la separación vertical entre las ramas

        # Ajustar el estilo de cada hoja (nodo)
        for node in tree.iter_leaves():

            # Cambiar el color de la etiqueta si coincide con la secuencia a resaltar
            if node.name == highlight_seq:
                # Crear una etiqueta en rojo para la secuencia destacada
                name_face = faces.TextFace(node.name, fsize=label_size, fgcolor="red")
            else:
                # Etiquetas normales en azul oscuro para el resto de las secuencias
                name_face = faces.TextFace(node.name, fsize=label_size, fgcolor="darkblue")
                
            # Crear un estilo de texto personalizado para cada hoja
            #name_face = faces.TextFace(node.name, fsize=label_size, fgcolor="darkblue")
            node.add_face(name_face, column=0, position="branch-right")  # Colocar el texto a la derecha de la rama

            # Personalizar el estilo de cada nodo (puntos de los nodos)
            node_style = NodeStyle()
            node_style["size"] = 0  # Tamaño del nodo (ajustar según sea necesario)
            node.set_style(node_style)

        # Ajustar el estilo de las ramas
        for node in tree.traverse():
            if not node.is_leaf():
                # Ajustar el tamaño de los valores de branch length
                branch_length_face = faces.TextFace(f"{node.dist:.1e}", fsize=label_size - 1, fgcolor="black")
                node.add_face(branch_length_face, column=1, position="branch-top")  # Mostrar debajo de la rama

                # Ajustar el tamaño de los valores de branch support (soporte)
                if node.support:
                    branch_support_face = faces.TextFace(f"{node.support:.2f}", fsize=label_size - 1, fgcolor="red")
                    node.add_face(branch_support_face, column=1, position="branch-bottom")  # Mostrar encima de la rama


        # Guardar el árbol como imagen (por ejemplo, PNG)
        tree.render(output_file, w=800, units="px", tree_style=ts)
        return output_file
    except Exception as e:
        st.error(f"Error al cargar el archivo Newick: {e}")
        return None

# Función para capturar los argumentos desde la línea de comandos
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_name', type=str, required=True)
    parser.add_argument('--id_to_show', type=str, default="all")
    return parser.parse_args()

# Leer los argumentos
if __name__ == "__main__":
    # Necesitamos evitar que Streamlit se interponga con los argumentos, por lo que usamos sys.argv para pasarlos correctamente
    args = parse_args()

    # Managing directories
    mmseqs_workdir = args.job_name
    

    best_matches_file = f'{mmseqs_workdir}/best_matches.m8'
    
    # Cargar el archivo de resultados como un DataFrame
    df = pd.read_csv(best_matches_file, sep='\t', header=None, skiprows=1)
    df.columns = ['qseqid', 'tseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'organism', 'genomeid', 'proteinid', 'geneid', 'description']

    if not args.id_to_show == "all":
        output = df[df['qseqid'] == args.id_to_show]
        st.title(f"Report of {mmseqs_workdir}: {args.id_to_show}")
    else:
        output = df
        st.title(f"Report of {mmseqs_workdir}: All")
        
    st.header("Homology searching results")
    # Mostrar la tabla
    st.dataframe(output)

    if not args.id_to_show == "all":
        tree_file = f'{mmseqs_workdir}/trees/{args.id_to_show}_tree.nwk'    
        if os.path.exists(tree_file):
            seq_name = output.iloc[0,0]
            # Cargar y dibujar el árbol filogenético
            st.header("Phylogenetic tree")
            
            label_size = 6
            highlight_seq = seq_name
            branch_width = 2
            vertical_margin = 6
            output_file = f'{mmseqs_workdir}/trees/{args.id_to_show}_tree.png'    
            output_image = draw_with_ete(tree_file=tree_file, output_file=output_file, label_size=label_size, highlight_seq=highlight_seq, vertical_margin=vertical_margin)
    
            # Verificar si se generó la imagen correctamente
            if output_image and os.path.exists(output_image):
                # Mostrar la imagen en Streamlit
                image_bytes = BytesIO(open(output_image, 'rb').read())
                st.image(image_bytes, use_column_width=True)
