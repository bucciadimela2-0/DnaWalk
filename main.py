
### DNA walks: deepening our knowledge on wild-type genes and mutated genes ###
###### An example on lung cancer genes ######

# Main libraries

import tkinter as tk
import turtle
import xml.etree.ElementTree as ET

import cairosvg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from IPython.display import SVG, display
from PIL import Image
from skimage import measure
from svg_turtle import SvgTurtle


########Reading and parsing of the fasta gene file########
def read_fasta(file_path):

    c  = 0
    for seq_record in SeqIO.parse(file_path, "fasta"):
        print(f"ID: {seq_record.id}")
        print(f"Name: {seq_record.name}")
        print(f"Description: {seq_record.description}")
        print(f"Annotations: {seq_record.annotations}")
        print(f"Features: {seq_record.features}")
        print(f"Sequence Data: {repr(seq_record.seq)}")
        print(f"Length: {len(seq_record)}")
        print("\n")
        c+= 1  

    print(f"Quantity of gene in the chosen database: {c}")

    return seq_record

########Extracting interested exons ########
def extract_exons(file_path, exons, direction):

    finalsequence = " "
    c = 0
    
    for seq_record in SeqIO.parse(file_path, "fasta"):
        c =  c + 1

        if direction == 'fd' and  c<= exons:
            finalsequence = finalsequence + seq_record.seq
                    
        elif direction == 'bw' and c >= exons :
            finalsequence = finalsequence + seq_record.seq
    
    return finalsequence    


######Defining the DNA Walk ######

def dna_walk(sequence):
    x, y = 0, 0
    walk = [(x, y)]
    for nucleotide in sequence:
        if nucleotide == 'A':
            x -= 1
        elif nucleotide == 'T':
            x += 1
        elif nucleotide == 'C':
            y -= 1
        elif nucleotide == 'G':
            y += 1
        walk.append((x, y))
    return walk

##### Plotting the DNA Walk using Turtle ######

def plot_dna_walk(walk, name):
    t = SvgTurtle(300, 300)
    output_file = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/Images/'+name
    window = turtle.Screen()
    window.bgcolor("white")
    pen = turtle.Turtle()
    turtle.hideturtle
    pen.speed(1000)
    for x, y in walk:
        pen.goto(x*1.5, y*1.5)  
    turtle.getscreen().getcanvas().postscript(file = output_file+'.ps')
    print(f"Immagine della camminata del DNA salvata in {output_file}")
    convert_to_png(output_file)
    window.exitonclick()

##### Converting the output file in PNG file #####

def convert_to_png(path):
    img = Image.open(path+'.ps')
    img.save(path+'.png')
    print(f"Immagine convertita da .ps a .png")

####### Parsing the ALK sequence ########

file_path = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/dataset/ALK.fasta'

alk_seq = read_fasta(file_path)
print(f"Sequence Data: {alk_seq}")

#### Printing the ALK walk #####

alk_walk = dna_walk(alk_seq)
print(f"this is the DNA walk: {alk_walk} ")

output_file = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/Images/ALK/ALK.ps'

#### Plotting the ALK walk and converting to png ####

#plot_dna_walk(alk_walk, 'alk_walk')


#####Building the Mutation string #####
### reading the files, extracting the interested exons ###

exons_elm4_path = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/dataset/EML4Exons.fa'
read_fasta(exons_elm4_path)
exons_eml4_13 = extract_exons(exons_elm4_path, 13, 'fd')
print(exons_eml4_13)
#exons_eml4_walk = dna_walk(exons_eml4_13)
#plot_dna_walk(exons_eml4_walk)

exons_alk_path = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/dataset/ALKexons.fa'
read_fasta(exons_alk_path)
exons_alk_20 = extract_exons(exons_alk_path, 20, 'bw')
print(exons_alk_20)
#exons_alk_walk = dna_walk(exons_alk_20)
#plot_dna_walk(exons_alk_walk)

### Building and printing the mutated gene ###
alk_eml4_mutation = exons_eml4_13 + exons_alk_20
print(alk_eml4_mutation)

alk_eml4_walk = dna_walk(alk_eml4_mutation)
plot_dna_walk(alk_eml4_walk, 'alk_eml4_walk')

