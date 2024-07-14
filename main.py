
### DNA walks: deepening our knowledge on wild-type genes ###
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

#Downloading and parsing of the fasta gene file

########Just to understand Bio #########
file_path = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/dataset/ALK.fasta'
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

####### Parsing the ALK sequence ########

alk_seq = seq_record.seq
print(f"Sequence Data: {alk_seq}")


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


#### Printing the ALK walk #####
alk_walk = dna_walk(alk_seq)
print(f"this is the DNA walk: {alk_walk} ")


##### Plotting the DNA Walk using Turtle ######

def plot_dna_walk(walk):
    t = SvgTurtle(300, 300)
    output_file = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/Images/ALK/ALK.ps'
    window = turtle.Screen()
    window.bgcolor("white")
    pen = turtle.Turtle()
    pen.speed(1000)

    for x, y in walk:
        pen.goto(x*1.5, y*1.5)  
    
    turtle.getscreen().getcanvas().postscript(file = output_file)
    print(f"Immagine della camminata del DNA salvata in {output_file}")


    
    window.exitonclick()

##### Converting the output file in PNG file #####

def convert_to_png(path):
    img = Image.open(path)
    img.save('/Users/kepler/Desktop/StatisticaFisica/DnaWalk/Images/ALK/ALK.png')
    print(f"Immagine convertita da .ps a .png")

##### Plotting The ALK WALK ########

output_file = '/Users/kepler/Desktop/StatisticaFisica/DnaWalk/Images/ALK/ALK.ps'

plot_dna_walk(alk_walk)
convert_to_png(output_file)




######







