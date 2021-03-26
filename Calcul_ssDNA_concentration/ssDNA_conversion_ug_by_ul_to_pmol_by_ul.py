#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:37:37 2020

@author: christophe
"""


# DNA Conversion ssDNA : ug/ul to pmol/ul

ssDNAmassiq = float(input("Enter the mass concentration of your sample (ug/ul) : \n"))

NumOligo = int(input("Enter the oligos length bases of the primary probes : \n"))

ssDNAConc = ssDNAmassiq * 1000000 / 330 / NumOligo
# pour obtenir l'arrondi à 2 chiffres après la virgule :
ssDNAConc = round(ssDNAConc,2)

print('-'*40, '\n')
print('formula = ')
print(' '*12, '  Cm(g/L)', '\n', ' C (mol/L) = ---------', '\n', ' '*12,
      'Mw (g/mol)', '\n')
print(' '*14, '  Cm(g/L)', '\n', ' C (mol/L) = --------------', '\n', ' '*12,
      '330*N (g/mol)', '\n')
print(' '*14, '  Cm(ug/uL)',' '*8, '1', '\n', ' C (pmol/uL) = ----------- x ',
      '-'*14, '\n', ' '*16,'330*N', ' '*6, '1e-6 (ug/pmol)', '\n')
print('N is the number of nucleotides \n')
print('330 g/mol is the average molecular weight of a nucleotide \n')
print('-'*40, '\n')
print (f'The concentration of your sample is {ssDNAConc} pmol/ul')
