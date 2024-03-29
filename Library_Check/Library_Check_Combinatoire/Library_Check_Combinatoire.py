#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:11:17 2021

@author: christophe


Analyse des séquences d'un custom array

Installation de Blast+ pour le blast des séquences des sondes primaires :
    
Aller sur : ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
Télécharger la dernière version de blast+
Mettre le fichier téléchargé dans le fichier désiré
ouvrir le terminal a partir du dossier contenant le fichier Blast+ a décompréser
lancer la ligne de commande : (remplacer la XXXXX par la version téléchargée)
tar zxvpf ncbi-blastXXXXXX-x64-linux.tar.gz

Il faut créer un répertoire qui accueillera la database (blastdb) et un fichier qui
acueillera les fichiers FASTA généré lors du scipt (blast_fasta 
$ mkdir /home/christophe/Documents/Informatique/ncbi-blast-2.10.1+/blastdb
$ mkdir /home/christophe/Documents/Informatique/ncbi-blast-2.10.1+/blast_fasta

Afin de pouvoir appeler les lignes de commandes (ex : blastn), il faut modifier le fichier bashrc :
Vérifier que vous êtes bien dans le répertoire $HOME (répertoire personnel)
exécuter : 
$ gedit .bashrc
insérer les lignes a la fin du fichier
export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin
export BLASTDB=$HOME/ncbi-blast-2.10.1+/blastdb/
Saugarder
Fermer et rouvrir le terminal.
Si vous taper blastn h, vous obtiendrez l'aide (tout fonctionne donc correctement)

Pour créer la database sur laquelle vous voulez blaster votre séquence :
Télécharger le fichier FASTA contenant la séquence du génome de Drosophila Meg:
Aller sur le FTP de Flybase :
ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/
Aller dans le dossier contenant la dernière version (exemple : dmel_r6.34_FB2020_03/)
Aller dans le dossier fasta/
Télécharger le fichier dmel-all-chromosome-r6.34.fasta
Mettre le fichier dans le répertoire $HOME/ncbi-blast-2.10.1+/blastdb
Afin de créer la database, exécuter la ligne de commande ci dessous dans le terminal:
$ makeblastdb -in $HOME/ncbi-blast-2.10.1+/blastdb/dmel-all-chromosome-r6.34.fasta -dbtype nucl -out $HOME/ncbi-blast-2.10.1+/blast/blastdb/droso_nt
Dans le dossier blastdb, des fichiers droso_nt.ndb, droso_nt.nhr, droso_nt.not.....apparaitront
Pour vérifier la database dans le terminal:
$ blastdbcmd -db droso_nt  -info
voici ce qui doit s'afficher :
-------------------------------------------------------------------
Database: /home/christophe/ncbi-blast-2.10.1+/blast/blastdb/dmel-all-chromosome-r6.34.fasta
	1,870 sequences; 143,726,002 total bases

Date: Jul 8, 2020  1:58 PM	Longest sequence: 32,079,331 bases

BLASTDB Version: 5

Volumes:
	/home/christophe/ncbi-blast-2.10.1+/blastdb/droso_nt
---------------------------------------------------------------------
"""

import os

#Folder ou se trouve la bibliothèque à tester :
rootFolder = '/home/christophe/Documents/Projets/Multiplexed/Librairies/Juin_2021_Librairie_test_et_combinatoire'
#Folder ou se trouve tous les fichiers pour l'utilisation du script :
rootFolderScript = '/home/christophe/Documents/Informatique/Python/Scripts/Python_Script_Tof_git/Library_Check'


libraryFile = '2_Full_sequence_Only_Final'
barcodeFile = 'seq_brcd_on_primary.csv'
pu_fwFile = 'seq_primer_univ_fw.csv'
pu_RevFile = 'seq_primer_univ_rev.csv'
merFile = 'seq_MER_on_primary.csv'


libraryFilePath = rootFolder + os.sep +libraryFile
primer_fPath = rootFolderScript + os.sep + pu_fwFile
primer_rPath = rootFolderScript + os.sep + pu_RevFile
barcode_filePath = rootFolderScript + os.sep + barcodeFile
mer_filePath = rootFolderScript + os.sep + merFile

###-------------------CREATION DU DOSSIER RESULTAT----------------------
resultFolder = os.path.expanduser('~/Python_Results')
Lib_Check_Folder = 'Library_Check'
pathResultFolder= resultFolder+ os.sep + Lib_Check_Folder
if not os.path.exists(resultFolder):
    os.mkdir(resultFolder)
if not os.path.exists(pathResultFolder):    
    os.mkdir(pathResultFolder)


#%% Calcul de la longueur des sondes primaires (ne doit pas excéder 10% entre
# la plus petite et la plus grande sonde primaire)

seq_probes_Raw = []
with open(libraryFilePath, 'r', encoding='utf-8') as probes_file:
    for line in probes_file :
        seq_probes_Raw.append((line.upper().replace('\n', '')))
    

# Elimination des espaces si présents dans les sondes primaires
seq_probes_full = []
for seq in seq_probes_Raw :
    seq_probes_full.append(seq.replace(' ', ''))
print(seq_probes_full[0:5])
print('-'*70)
print (f"Au total, il y a {len(seq_probes_Raw)} sondes primaires")
print('-'*70)

# Calcul de la longueur minimale et maximale
len_min = []
len_max = []
len_total = []
len_total = list(set([len(x) for x in seq_probes_full]))
len_total.sort()
len_min = len_total [0]
len_max = len_total [-1]

# Verification que la difference de taille entre la sequence de la sonde primaire
# la plus petite et la plus grande n'est pas suppérieur a 10%
pourcent = 100-(len_min*100/len_max)
if pourcent < 10 :
    print('-'*70)
    print ("Après vérification, la différence de taille de séquence entre vos \
sondes primaires n'excède pas 10% ")
    print(f"Cette différence est de {len_max-len_min} bases au maximum entre \
vos différentes séquences")
    print('-'*70)
elif pourcent > 10 :
    print('-'*70)
    print(f"ATTENTION : il y a une trop grande différence de longueur de séquence \
entre vos sondes primaires, différence = {pourcent}%")
    print('Cette différrence ne doit pas excéder 10% ')
    print(f"Cette différence pour votre bibliothèque est de {len_max-len_min} bases")
    print('-'*70)

#%% Cellule permettant de créer les dictionnaires des primers universels
# des barcodes et des MER (readout probes)

pu_fw = {}
pu_rev = {}
dict_barcode = {}
dict_MER = {}

# Def qui permet de mettre sous forme de dictionnaire les fichiers
# contenant les noms et les séquences des primers universels Fw et Rev

def Open_file_dict(pathw_file, file_name) :
    with open (pathw_file, 'r', encoding ='utf-8') as file :
        temp = []
        for item in file :
            temp = item.replace('\n', '').split(',')
            file_name[temp[0]] = temp[1]
        return

# Bloc qui permet de mettre en forme les fichiers contenant les barcodes et
# les primers universels

Open_file_dict(primer_fPath, pu_fw)

Open_file_dict(primer_rPath, pu_rev)

Open_file_dict(barcode_filePath, dict_barcode)

Open_file_dict(mer_filePath,dict_MER)

#Récupération des noms des primers fw et rev pour le trie des différentes 
# librairies à la fin en fonction des couples de primers utilisés
name_fw =[]
for k in pu_fw.keys():
    name_fw.append(k)

name_rev =[]
for k in pu_rev.keys():
    name_rev.append(k)

compil_name_pu = []
for fw in name_fw :
    for rev in name_rev :
        compil_name_pu.append(fw +':'+ rev)


# Bloc qui permet de vérifier que tous c'est bien passé avec les
# fichier primers universels Fw et Rev et avec les barcodes
# en les imprimant pour vérification

print('-'*40,'\n',
      'liste des primers forward :\n')
for key, value in pu_fw.items() :
    print (key, ' : ', value)
    
print('-'*40,'\n',
      'liste des primers reverse :\n')
for key, value in pu_rev.items() :
    print (key, ' : ', value)
    
print('-'*40,'\n',
      'liste des barcodes :\n')
i = 0
for key, value in dict_barcode.items() : 
    if i < 8 :
        i += 1
        print (key, ' : ',value)
print('-'*40,'\n',
      'liste des Readout (MER) :\n')
i = 0
for key, value in dict_MER.items() : 
    if i < 8 :
        i += 1
        print (key, ' : ',value)
#%% Recherche des différentes bibliothèques dans un même custom array en
# fonction des primers universels


#Fonction permettant de sous-trier un dictionnaire contenant des séquences 
#en fonction des sequences primers universal fw, rev
#Dictionnaire sous la forme :
#{('BBx.Fw', 'BBx.Rev', 'Bcd_01', Bcd_02'):['SeqA','SeqB','SeqC', ....]}
def Sorting (library, primer, library_sorted):
    for key_library, value_library in library.items():
        for seq_library in value_library :
            temp_Bcd=[]
            for key_primer, seq_primer in primer.items() :
                if seq_primer in seq_library :
                    temp_Bcd.append(key_primer)
            temp_Bcd = tuple(temp_Bcd)
            if key_library+temp_Bcd not in library_sorted.keys() :
                library_sorted[key_library+temp_Bcd] = []
                library_sorted[key_library+temp_Bcd].append(seq_library)
            else :
                library_sorted[key_library+temp_Bcd].append(seq_library)
    return           


# Mise en forme de la variable seq_probes_full pour utiliser la fonction Sorting
# qui prend en entrée un dictionnaire :
#dic_seq_probes_full= {'seq_probes_full': ['SeqA', SeqB,'SeqC',.....]}
dic_seq_probes_full=dict()
dic_seq_probes_full[('seq_probes_full',)]= seq_probes_full

# 1er tri des séquences des sondes primaires en fonction du PU.fw
# Création d'un dictionnaire sous la forme {'primer.Fw' : [seqA, seqB, seqC..]} 
libraries_F = {} 
Sorting(dic_seq_probes_full, pu_fw, libraries_F)
# On obtient un dictionnaire  
#libraries_F = {('seq_probes_full','BBx.Fw') : ['SeqA', 'SeqB', ....}
# On élimine dans la clé sous forme de tulpe 'seq_probes_full'
keys_list = list(libraries_F.keys())
for i in keys_list:
        libraries_F[(i[1],)] = libraries_F.pop(i)
# On obtient le dictionnaire libraries_F
# {'BBx.Fw53': ['SeqA', 'SeqB',...]}list(key)


# 2eme tri des séquences des sondes primaires en fonction du PU.rev
# Création d'un dictionnaire sous la forme {'primer.Fw:primer.rev ': [seqA, seqB, seqC..]}
libraries_F_R = {} 
Sorting(libraries_F, pu_rev, libraries_F_R)


nbr_lib = len(libraries_F_R)
print('-'*80)
print(f"Après analyse, {nbr_lib} librairie(s) a (ont) été détectée(s) dans le custom array.")
for i, k in enumerate(libraries_F_R.keys(),1) :
    print('n°',i,'=' ,k[0],'/',k[1] ,'-'*5, 'nbre de sondes totales :', len(libraries_F_R[k]))
print('-'*80)




#%%---------------------ATTENTION section pour Barcodes-----------------------

libraries_F_R_RT = {}      
Sorting(libraries_F_R, dict_barcode, libraries_F_R_RT)
for i, (k,v) in enumerate(libraries_F_R_RT.items(),1) :
    print(i, k, ':', len(v), 'sondes primaires')

#%%--------------------ATTENTION section pour MER probes----------------------
 
libraries_F_R_RT = {}      
Sorting(libraries_F_R, dict_MER, libraries_F_R_RT)
for i, (k,v) in enumerate(libraries_F_R_RT.items(),1) :
    print(i, k, ':', len(v), 'sondes primaires')
    
#%% Résumé de l'analyse des librairies dans un fichier csv
summaryFile = '1_Library_Summary.csv'
SummaryFilePath = pathResultFolder + os.sep + summaryFile
with open(SummaryFilePath, 'w') as file :
    file.write('N°'+','+'PU.fw'+','+'PU.rev'+','+'Barcode/MER'+','+'Nbr probes'+'\n')
    for i, (key,value) in enumerate(libraries_F_R_RT.items(),1) :
        list_key = [x for x in key]
        file.write(str(i)+','+list_key[0]+','+list_key[1]+','+' - '.join(list_key[2:])+','+str(len(value))+'\n')

print('-'*70)
print(f"Les fichiers résultats ont été saucegardé dans le dossier {resultFolder}")

#%% recupère la première et dernière sequence des sondes primaires de chaque 
# locus pour les mettre dans un dictionnaire (first_last_seq_probe)
first_last_seq_probe = {}
for key,value in libraries_F_R_RT.items() :
    first_last_seq_probe[key] = [value[0],value[-1]]
print (first_last_seq_probe)

#%% Convertion du dictionnaire contenant les séquences des premieres et 
# dernieres sondes primaires pour chaque locus dans un fichier FASTA pour
# pourvoir les blaster
def Fasta_convert(dictionnaire, name_file) :
    with open(name_file, 'w', encoding='utf8') as file :
# le module texwrap permet de découper le texte en un nbre de caractère précis
# pour que la sequence ne dépasse pas 70pb pour le format fasta
        import textwrap
        for key, value in dictionnaire.items():
            seq1 = textwrap.wrap(value[0],70)
            seq2 = textwrap.wrap(value[1],70)
            file.write('>' + key + '\n')
            for x1 in seq1 :
                file.write(x1 + '\n')
            file.write('>' + key + '\n')
            for x2 in seq2 :
                file.write(x2 + '\n')
                
Fasta_convert(first_last_seq_probe, 'first_last_seq_probe.fasta')