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

from pathlib import Path
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()

    parser.add_argument(
        '-i',
        '--input',
        required = True,
        help = 'input librairie file'
    )

    parser.add_argument(
        '-o',
        '--output',
        default = Path.cwd(),
        help = 'Result file destination folder.\nDEFAULT: current directory'
    )
    parser.add_argument(
        '-l',
        '--label',
        required = True,
        choices = ['bcd', 'mer'],
        help = 'Type of labelling : direct labelling = mer indirect labelling = bcd'
    )

    args = parser.parse_args()



    #Folder ou se trouve ce sript, avec avec les fichiers barcode, primer_F, primer_R.....
    script_folder = Path(__file__).absolute()
    resources_folder = script_folder.parents[1]

    result_folder = Path(args.output)


    library_file_path = Path(args.input)
    barcode_file_path = resources_folder.joinpath('seq_brcd_on_primary.csv')
    primer_fw_path = resources_folder.joinpath('seq_primer_univ_fw.csv')
    primer_rev_path = resources_folder.joinpath('seq_primer_univ_rev.csv')
    mer_file_path = resources_folder.joinpath('seq_MER_on_primary.csv')


    #### Calcul de la longueur des sondes primaires (ne doit pas excéder 10% entre
    # la plus petite et la plus grande sonde primaire)

    seq_probes_Raw = []
    with open(library_file_path, 'r', encoding='utf-8') as probes_file:
        for line in probes_file :
            seq_probes_Raw.append((line.upper().replace('\n', '')))
        print(seq_probes_Raw[0:5])
        

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

 
    # Def qui permet de mettre sous forme de dictionnaire les fichiers
    # contenant les noms et les séquences des primers universels Fw et Rev

    def create_dict(pathw_file) :
        with open (pathw_file, 'r', encoding ='utf-8') as file :
            temp = []
            dic = {}
            for item in file :
                temp = item.replace('\n', '').split(',')
                dic[temp[0]] = temp[1]
            return dic

    # Mise en forme des listes de primers, barcodes et autres... sous forme d'un dictionnaire

    pu_fw = create_dict(primer_fw_path)

    pu_rev = create_dict(primer_rev_path)

    dict_barcode =  create_dict(barcode_file_path)

    dict_MER =  create_dict(mer_file_path)

    # Récupération des noms des primers fw et rev pour le trie des différentes 
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


    # Vérification sous forme de sortie concole que tous c'est bien passé avec les
    # fichier primers universels Fw et Rev et avec les barcodes

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

    # Analyse des différentes bibliothèques dans un même custom array en
    # fonction des primers universels forward
        
    # Création d'un dictionnaire qui rassemble les sequences de sondes primaires
    # en fonction des primers universels forward
    temp_libraries = {}
    libraries = {}
    for key, value in pu_fw.items() :
        temp_libraries[key] = list()
        for seq_probe in seq_probes_full:
            if value in seq_probe :
                temp_libraries[key].append(seq_probe)
                
    # Permet d'éliminer tout ce qui est vide, en ajoutant seulement
    # les résultats positifs dans un nouveau dictionnaire : Libraries 
    for key, value in temp_libraries.items() :
        if len(value) != 0 :
            libraries[key] = value
    

                    
    # Parametrage d'une fonction permettant de sous-trier un dictionnaire contenant 
    # des séquences en fonction des primers Fw, en fonction des sequences primers 
    # universal rev ou barcodes
    def Sorting(library, primer, library_sorted) :
        for key_library, value_library in library.items():
            for seq_library in value_library :
                for key_primer, seq_primer in primer.items() :
                    if seq_primer in seq_library :
                        if (key_library + ':' + key_primer) not in library_sorted.keys() :
                            library_sorted[key_library + ':' + key_primer] = []
                            library_sorted[key_library + ':' + key_primer].append(seq_library)
                        else :
                            library_sorted[key_library + ':' + key_primer].append(seq_library)
        return
    # Trie du dictionnaire (primer.Fw : sequence de la librairie) en fonction des
    # primers Reverse
    libraries_F_R = {}    
    Sorting(libraries, pu_rev, libraries_F_R)
    nbr_lib = len(libraries_F_R)
    print('-'*80)
    print(f"Après analyse, {nbr_lib} librairie(s) ont été détectée(s) dans le custom array.")
    for i, k in enumerate(libraries_F_R.keys(),1) :
        print('n°',i,'=' ,k)
    print('-'*80)  

     
    if args.label == 'bcd':    #ATTENTION section pour barcodes
        #Trie du dictionnaire (primer.Fw : sequence de la librairie) en fonction des Barcodes
        libraries_F_R_B = {}      
        Sorting(libraries_F_R, dict_barcode, libraries_F_R_B)
        for i, (k,v) in enumerate(libraries_F_R_B.items(),1) :
            print(i, k, ':', len(v), 'sondes primaires')
    elif args.label == 'mer':   #ATTENTION section pour Readout (MER)
        # Trie du dictionnaire (primer.Fw : sequence de la librairie) en fonction des Readout (MER)
        libraries_F_R_B = {}      
        Sorting(libraries_F_R, dict_MER, libraries_F_R_B)
        for i, (k,v) in enumerate(libraries_F_R_B.items(),1) :
            print(i, k, ':', len(v), 'sondes primaires')
        # Trie et affiche les barcodes/MER pour chaque librairie
        details_lib = {}
        for k in libraries_F_R_B.keys() :
            for item in compil_name_pu :
                if k.startswith(item) :
                    if item not in details_lib.keys() :
                        details_lib[item] = []
                        details_lib[item].append(k.split(':')[2])
                    else :
                        details_lib[item].append(k.split(':')[2])
        print(details_lib)
                    


    # Résumé de l'analyse des librairies dans un fichier csv
    # Création d'un dictionnaire regroupant toutes les infos
    summary = dict()
    for k,v in libraries_F_R.items() :
        summary[k] = [['Nbr_probes_library',len(v)]]
    for k,v in details_lib.items() :
        summary[k].extend([v,['Nbr_Bcd/MER', len(v)]])


    summary_file = '1_Library_Summary.csv'
    summary_file_path = result_folder.joinpath(summary_file)
    with open(summary_file_path, 'w') as file :
        file.write('Summary :'+'\n')
        file.write('N°'+','+'PU.fw'+','+'PU.rev'+','+'Nbr_Barcode/MER'+','+'Total nbr probes'+','+'Barcode/MER'+'\n')
        for i, (key,value) in enumerate(summary.items(),1) :
            list_key = key.split(':')
            file.write(str(i)+','+list_key[0]+','+list_key[1]+','+str(value[2][1])+','+str(value[0][1])+','+" ".join(value[1])+'\n')
        file.write('\n'*3)
        file.write('Details'+'\n')
        file.write('N°'+','+'PU.fw'+','+'PU.rev'+','+'Barcode/MER'+','+'Nbr probes'+'\n')
        for i, (key,value) in enumerate(libraries_F_R_B.items(),1) :
            list_key = key.split(':')
            file.write(str(i)+','+list_key[0]+','+list_key[1]+','+list_key[2]+','+str(len(value))+' sondes primaires'+'\n')
    print('-'*40, '\n')
    print(f'The summary file of the library analysis has been saved here: {summary_file_path}')

# # recupère la première et dernière sequence des sondes primaires de chaque 
# # locus pour les mettre dans un dictionnaire (first_last_seq_probe)
# first_last_seq_probe = {}
# for key,value in libraries_F_R_B.items() :
#     first_last_seq_probe[key] = [value[0],value[-1]]
# print (first_last_seq_probe)

# # Convertion du dictionnaire contenant les séquences des premieres et 
# # dernieres sondes primaires pour chaque locus dans un fichier FASTA pour
# # pourvoir les blaster
# def Fasta_convert(dictionnaire, name_file) :
#     with open(name_file, 'w', encoding='utf8') as file :
# # le module texwrap permet de découper le texte en un nbre de caractère précis
# # pour que la sequence ne dépasse pas 70pb pour le format fasta
#         import textwrap
#         for key, value in dictionnaire.items():
#             seq1 = textwrap.wrap(value[0],70)
#             seq2 = textwrap.wrap(value[1],70)
#             file.write('>' + key + '\n')
#             for x1 in seq1 :
#                 file.write(x1 + '\n')
#             file.write('>' + key + '\n')
#             for x2 in seq2 :
#                 file.write(x2 + '\n')
                
# Fasta_convert(first_last_seq_probe, 'first_last_seq_probe.fasta')

if __name__ == "__main__":
    main()

