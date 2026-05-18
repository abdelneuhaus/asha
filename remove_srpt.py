import os
import shutil

def clean_sr_pt_folders(root_path):
    # On normalise le chemin pour éviter les problèmes de \
    root_path = os.path.abspath(root_path)
    
    print(f"Analyse du dossier : {root_path}")
    
    # topdown=False est indispensable pour ne pas perturber la boucle 
    # os.walk quand on supprime un élément de la liste
    for dirpath, dirnames, filenames in os.walk(root_path, topdown=False):
        for dirname in dirnames:
            if dirname == "SR.PT":
                # Construction du chemin complet
                full_path = os.path.join(dirpath, dirname)
                
                try:
                    # Suppression récursive (dossier + contenu)
                    shutil.rmtree(full_path)
                    print(f"SUCCÈS : {full_path} a été supprimé.")
                except Exception as e:
                    print(f"ERREUR : Impossible de supprimer {full_path}. Raison : {e}")

# Utilisation avec votre chemin :
path_to_clean = r"D:/ANALYSIS_PAPER/240417_W1_FPs"
clean_sr_pt_folders(path_to_clean)