=== Toolbox - Analyse Spectrale Paramétrique ===


Le fichier toolbox.m contient la fonction main. C'est ce fichier qui doit être lancé avec Matlab.


---- Paramètres ----
Les différents paramètres du signal d'entrée sont à rensigner au début du fichier toolbox.m dans les variables wn et an. Du bruit peut être ajouté en mettant la valeur de la variable noise à 1;


---- Choix de la méthode -----
Pour une estimation ponctuelle, le choix de la méthode se fait à l'aide de la variable method.


---- Graphiques de comparaison ----
La variable displayErrorGraph permet, si mise à 1, de calculer le graphe d'erreur pour les différentes méthodes en fonction du SNR. Il est calculé pour le signal d'entrée renseigné plus haut dans la toolbox.


---- Méthodes de résolution ----
Chaque méthode implémentée prend en argument deux paramètres pour la reconstruction des fréquences : le signal bruité ainsi que le nombre de pulsations recherchées.
Le calcul des amplitudes utilise ces pulsations et est analogue pour toutes les méthodes implémentées.
