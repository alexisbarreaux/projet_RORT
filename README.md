# projet_RORT

## Utilisation du code
Pour utiliser notre code, dans un terminal à la racine du projet dans lequel julia est lancé:
- pour faire tourner la version exacte sur une instance de grille fournie, utiliser les commandes:
    - include("./models/modelPath.jl")
    - pathSolve("taxe_grille_2x3.txt")

    Les valeurs renvoyées sont dans l'ordre : est ce que la solution est optimale, temps de résolution, valeur, meilleure borne, gap.

- pour faire tourner l'heuristique sur une instance de grille fournie:
    - include("./heuristics/heuristic2.jl")
    - heurSolve("taxe_grille_2x3.txt")

    Les valeurs renvoyées sont de même que précédemment : est ce que la solution est optimale, temps de résolution, valeur, meilleure borne, gap.


Pour faire tourner sur une autre instance, il faudrait la placer dans le dossier data avec le même formalisme que nous avons employé et qui permet à julia de facilement charger les variables.