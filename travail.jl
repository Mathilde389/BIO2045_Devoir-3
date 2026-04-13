# ---
# title: Titre du travail
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Hammel-Monzon
#      prenom: Valentina 
#      matricule: 20274033
#      github: premierAuteur
#    - nom: Nguyen
#      prenom: Mathilde
#      matricule: 20267325 
#      github: Mathilde389
#    - nom: Rochon
#      prenom: Marilou
#      matricule: XXXXXXXX 
#      github: DeuxiAut
# ---

# # Introduction
# Lors de propagation rapide de maladies infectieuses, comme les épidémies, il est essentiel
# de comprendre comment différentes stratégies d’intervention permettent de limiter la transmission 
# et la mortalité. En effet, cela représente un enjeu majeur en santé publique, surtout lorsque 
# les individus infectés peuvent transmettre la maladie sans présenter de symptômes. 
# Par exemple, dans le cas de la covid, le contrôle de la propagation de la maladie a été 
# compliqué du aux individus porteurs asymptomatiques (Zhang et al, 2021).

# Afin de contrôler les épidémies, il est important de mettre en place des outils comme
# le dépistage et la vaccination. Les tests de dépistages, comme les tests antigéniques
# rapides, permettent d’identifier rapidement les individus infectieux (Mina et al, 2020). 
# Les individus infectieux peuvent alors prendre des mesures de sécurités afin de pas 
# infecter le reste de la population. Les vaccins sont des mesures de prévention qui
# permettent de donner de l’immunité protectrice comme celui pour l’hépatite B (Schillie et al, 2018).
# Cela empêche le développement de la maladie lors d’un contact par la suite. 

# Cependant l’efficacité de ces stratégies dépend de plusieurs contraintes biologiques et
# logistiques. Les tests de dépistage ne sont pas totalement fiable et peuvent donner 
# des faux positives ou des faux négatifs. Les faux négatifs sont ceux qui mettent 
# à risque le reste de la population puisqu’un individu infecté est identifié comme 
# sain et peux donc infecter d’autres personnes. De plus les vaccins n’offrent pas 
#toujours une protection immédiate, et il y a donc une période pendant laquelle les
# individus vaccinés sont vulnérable. Le budget attribué aux stratégies de vaccination
# est également une limite.

# L’objectif de ce travail est donc de tester une stratégie de contrôle d’une maladie
# infectieuse dans une population en agissant sur le dépistage et la vaccination des 
# individus. 

# # Présentation de la stratégie d'intervention

# Dans le contrôle de maladies infectieuses, plusieurs stratégies peuvent être mise
# en place, notamment en ce qui concerne la vaccination et le dépistage.

# Pour la vaccination, il existe différentes stratégies, comme la vaccination de masse,
# la vaccination ciblée des groupes à risques ainsi que la vaccination en anneau.

# La vaccination de masse consiste à immuniser une grande partie de la population afin
# de réduire la transmission de la maladie (Bullen et al, 2023). Cette stratégie est
# très efficace, cependant, elle nécessite beaucoup de ressources comme une disponibilité
# suffisante de vaccins (Wouters et al, 2021). Cela rend son application compliquée
# surtout en début d’épidémie.  Cette stratégie a été utilisé lors de la COVID-19, et un
# des défis était l’accès aux ressources comme les vaccins et le budget (Wouters et al, 2021).

# La vaccination ciblée des groupes a risques, comme les personnes immunodéprimées ou
# les personnes âgées, permet de réduire la mortalité (Self et al, 2021). Cependant
# elle ne permet pas de casser la chaine de transmission dans la population. Cette stratégie
# a été également été utilisée lors de la COVID-19 (Self et al, 2021). En effet, au début de
# l’épidémie, les personnes à risques été vacciner en priorité.

# La vaccination en anneau consiste à vacciner les contacts proches des individus infectés
# (Henao-Restrepo et al, 2016). Cela permet de créer une barrière immunitaire. Cette stratégie
# permet de limiter localement la propagation de la maladie. Cette stratégie à été utilisé lors
# des épidémies d’Ébola et a été efficace pour contrôler la transmission de la maladie
# (Henao-Restrepo et al, 2016)

# Pour le dépistage, il existe également différentes stratégies, comme le dépistage massif
# et le dépistage ciblé.

# Le dépistage massif consiste à tester une grande partie de la population, ce qui permet de
# détecter un grand nombre de cas (Shen et al, 2021). Cependant, cette stratégie nécessite
# beaucoup de ressource et peut donc être difficile à garder en place sur le long terme
# (Shen et al, 2021). 

# Le dépistage ciblé consiste à tester en priorité certains individus, par exemple les contacts
# de cas confirmés (Kretzschmar et al, 2020). Cette stratégie utilise plus efficacement les
# ressources sur le long terme.

# Dans le modèle, une stratégie combinant un dépistage relativement massif ainsi qu’une
# vaccination en anneau à été choisie. En effet, le dépistage relativement massif permet
# d’identifier rapidement les individus infectieux. Une fois les cas détectés, la vaccination 
#des contacts permet de limiter la propagation locale, ce qui créer une barrière immunitaire.
# Cette stratégie a été choisie puisqu’elle permet de mieux utiliser les ressources. En effet,
# les tests de dépistages coutent moins cher, et la vaccination en anneau permet de vacciner
# moins de personnes. Cela permet donc une optimisation des ressources. La combinaison du
# dépistage et de la vaccination en anneau constitue donc une stratégie efficace et réaliste
# pour limiter la propagation d’une maladie infectieuse.





# # Présentation du modèle

# # Implémentation

# ## Packages nécessaires

import Random
Random.seed!(123456)
using CairoMakie

# Code du cours: 
using CairoMakie
CairoMakie.activate!(px_per_unit=6.0)
using StatsBase
import Random

# Puisque nous allons identifier des agents, nous utiliserons des UUIDs pour
# leur donner un indentifiant unique:

import UUIDs
UUIDs.uuid4()

# ## Création des types

# Le premier type que nous avons besoin de créer est un agent. Les agents se
# déplacent sur une lattice, et on doit donc suivre leur position. On doit
# savoir si ils sont infectieux, et dans ce cas, combien de jours il leur reste:

Base.@kwdef mutable struct Agent
    x::Int64 = 0
    y::Int64 = 0
    clock::Int64 = 21 # quand infectés, durée de vie est de 21 jours 
    infectious::Bool = false
    id::UUIDs.UUID = UUIDs.uuid4()
end

# On peut créer un agent pour vérifier:

Agent()

# La deuxième structure dont nous aurons besoin est un paysage, qui est défini
# par les coordonnées min/max sur les axes x et y:

Base.@kwdef mutable struct Landscape
    xmin::Int64 = -25
    xmax::Int64 = 25
    ymin::Int64 = -25
    ymax::Int64 = 25
end

# Nous allons maintenant créer un paysage de départ:

L = Landscape(xmin=-50, xmax=50, ymin=-50, ymax=50)

# ## Création de nouvelles fonctions

# On va commencer par générer une fonction pour créer des agents au hasard. Il
# existe une fonction pour faire ceci dans _Julia_: `rand`. Pour que notre code
# soit facile a comprendre, nous allons donc ajouter une méthode à cette
# fonction:

Random.rand(::Type{Agent}, L::Landscape) = Agent(x=rand(L.xmin:L.xmax), y=rand(L.ymin:L.ymax))
Random.rand(::Type{Agent}, L::Landscape, n::Int64) = [rand(Agent, L) for _ in 1:n]

# Cette fonction nous permet donc de générer un nouvel agent dans un paysage:

rand(Agent, L)

# Mais aussi de générer plusieurs agents:

rand(Agent, L, 3)

# On peut maintenant exprimer l'opération de déplacer un agent dans le paysage.
# Puisque la position de l'agent va changer, notre fonction se termine par `!`:

function move!(A::Agent, L::Landscape; torus=true)
    A.x += rand(-1:1)
    A.y += rand(-1:1)
    if torus
        A.y = A.y < L.ymin ? L.ymax : A.y
        A.x = A.x < L.xmin ? L.xmax : A.x
        A.y = A.y > L.ymax ? L.ymin : A.y
        A.x = A.x > L.xmax ? L.xmin : A.x
    else
        A.y = A.y < L.ymin ? L.ymin : A.y
        A.x = A.x < L.xmin ? L.xmin : A.x
        A.y = A.y > L.ymax ? L.ymax : A.y
        A.x = A.x > L.xmax ? L.xmax : A.x
    end
    return A
end

# Nous pouvons maintenant définir des fonctions qui vont nous permettre de nous
# simplifier la rédaction du code. Par exemple, on peut vérifier si un agent est
# infectieux:

isinfectious(agent::Agent) = agent.infectious

# Et on peut donc vérifier si un agent est sain:

ishealthy(agent::Agent) = !isinfectious(agent)

# On peut maintenant définir une fonction pour prendre uniquement les agents qui
# sont infectieux dans une population. Pour que ce soit clair, nous allons créer
# un _alias_, `Population`, qui voudra dire `Vector{Agent}`:

const Population = Vector{Agent}
infectious(pop::Population) = filter(isinfectious, pop)
healthy(pop::Population) = filter(ishealthy, pop)

# Nous allons enfin écrire une fonction pour trouver l'ensemble des agents d'une
# population qui sont dans la même cellule qu'un agent:

incell(target::Agent, pop::Population) = filter(ag -> (ag.x, ag.y) == (target.x, target.y), pop)

# ## Paramètres initiaux

# Notez qu'on peut réutiliser notre _alias_ pour écrire une fonction beaucoup plus
# expressive pour générer une population:

function Population(L::Landscape, n::Integer)
    return rand(Agent, L, n)
end

# On en profite pour simplifier l'affichage de cette population:

Base.show(io::IO, ::MIME"text/plain", p::Population) = print(io, "Une population avec $(length(p)) agents")

# Et on génère notre population initiale:

population = Population(L, 3750)

# Pour commencer la simulation, il faut identifier un cas index, que nous allons
# choisir au hasard dans la population:

rand(population).infectious = true

# Nous initialisons la simulation au temps 0, et nous allons la laisser se
# dérouler au plus 1000 pas de temps:

tick = 0
maxlength = 2000

# Pour étudier les résultats de la simulation, nous allons stocker la taille de
# populations à chaque pas de temps:

S = zeros(Int64, maxlength);
I = zeros(Int64, maxlength);

# Mais nous allons aussi stocker tous les évènements d'infection qui ont lieu
# pendant la simulation:

struct InfectionEvent
    time::Int64
    from::UUIDs.UUID
    to::UUIDs.UUID
    x::Int64
    y::Int64
end

events = InfectionEvent[]

# Notez qu'on a contraint notre vecteur `events` a ne contenir _que_ des valeurs
# du bon type, et que nos `InfectionEvent` sont immutables.

# ## Simulation

while (length(infectious(population)) != 0) & (tick < maxlength)

    ## On spécifie que nous utilisons les variables définies plus haut
    global tick, population

    tick += 1

    ## Movement
    for agent in population
        move!(agent, L; torus=false)
    end

    ## Infection
    for agent in Random.shuffle(infectious(population))
        neighbors = healthy(incell(agent, population))
        for neighbor in neighbors
            if rand() <= 0.4
                neighbor.infectious = true
                push!(events, InfectionEvent(tick, agent.id, neighbor.id, agent.x, agent.y))
            end
        end
    end

    ## Change in survival
    for agent in infectious(population)
        agent.clock -= 1
    end

    ## Remove agents that died
    population = filter(x -> x.clock > 0, population)

    ## Store population size
    S[tick] = length(healthy(population))
    I[tick] = length(infectious(population))

end

# ## Analyse des résultats

# ### Série temporelle

# Avant toute chose, nous allons couper les séries temporelles au moment de la
# dernière génération:

S = S[1:tick];
I = I[1:tick];

#-

f = Figure()
ax = Axis(f[1, 1]; xlabel="Génération", ylabel="Population")
stairs!(ax, 1:tick, S, label="Susceptibles", color=:black)
stairs!(ax, 1:tick, I, label="Infectieux", color=:red)
axislegend(ax)
current_figure()

# ### Nombre de cas par individu infectieux

# Nous allons ensuite observer la distribution du nombre de cas créés par chaque
# individus. Pour ceci, nous devons prendre le contenu de `events`, et vérifier
# combien de fois chaque individu est représenté dans le champ `from`:

infxn_by_uuid = countmap([event.from for event in events]);

# La commande `countmap` renvoie un dictionnaire, qui associe chaque UUID au
# nombre de fois ou il apparaît:

# Notez que ceci nous indique combien d'individus ont été infectieux au total:

length(infxn_by_uuid)

# Pour savoir combien de fois chaque nombre d'infections apparaît, il faut
# utiliser `countmap` une deuxième fois:

nb_inxfn = countmap(values(infxn_by_uuid))

# On peut maintenant visualiser ces données:

f = Figure()
ax = Axis(f[1, 1]; xlabel="Nombre d'infections", ylabel="Nombre d'agents")
scatterlines!(ax, [get(nb_inxfn, i, 0) for i in Base.OneTo(maximum(keys(nb_inxfn)))], color=:black)
f

# ### Hotspots

# Nous allons enfin nous intéresser à la propagation spatio-temporelle de
# l'épidémie. Pour ceci, nous allons extraire l'information sur le temps et la
# position de chaque infection:

t = [event.time for event in events];
pos = [(event.x, event.y) for event in events];

#

f = Figure()
ax = Axis(f[1, 1]; aspect=1, backgroundcolor=:grey97)
hm = scatter!(ax, pos, color=t, colormap=:navia, strokecolor=:black, strokewidth=1, colorrange=(0, tick), markersize=6)
Colorbar(f[1, 2], hm, label="Time of infection")
hidedecorations!(ax)
current_figure()

# # Modifications possibles

# Pendant le cours, formulez des hypothèses sur l'effet de 

# - la taille du paysage
# - la taille de la population
# - la dispersion sur une lattice toroïdale
# - la durée de l'épidémie
# - la survie de la population

# Étudiez le code en profondeur avant de commencer. Est-ce que certains
# paramètres sont représentés par des _magic numbers_ qui devraient être rendu
# explicites?

# Testez ces hypothèses en variant les paramètres du modèle. Est-ce qu'il existe
# des situations dans lesquelles la population est protégée contre l'épidémie?
# Des situations dans laquelle la structure spatiale de l'épidémie change?

# # Figures supplémentaires

# Visualisation des infections sur l'axe x

scatter(t, first.(pos), color=:black, alpha=0.5)

# et y

scatter(t, last.(pos), color=:black, alpha=0.5)


# ## Inclure du code

# Tous les fichiers dans le dossier `code` peuvent être ajoutés au travail
# final. C'est par exemple utile pour déclarer l'ensemble des fonctions du
# modèle hors du document principal.

# Le contenu des fichiers est inclus avec `include("code/nom_fichier.jl")`.

# Attention! Il faut que le code soit inclus au bon endroit (avant que les
# fonctions déclarées soient appellées).

include("code/01_test.jl")

# ## Une autre section

"""
    foo(x, y)

Cette fonction ne fait rien.
"""
function foo(x, y)
    ## Cette ligne est un commentaire
    return nothing
end

# # Présentation des résultats

# La figure suivante représente des valeurs aléatoires:

hist(randn(1000), color=:grey80)

# # Discussion
#
# Malgré les résultats obtenus, ce modèle présente plusieurs limites qui devraient 
# être prises en compte lors de l’interprétation des résultats. Tout d’abord, la 
# dynamique de transmission de la maladie utilisé dans cette simulation est très 
# simplifié. Par exemple, la probabilité d’infection est fixé à 0,4 pour tous les
#  individus et toutes les interactions. Cependant, dans les situations réelles la 
# transmission d’une maladie dépend de plusieurs facteurs, comme la charge virale, 
# la durée de contact et les caractéristiques individuelles (âge, état de santé etc)
# (Blaser et al, 2014). De plus, ces facteurs là dépendent également de la maladie transmisse.
# Par exemple, pour le VIH, la transmission dépend fortement de la charge virale 
#(Blaser et al, 2014), alors que pour la COVID-19, la transmission dépend plus de 
# la proximité ainsi que la durée de contact (Karia et al, 2020). Puisque le modèle
# présente la transmission de façon simpliste, cela limite le modèle à refléter
# fidèlement la complexité de la transmission de façon réelle. 
#
# Une autre limite du modèle est que les agents se déplace de manière aléatoire dans
# un environnement homogène. Cette modélisation ne prend pas en compte les interactions 
# sociales des individus. En effet, dans des situations réelles, les individus ont
# tendance à interagir avec certains groupes plus que d’autres. Par exemple, les enfants
# vont à l’école, les adultes au travail, les individus ont des familles etc (Mossong et al, 2008).
# Cela influence donc la transmission des maladies. Sans modélisation des interactions
# réelles entre les individus, la vitesse de propagation de la maladie peut être
# sous-estimée. Afin d’améliorer la simulation, il serait intéressant d’intégrer
# des réseaux sociaux en créant des zones de fortes densités. 

# De plus, dans ce modèle, la maladie est toujours fatale, ce qui est également une 
# simplification de la réalité. La plupart des maladies infectieuses ont des taux 
# de mortalité variable, avec un certain pourcentage d’individus qui guérissent et
# développent une immunité naturelle (Szomolanyi et al, 2010). Puisqu’il n’y a pas d’individus
# « rétablis », les dynamiques de rétablissement réelles ne sont pas représentés. 
# L’ajout de cet état permettrait de mieux représenter l’effet de la maladie sur la
# population. 
#
#En ce qui concerne les types d’intervention, des simplifications ont été faites pour
# la modélisation de la simulation. En effet, les tests de dépistages sont administrés
# de manière aléatoire et sans priorisation. Dans des scénarios réels, certains groupes
# de populations sont plus à risque que d’autres, par exemple les personnes âgées ou
# les personnes immuno-déficientes lors du COVID-19 (Liu et al, 2021) (Chapman et al, 2025).
# De plus, dans ce scénario, le vaccin est parfaitement efficace après deux générations
# pour tous les individus. Cependant, l’efficacité vaccinale peut être partielle et
# n’est pas efficace après le même temps d’incubation pour tout le monde (Zachreson et al, 2023).
# De plus, ce scénario ne prend pas en compte que le virus peut muter et qu’il est
# donc nécessaire de faire plusieurs doses de vaccins (Behl et al, 2022), comme lors de
# la covid. Afin d’améliorer le modèle, il serait intéressant d’introduire des 
#stratégies de dépistage et une probabilité de l’efficacité vaccinale.
#
# Par ailleurs, le déclanchement de l’intervention débute seulement après le premier
# décès. La réponse sanitaire est donc déclenchées assez tardivement. En effet, elle
# aurait pu être déclenchée dès le premier cas détecté. Par exemple, ((Kucharski et al, 2020).
# Il serait intéressant de faire plusieurs simulations avec différents temps de
# déclanchement de l’intervention afin de sélectionner celle qui induit le moins
# de décès.

# Finalement, le modèle est stochastique, ce qui est important pour représenter
# l’incertitude. Cependant, cela introduit une forte variabilité des résultats.
# Cette variabilité est prise en compte par la répétition des simulations mais
# une analyse plus approfondie permettrait de déterminer quels paramètres influencent
# le plus les résultats.



# On peut aussi citer des références dans le document `references.bib`, qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# `@` suivi de la clé de citation. Par exemple: @ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.
