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
# ---

# # Introduction

# # Présentation du modèle
# Nous utilisons un modèle de simulation individu-centré pour représenter la propagation d’une maladie infectieuse dans une population. 
# Les modèles individu-centrés permettent de représenter explicitement les interactions entre individus et d’évaluer l’effet de différentes 
# stratégies de contrôle, notamment le dépistage et la vaccination (Adam & Arduin, 2023). Chaque individu est représenté par un agent 
# évoluant dans un espace discret bidimensionnel (lattice) défini par des coordonnées allant de -50 à 50 sur les axes x et y. La population 
# initiale est composée de 3750 agents distribués aléatoirement dans cet espace.

# Chaque agent possède plusieurs caractéristiques : une position spatiale, un état sanitaire (infectieux ou sain), un compteur de durée 
# d’infection (fixé à 21 jours), ainsi qu’un identifiant unique. Lorsqu’un agent devient infectieux, il le reste pendant toute la durée 
# de la maladie, qui est toujours fatale dans ce modèle. La population est initialement entièrement naïve, c’est-à-dire qu’aucun individu 
# ne possède d’immunité au départ.

# À chaque pas de temps (génération), les agents se déplacent aléatoirement dans la grille. La transmission de la maladie se produit lorsqu’un 
# agent infectieux partage la même cellule qu’un agent sain. Dans ce cas, la probabilité de transmission est fixée à 0,4. Les agents infectieux 
# voient leur durée d’infection diminuer à chaque génération, et sont retirés de la population lorsque cette durée atteint zéro, ce qui correspond 
# à leur décès.

# Afin de simuler une intervention sanitaire, un mécanisme combinant dépistage et vaccination a été intégré au modèle. L’intervention est
# déclenchée uniquement après la détection du premier décès dans la population. Une stratégie de dépistage actif est alors mise en place à l’aide 
# de tests antigéniques rapides (RAT), qui permettent d’identifier les individus infectieux avec une probabilité élevée (Lasser et al., 2021 ; Tran et al., 2023). 
# Lorsqu’un cas est détecté, une intervention ciblée est appliquée à ses contacts, ce qui reflète les approches modernes de contrôle en épidémiologie 
# (Lasser et al., 2021).
# Dans ce modèle, les tests RAT présentent une probabilité de détection de 95 %, pour un coût unitaire de 4 $, avec un budget total d’intervention 
# fixé à 21 000 $. Ils offrent donc un moyen rapide d’identifier les individus infectieux, bien que leur sensibilité imparfaite puisse en limiter 
# l’efficacité dans le contrôle global de l’épidémie (Lasser et al., 2021).

# Lorsqu’un individu est testé positif, les agents présents dans la même cellule spatiale (considérés comme ses contacts proches) peuvent être vaccinés, 
# à condition que le budget le permette. Le coût d’une vaccination est de 17$ par individu. Le vaccin est entièrement efficace, mais son effet n’est actif 
# qu’après un délai de deux générations. Une fois le vaccin actif, l’agent ne peut plus être infecté ni transmettre la maladie.

#La simulation se poursuit jusqu’à l’extinction de la maladie (absence d’individus infectieux) ou jusqu’à un maximum de 2000 générations. Les principales 
# variables suivies sont le nombre d’individus sains et infectieux au cours du temps, ainsi que le nombre total de décès, permettant d’évaluer l’efficacité 
# de la stratégie de vaccination en comparaison avec un scénario sans intervention.


# # Implémentation

# ## Packages nécessaires
using CairoMakie
using StatsBase
import Random
# METTRE AU BON ENDROIT (PROJECT.TOML)

# Initialisation
Random.seed!(123456)
CairoMakie.activate!(px_per_unit=6.0)

# Identifiants uniques pour chaque agent
import UUIDs
UUIDs.uuid4()

# ## Création des types

# Structure d'un agent (individu)
Base.@kwdef mutable struct Agent
    x::Int64 = 0 # position x
    y::Int64 = 0 # position y
    clock::Int64 = 21 # durée de la maladie est de 21 jours
    infectious::Bool = false # état infectieux
    ## Vaccination
    vaccinated::Bool = false        # vacciné ou non
    vaccine_delay::Int64 = 0        # délai avant efficacité

    id::UUIDs.UUID = UUIDs.uuid4() # identifiant unique
end

# Création d'un agent pour vérifier:

Agent()

# Structure du paysage (espace de simulation)
Base.@kwdef mutable struct Landscape
    xmin::Int64 = -50
    xmax::Int64 = 50
    ymin::Int64 = -50
    ymax::Int64 = 50
end

# Création du paysage

L = Landscape(xmin=-50, xmax=50, ymin=-50, ymax=50)

# ## Génération d'agents aléatoires

# Crée un agent à une position aléatoire
Random.rand(::Type{Agent}, L::Landscape) = Agent(x=rand(L.xmin:L.xmax), y=rand(L.ymin:L.ymax))
# Crée plusieurs agents
Random.rand(::Type{Agent}, L::Landscape, n::Int64) = [rand(Agent, L) for _ in 1:n]

# Cette fonction nous permet donc de générer un nouvel agent dans un paysage:

rand(Agent, L)

# Mais aussi de générer plusieurs agents:

rand(Agent, L, 3)

# Déplacement des agents

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

# ## Fonctions utiles

# Déterminer si l'agent est infectueux 
isinfectious(agent::Agent) = agent.infectious

# Déterminer si l'agent est sain
ishealthy(agent::Agent) = !isinfectious(agent)

# Type population = liste d'agents
const Population = Vector{Agent}

# Filtrage des agents
infectious(pop::Population) = filter(isinfectious, pop)
healthy(pop::Population) = filter(ishealthy, pop)

# Agents dans la même cellule (même position)
incell(target::Agent, pop::Population) = filter(ag -> (ag.x, ag.y) == (target.x, target.y), pop)

# ## Initialisation

# Génération de la population
function Population(L::Landscape, n::Integer)
    return rand(Agent, L, n)
end

# Simplification de l'affichage de cette population
Base.show(io::IO, ::MIME"text/plain", p::Population) = print(io, "Une population avec $(length(p)) agents")

# Génération de la population initiale
population = Population(L, 3750)

# Cas index (premier infecté)
rand(population).infectious = true

# ## Paramètres intervention
budget = 21000 # budget total
cost_vaccine = 17 # coût vaccinated
cost_test = 4 # coût test

intervention_started = false   # commence après premier décès
first_death = false
deaths = 0                     # compteur de morts

# Fonction de test RAT (détection infection)
function test_agent(agent::Agent)
    if agent.infectious
        return rand() <= 0.95
    else
        return false
    end
end

# Fonction de vaccination
function vaccinate!(agent::Agent)
    if !agent.vaccinated
        agent.vaccinated = true
        agent.vaccine_delay = 2   # actif après 2 jours
    end
end

# JUSTIFIER CE BLOQUE
function run_simulation(L, n, budget_total; with_intervention=true)

    population = Population(L, 3750)
    rand(population).infectious = true

    budget = 21000
    intervention_started = false
    first_death = false
    deaths = 0

    tick = 0
    maxlength = 2000

    while (length(infectious(population)) != 0) & (tick < maxlength)

        tick += 1

        ## Déplacement
        for agent in population
            move!(agent, L; torus=false)
        end

        ## Mise à jour délai vaccin
        for agent in population
            if agent.vaccine_delay > 0
                agent.vaccine_delay -= 1
            end
        end

        ## Infection
        for agent in Random.shuffle(infectious(population))
            neighbors = healthy(incell(agent, population))
            for neighbor in neighbors
                if rand() <= 0.4 && !(neighbor.vaccinated && neighbor.vaccine_delay == 0)
                    neighbor.infectious = true
                end
            end
        end

        ## Progression maladie
        for agent in infectious(population)
            agent.clock -= 1
        end

        ## Décès
        before = length(population)
        population = filter(x -> x.clock > 0, population)
        after = length(population)

        ## Déclenche intervention
        if !first_death && after < before
            intervention_started = true
            first_death = true
        end

        deaths += (before - after)

        ## Intervention (activable)
        if with_intervention && intervention_started && budget > 0
            for agent in population
                if budget >= cost_test
                    budget -= cost_test

                    if test_agent(agent)
                        neighbors = incell(agent, population)
                        for n in neighbors
                            if budget >= cost_vaccine && !n.vaccinated
                                vaccinate!(n)
                                budget -= cost_vaccine
                            end
                        end
                    end
                end

                if budget <= 0
                    break
                end
            end
        end
    end

    return deaths
end

# ## Simulation

tick = 0
maxlength = 2000

# Stockage des résultats
S = zeros(Int64, maxlength); # sains
I = zeros(Int64, maxlength); # infectieux
D = zeros(Int64, maxlength)  # morts cumulés

# Événements d'infection
struct InfectionEvent
    time::Int64
    from::UUIDs.UUID
    to::UUIDs.UUID
    x::Int64
    y::Int64
end

events = InfectionEvent[]

run_simulation()

# ## Analyse des résultats
using Statistics

n_runs = 20

deaths_with = [run_simulation(with_intervention=true) for _ in 1:n_runs]
deaths_without = [run_simulation(with_intervention=false) for _ in 1:n_runs]

mean_with = mean(deaths_with)
mean_without = mean(deaths_without)

var_with = var(deaths_with)
var_without = var(deaths_without)

println("\n RÉSULTATS SUR ", n_runs, " SIMULATIONS")

println("\nAVEC intervention :")
println("Moyenne des morts = ", mean_with)
println("Variance = ", var_with)

println("\nSANS intervention :")
println("Moyenne des morts = ", mean_without)
println("Variance = ", var_without)

println("\nGAIN (morts évités) = ", mean_without - mean_with)

f2 = Figure()

ax = Axis(f2[1, 1];
    xlabel="Simulation",
    ylabel="Nombre de morts",
    title="Comparaison des décès avec et sans intervention"
)

scatter!(ax, 1:n_runs, deaths_with, label="Avec intervention", color=:blue)
scatter!(ax, 1:n_runs, deaths_without, label="Sans intervention", color=:red)

axislegend(ax)

current_figure()

println("Moyenne morts avec intervention = ", mean(deaths_with))
println("Moyenne morts sans intervention = ", mean(deaths_without))

# ### Série temporelle

# Avant toute chose, nous allons couper les séries temporelles au moment de la
# dernière génération:

S = S[1:tick];
I = I[1:tick];

# ## Analyse des résultats 

f = Figure()
ax = Axis(f[1, 1]; xlabel="Génération", ylabel="Population")
stairs!(ax, 1:tick, S, label="Sains", color=:black)
stairs!(ax, 1:tick, I, label="Infectieux", color=:red)
stairs!(ax, 1:tick, D, label="Morts", color=:gray)
axislegend(ax)
current_figure()

# Ajout: résultats finaux

println("Nombre total de morts : ", deaths)
println("Budget restant : ", budget)

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

# On peut aussi citer des références dans le document `references.bib`, qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# `@` suivi de la clé de citation. Par exemple: @ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.
