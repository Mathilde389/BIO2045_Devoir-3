# ---
# title: Titre du travail
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Hammel-Monzon
#      prenom: Valentina 
#      matricule: 20274033
#      github: valentina9000
#    - nom: Nguyen
#      prenom: Mathilde
#      matricule: 20267325 
#      github: Mathilde389

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
    tested::Bool = false 
      # Vaccination
    vaccinated::Bool = false        # vacciné ou non
    days_after_vax::Int64 = 0        # délai avant efficacité
    id::UUIDs.UUID = UUIDs.uuid4() # identifiant unique
end

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
function make_population(L::Landscape, n::Int)
    return rand(Agent, L, n)
end

# Simplification de l'affichage de cette population
Base.show(io::IO, ::MIME"text/plain", p::Population) = print(io, "Une population avec $(length(p)) agents")

# Génération de la population initiale
population = make_population(L, 3750)

# Cas index (premier infecté)
rand(population).infectious = true

# ## Paramètres intervention
budget = 21000 # budget total
cost_vaccine = 17 # coût vaccinated
cost_test = 4 # coût test

intervention_started = false   # commence après premier décès
first_death = false
deaths = 0                     # compteur de morts

# Fonction de vaccination
function vaccinate!(agent::Agent)
    if !agent.vaccinated
        agent.vaccinated = true
        agent.days_after_vax = 0   # délai avant activation = 2 jours (immunité après 2 ticks)
    end
end

# ## Simulation principale
Base.@kwdef struct InfectionEvent
    tick::Int
    from::UUIDs.UUID
    to::UUIDs.UUID
    x::Int
    y::Int
end

function run_simulation(L::Landscape, n::Int, budget_total; with_intervention=true)

    population = make_population(L, n)
    rand(population).infectious = true

    events = InfectionEvent[]  

    budget = budget_total
    intervention_started = false
    first_death = false
    deaths = 0

    tick = 0
    maxlength = 2000
    S = Int[]
    I = Int[]
    D = Int[]

    while (length(infectious(population)) > 0) && (tick < maxlength)

        tick += 1

        # Déplacement
        for agent in population
            move!(agent, L)
        end

        # Vaccination delay
        for agent in population
            if agent.vaccinated
                agent.days_after_vax += 1
            end
        end

        # Infection
        for agent in Random.shuffle(infectious(population))
            for neighbor in healthy(incell(agent, population))

                immune = neighbor.vaccinated && neighbor.days_after_vax >= 2

                if rand() <= 0.4 && !immune
                    neighbor.infectious = true
                    push!(events, InfectionEvent(tick, agent.id, neighbor.id, neighbor.x, neighbor.y))
                end
            end
        end

        # Disease progression
        for agent in infectious(population)
            agent.clock -= 1
        end

        # Death
        before = length(population)
        population = filter(a -> a.clock > 0, population)
        after = length(population)

        if !first_death && after < before
            intervention_started = true
            first_death = true
        end

        deaths += (before - after)

        # Intervention
        if with_intervention && intervention_started && budget > 0

            for agent in population

                if !agent.tested && budget >= 4
                    budget -= 4
                    agent.tested = true

                    positive = agent.infectious ? rand() <= 0.95 : rand() <= 0.05

                    if positive
                        for neighbor in incell(agent, population)
                            if budget >= 17 && !neighbor.vaccinated
                                vaccinate!(neighbor)
                                budget -= 17
                            end
                        end
                    end
                end

                if budget <= 0
                    break
                end
            end
        end
        # Suivi des états
        push!(S, length(healthy(population)))
        push!(I, length(infectious(population)))
        push!(D, deaths)
    end

    return deaths, budget, events, S, I, D
end

# ## Analyse des résultats
using Statistics

n_runs = 50

results_with = [run_simulation(L, 3750, 21000; with_intervention=true) for _ in 1:n_runs]
results_without = [run_simulation(L, 3750, 21000; with_intervention=false) for _ in 1:n_runs]

# Une simulation pour visualisation (courbes temporelles)
sim_with = run_simulation(L, 3750, 21000; with_intervention=true)
sim_without = run_simulation(L, 3750, 21000; with_intervention=false)

S_with, I_with, D_with = sim_with[4], sim_with[5], sim_with[6]
S_without, I_without, D_without = sim_without[4], sim_without[5], sim_without[6]

deaths_with = [r[1] for r in results_with]
budget_with = [r[2] for r in results_with]

deaths_without = [r[1] for r in results_without]
budget_without = [r[2] for r in results_without]

println("Budget moyen restant (avec intervention) = ", mean(budget_with))

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

# Graphique comparatif
f2 = Figure()

ax = Axis(f2[1, 1];
    xlabel="Simulation",
    ylabel="Nombre de morts",
    title="Comparaison des décès avec et sans intervention"
)

scatter!(ax, 1:length(deaths_with), deaths_with, label="Avec intervention", color=:blue)
scatter!(ax, 1:length(deaths_without), deaths_without, label="Sans intervention", color=:red)

axislegend(ax)
current_figure()

# **Figure 1: Comparaison des décès**
# Cette figure présente, pour chacune des 50 simulations, le nombre total de décès observés 
# avec et sans intervention. Chaque point correspond à une simulation indépendante. Cette 
# figure permet de visualiser la variabilité des résultats entre les répétitions et de comparer 
# directement l’effet global de l’intervention sur la mortalité. Elle met en évidence la dispersion 
# des résultats ainsi que la tendance générale à une réduction du nombre de décès en présence 
# d’une intervention.

println("Moyenne morts avec intervention = ", mean(deaths_with))
println("Moyenne morts sans intervention = ", mean(deaths_without))

f3 = Figure()

# Graphiqe avec intervention
ax1 = Axis(f3[1, 1],
    title = "Évolution AVEC intervention",
    xlabel = "Temps",
    ylabel = "Population"
)

lines!(ax1, 1:length(S_with), S_with, label="Sains")
lines!(ax1, 1:length(I_with), I_with, label="Infectieux")
lines!(ax1, 1:length(D_with), D_with, label="Décédés")

axislegend(ax1)

# Graphique sans intervention
ax2 = Axis(f3[2, 1],
    title = "Évolution SANS intervention",
    xlabel = "Temps",
    ylabel = "Population"
)

lines!(ax2, 1:length(S_without), S_without, label="Sains")
lines!(ax2, 1:length(I_without), I_without, label="Infectieux")
lines!(ax2, 1:length(D_without), D_without, label="Décédés")

axislegend(ax2)

current_figure()

# **Figure 2: Dynamique  temporelle de l’épidémie avec intervention**
# Cette figure illustre l’évolution du nombre d’individus sains, infectieux et décédés au cours du 
# temps pour une simulation représentative avec intervention. Elle permet de suivre la progression 
# de l’épidémie génération par génération, d’identifier le pic d’infection et d’observer l’effet de 
# l’intervention sur la diminution du nombre d’individus infectieux et la limitation des décès.

# **Figure 3: Dynamique temporelle de l’épidémie sans intervention**
# Cette figure présente l’évolution du nombre d’individus sains, infectieux et décédés au cours du 
# temps pour une simulation sans intervention. Elle sert de référence pour comparer la dynamique 
# naturelle de l’épidémie en absence de contrôle. Elle permet notamment d’observer une propagation 
# plus rapide de la maladie, un pic d’infection plus élevé et une augmentation plus importante du 
# nombre de décès.

# Les Figures 2 et 3 permettent une comparaison directe de la dynamique temporelle de l’épidémie entre 
# les deux scénarios, tandis que la Figure 1 synthétise l’effet global de l’intervention sur plusieurs 
# simulations. Ensemble, ces figures permettent d’évaluer à la fois l’effet moyen de l’intervention et 
# la variabilité des trajectoires épidémiques.

# Les résultats du modèle sont obtenus à partir de 50 simulations indépendantes réalisées avec et sans 
# intervention. Cette réplication permet de tenir compte de la variabilité inhérente au caractère stochastique 
# du modèle.

# En moyenne, le nombre de décès observé en présence d’une intervention est de 2425.32 individus, alors qu’il 
# atteint 2739.9 individus en absence d’intervention. Cela correspond à une réduction moyenne d’environ 315 
# décès lorsque le dépistage et la vaccination sont appliqués. Ainsi, l’intervention permet de diminuer la 
# mortalité globale au sein de la population simulée.

# La variabilité des résultats diffère également entre les deux scénarios. La variance du nombre de décès 
# est de 389024.06 avec intervention, contre 669576.95 sans intervention. Les simulations avec intervention 
# présentent donc une dispersion plus faible, indiquant des résultats plus homogènes d’une répétition à l’autre.

# Le budget moyen restant à la fin des simulations avec intervention est de 1578.28, ce qui indique que la 
# majorité des ressources disponibles est utilisée au cours de l’épidémie. En revanche, dans le scénario sans 
# intervention, le budget demeure constant, puisqu’aucune dépense n’est effectuée.

# Les dynamiques temporelles de l’épidémie sont illustrées par les courbes représentant l’évolution du nombre 
# d’individus sains, infectieux et décédés. Sans intervention, le nombre d’individus infectieux augmente rapidement, 
# entraînant une diminution rapide du nombre d’individus sains et une accumulation importante de décès. 
# Avec intervention, l’augmentation du nombre d’infectieux est plus progressive, et la diminution des individus 
# sains est ralentie, ce qui se traduit par un nombre total de décès plus faible.

# L’analyse des événements d’infection révèle un total de 121 216 transmissions enregistrées au cours des simulations 
# avec intervention. De plus, 73 344 agents infectieux uniques sont impliqués dans ces événements. La distribution 
# du nombre d’infections causées par individu montre que la majorité des agents infectent un petit nombre d’individus, 
# tandis qu’une minorité contribue à un nombre plus élevé de transmissions.

# Dans l’ensemble, ces résultats mettent en évidence l’impact mesurable de l’intervention sur la réduction de la 
# mortalité et sur la dynamique de propagation de la maladie, tout en illustrant la variabilité importante des 
# trajectoires épidémiques entre simulations.


# ## Analyse events
all_events = vcat([r[3] for r in results_with]...)

infxn_by_uuid = countmap([e.from for e in all_events])
nb_inxfn = countmap(values(infxn_by_uuid))

println("Nombre d'agents infectieux uniques : ", length(infxn_by_uuid))

# ## Inclure du code

# Tous les fichiers dans le dossier `code` peuvent être ajoutés au travail
# final. C'est par exemple utile pour déclarer l'ensemble des fonctions du
# modèle hors du document principal.

# Le contenu des fichiers est inclus avec `include("code/nom_fichier.jl")`.

# Attention! Il faut que le code soit inclus au bon endroit (avant que les
# fonctions déclarées soient appellées).

include("code/01_test.jl")


# # Discussion

# On peut aussi citer des références dans le document `references.bib`, qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# `@` suivi de la clé de citation. Par exemple: @ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.