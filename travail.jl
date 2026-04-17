# ---
# title: Simulation individu-centrée de la propagation d’une maladie infectieuse et évaluation d’une stratégie de dépistage et vaccination
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
# ---

# # Introduction
# Lors de propagation rapide de maladies infectieuses, comme les épidémies, il est essentiel
# de comprendre comment différentes stratégies d’intervention permettent de limiter la transmission 
# et la mortalité. En effet, cela représente un enjeu majeur en santé publique, surtout lorsque 
# les individus infectés peuvent transmettre la maladie sans présenter de symptômes. 
# Par exemple, dans le cas de la covid, le contrôle de la propagation de la maladie a été 
# compliqué du aux individus porteurs asymptomatiques (Zhang et al, 2022).

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
# (Henao-Restrepo et al, 2017). Cela permet de créer une barrière immunitaire. Cette stratégie
# permet de limiter localement la propagation de la maladie. Cette stratégie à été utilisé lors
# des épidémies d’Ébola et a été efficace pour contrôler la transmission de la maladie
# (Henao-Restrepo et al, 2017)

# Pour le dépistage, il existe également différentes stratégies, comme le dépistage massif
# et le dépistage ciblé.

# Le dépistage massif consiste à tester une grande partie de la population, ce qui permet de
# détecter un grand nombre de cas (Shen et al, 2021). Cependant, cette stratégie nécessite
# beaucoup de ressource et peut donc être difficile à garder en place sur le long terme
# (Shen et al, 2021). 

# Le dépistage ciblé consiste à tester en priorité certains individus, par exemple les contacts
# de cas confirmés (Kretzschmar et al, 2020). Cette stratégie utilise plus efficacement les
# ressources sur le long terme.

# Dans le modèle, une stratégie combinant un dépistage ciblé progressif ainsi qu’une
# vaccination en anneau a été choisie. Une fois les personnes identifiées comme infectieux,
# les tests de dépistages sont réaliser. Le dépistage ciblé permet de tester les individus
# porches des cas confirmés, ce qui permet une meilleure gestion des ressources. Une fois les
# cas détectés, la vaccination des contacts permet de limiter la propagation locale, ce qui
# créer une barrière immunitaire. Cette stratégie a été choisie puisqu’elle permet de mieux
# utiliser les ressources. En effet, en testant uniquement les individus proches des cas
# confirmés, les gens les plus probable d’être infectés sont testé, ce qui permet une meilleure
# gestion des ressources. La combinaison du dépistage ciblé progressif et de la vaccination
# en anneau constitue donc une stratégie efficace et réaliste pour limiter la propagation
# d’une maladie infectieuse.

# # Présentation du modèle

# Nous utilisons un modèle de simulation individu-centré pour représenter la propagation d’une maladie infectieuse 
# dans une population. Les modèles individu-centrés permettent de représenter explicitement les interactions entre 
# individus et d’évaluer l’effet de différentes stratégies de contrôle, notamment le dépistage et la vaccination 
# (Adam & Arduin, 2023). Chaque individu est représenté par un agent évoluant dans un espace discret bidimensionnel 
# (lattice) défini par des coordonnées allant de -50 à 50 sur les axes x et y, avec des conditions de bord de type tore. 
# La population initiale est composée de 3750 agents distribués aléatoirement dans cet espace.

# Chaque agent possède plusieurs caractéristiques : une position spatiale, un état sanitaire (infectieux ou sain), un 
# compteur de durée d’infection (fixé à 21 générations), un statut de test, ainsi qu’un statut vaccinal avec un délai 
# avant efficacité. Lorsqu’un agent devient infectieux, il le reste jusqu’à la fin de sa période d’infection, après quoi 
# il est retiré de la population, ce qui correspond à un décès. La population est initialement entièrement naïve, avec 
# un seul cas index infectieux.

# À chaque pas de temps (tick), les agents se déplacent aléatoirement dans la grille. La transmission de la maladie se 
# produit lorsqu’un agent infectieux partage la même cellule qu’un agent sain. Dans ce cas, la probabilité de transmission 
# est fixée à 0,4, sauf si l’agent sain est vacciné et que le délai d’efficacité (2 générations) est atteint.

# L’intervention sanitaire combine dépistage ciblé et vaccination en anneau, et est déclenchée après le premier décès 
# observé dans la population. Elle est soumise à un budget total de 21 000 unités, avec un coût de 4 unités par test et de 
# 17 unités par vaccination.

# Le dépistage est réalisé de manière périodique, tous les 5 ticks, une fois l’intervention déclenchée. À chaque cycle, 
# jusqu’à 150 agents sont sélectionnés aléatoirement dans la population et testés à l’aide de tests antigéniques rapides. 
# Ces tests présentent une probabilité de détection de 95 % pour les individus infectieux, ainsi qu’un faible taux de faux 
# positifs (5 %) pour les individus sains (Lasser et al., 2021 ; Tran et al., 2023). Chaque agent ne peut être testé qu’une 
# seule fois.

# Lorsqu’un individu est testé positif, une intervention locale est mise en place : tous les agents présents dans la même 
# cellule spatiale sont considérés comme des contacts proches et peuvent être vaccinés, sous contrainte budgétaire. Cette 
# approche correspond à une stratégie de dépistage couplée à une vaccination ciblée de type « vaccination en anneau », 
# visant à limiter la propagation autour des cas détectés.

# Le vaccin est administré uniquement aux contacts des cas positifs, si le budget le permet. Il est supposé parfaitement 
# efficace après un délai de deux générations. Une fois ce délai écoulé, les individus vaccinés deviennent complètement 
# immunisés et ne peuvent plus être infectés ni transmettre la maladie.

# La simulation se poursuit jusqu’à l’extinction de la maladie (absence d’individus infectieux) ou jusqu’à un maximum de 
# 2000 générations. Les principales variables suivies sont le nombre d’individus sains, infectieux et décédés, ainsi que 
# l’évolution du budget, permettant d’évaluer l’efficacité de l’intervention.


# # Implémentation

# ## Packages nécessaires

# Initialisation
using Random
using CairoMakie
using UUIDs ## nous sommmes conscientes que cela doit se trouver dans "Project.toml", mais nous n'avons pas eu de succès avec le code en l'enlevant du "main"

Random.seed!(123456)
CairoMakie.activate!(px_per_unit=6.0) ## Permet de configurer l'affichage des figures

# Identifiants uniques pour chaque agent

UUIDs.uuid4()

# ## Création des types

# Structure d'un agent (individu)

Base.@kwdef mutable struct Agent ## Défini un individu modifiable
    x::Int64 = 0 ## position x
    y::Int64 = 0 ## position y
    clock::Int64 = 21 ## durée de la maladie est de 21 jours
    infectious::Bool = false ## état infectieux
    tested::Bool = false ## testé ou non
    ## Vaccination
    vaccinated::Bool = false        ## vacciné ou non
    days_after_vax::Int64 = 0        ## délai avant efficacité
    id::UUIDs.UUID = UUIDs.uuid4() ## identifiant unique
end

# Structure du paysage (espace de simulation)

Base.@kwdef mutable struct Landscape ## Définit les limites de l'espace
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

# Crée plusieurs agents (n)

Random.rand(::Type{Agent}, L::Landscape, n::Int64) = [rand(Agent, L) for _ in 1:n]

# Cette fonction nous permet donc de générer un nouvel agent dans un paysage:

rand(Agent, L)

# Mais aussi de générer plusieurs agents:

rand(Agent, L, 3)

# Déplacement des agents de -1, 0 ou +1 en x et y

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

budget = 21000      ## budget total
cost_vaccine = 17   ## coût vaccinated
cost_test = 4       ## coût test

intervention_started = false   ## commence après premier décès
first_death = false
deaths = 0                     ## compteur de morts

# Fonction de vaccination

function vaccinate!(agent::Agent) ## Rend un agent vacciné 
    if !agent.vaccinated
        agent.vaccinated = true
        agent.days_after_vax = 0   ## Juste après la vaccination, l'agent n'est pas encore protégé
    end
end

# ## Simulation principale
Base.@kwdef struct InfectionEvent
    tick::Int                ## Moment (temps) où l'infection a lieu
    from::UUIDs.UUID         ## identifiant de l'agent infectant
    to::UUIDs.UUID           ## idendifiant de l'agent infecté
    x::Int                   ## position x où l'infection se produit
    y::Int                   ## position y où l'infection se produit
end

# Fonction principale de la simulation
function run_simulation(L::Landscape, n::Int, budget_total; with_intervention=true)

    population = make_population(L, n) ## Création d'une population de n agents dans le paysage
    rand(population).infectious = true ## Choisit un agent au hasard et le rend infectieux (cas index)

    events = InfectionEvent[]  ## Liste pour stocker tous les événements d'infection

    budget = budget_total ## Initialisation du budget disponible pour les interventions
    intervention_started = false ## Booléen indiquant si l'intervention a commencé
    first_death = false ## Booléen pour détecter le premier décès
    deaths = 0 ## Compteur total de décès

    tick = 0 ## Temps initial
    maxlength = 2000 ## Durée maximale de la simulation (2000 générations)
    ## Liste pour stocker l'évolution des états
    S = Int[] ## nombre de sains
    I = Int[] ## nombre d'infectieux
    D = Int[] ## nombre de morts cumulés
    B = Int[] ## historique du budget

    while (length(infectious(population)) > 0) && (tick < maxlength) ## Boucle principale: continue tant qu'il y a des infectés et que la durée max n'est pas atteinte

        tick += 1 ## Incrément du temps

        ## Déplacement : chaque agent se déplace aléatoirement dans le paysage
        for agent in population
            move!(agent, L)
        end

        ## Mise à jour délai vaccin : si un agent est vacciné, on augmente le temps depuis la vaccination
        for agent in population
            if agent.vaccinated
                agent.days_after_vax += 1
            end
        end

        ## Infection 
        for agent in Random.shuffle(infectious(population)) ## On parcourt les agents infectieux dans un ordre aléatoire
            for neighbor in healthy(incell(agent, population)) ## On regarde les voisins sains dans la même cellule 

                immune = neighbor.vaccinated && neighbor.days_after_vax >= 2 ## Vérifie si le voisin est protégé par le vaccin (après 2 jours)

                if rand() <= 0.4 && !immune ## Probabilité de transmission de 40% si non immunisé
                    neighbor.infectious = true ## Le voisin devient infectieux
                    push!(events, InfectionEvent(tick, agent.id, neighbor.id, neighbor.x, neighbor.y)) ## On enregistre l'événement d'infection
                end
            end
        end

        ## Progression maladie : Chaque agent infecté voit son compteur de maladie diminuer
        for agent in infectious(population)
            agent.clock -= 1
        end

        ## Décès
        before = length(population) ## Nombre d'agents avant suppression des morts
        population = filter(a -> a.clock > 0, population) ##  # On retire les agents dont le temps de vie (clock) est écoulé
        after = length(population) ## Nombre d'agents après suppression

        ## Déclenche intervention : Si c'est le premier décès, on active l'intervention
        if !first_death && after < before
            intervention_started = true
            first_death = true
        end

        deaths += (before - after) ## Mise à jours du nombre total de décès

        ## Intervention (activable)
        if with_intervention && intervention_started && budget > 0 && tick % 5 == 0 ## Vérifie que l'intervention est activée, qu'elle commence après un décès, qu'il reste du budget et que l'on est à un tuck multiple de 5 (tests prériodiques)

            ## Faire une liste vide pour stocker les agents à risque afin de tester ceux qui sont autour des infectés
            a_risque = Agent[]

            for inf in infectious(population) ## Parcourt tous les agents infectés
                for neighbor in incell(inf, population) ## Parcourt tous les agents dans la même cellule que l'infecté
                    push!(a_risque, neighbor) ## Ajoute chaque voisin dans la liste des agetns à risque
                end
            end

            a_risque = unique(a_risque) # Vient supprimer les doublons (un agent peut apparaître plusieurs fois)
            candidates = isempty(a_risque) ? population : a_risque # Si aucun agent à risque, on prend toute la population, sinon on teste que les agents à risque
            n_tests_per_tick = 150   ## 150 agents testés à ce tick
            candidates = Random.shuffle(population)[1:min(n_tests_per_tick, length(population))] ## Sélection random d'agents et en garde maximum 50

            for agent in candidates ## Parcours les 150 agents sélectionnés aléatoirement

                if !agent.tested && budget >= 4 ## Si l'agent n'a pas encore été testé et que le budget le permet
                    budget -= 4 ## Coût du test
                    agent.tested = true ## Marque l'agent comme testé

                    positive = agent.infectious ? rand() <= 0.95 : rand() <= 0.05 ## Résultat du test (avec faux positifs et faux négatifs)

                    if positive ## Si test positif
                        agent.tested = true # Isolement des cas positifs pour empêcher la propagation de la maladie
                        for neighbor in incell(agent, population) ## On vaccine tous les agents dans la même cellule
                            if budget >= 17 && !neighbor.vaccinated ## Si budget suffisant et agent non vacciné
                                vaccinate!(neighbor) ## Vacciation du voisin
                                budget -= 17 ## Coût du vaccin
                            end
                        end
                    end
                end

                if budget <= 0 ## Arrêt si budget épuisé
                    break
                end
            end
        end
        ## Suivi des états
        push!(S, length(healthy(population))) ## Enregistre le nombre d'agents sains
        push!(I, length(infectious(population))) ## Enregistre le nombre d'agents infectieux
        push!(D, deaths) ## Enregistre le nombre total de morts
        push!(B, budget) ## Enregistre la progression du budget
    end

    ## Retourne: nombre total de morts, budget restant, liste des infections, historique S/I/D/B
    return deaths, budget, events, S, I, D, B
end

# ## Analyse des résultats
using Statistics ## nous sommmes conscientes que cela doit se trouver dans "Project.toml", mais nous n'avons pas eu de succès avec le code en l'enlevant du "main"

n_runs = 50;

results_with = [run_simulation(L, 3750, 21000; with_intervention=true) for _ in 1:n_runs];
results_without = [run_simulation(L, 3750, 21000; with_intervention=false) for _ in 1:n_runs];

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

# ## Statistiques

mean_with = mean(deaths_with) ## moyenne nombre décès AVEC intervention
mean_without = mean(deaths_without) ## moyenne nombre décès SANS intervention

var_with = var(deaths_with) ## Calcule la variance des décès AVEC intervention
var_without = var(deaths_without) ## Calcule la variance des décès SANS intervention

std_with = std(deaths_with)       # écart-type AVEC intervention
std_without = std(deaths_without) # écart-type SANS intervention

println("\n RÉSULTATS SUR ", n_runs, " SIMULATIONS")

println("\nAVEC intervention :")
println("Moyenne des morts = ", mean_with)
println("Variance = ", var_with)

println("\nSANS intervention :")
println("Moyenne des morts = ", mean_without)
println("Variance = ", var_without)

println("\nGAIN (morts évités) = ", mean_without - mean_with)

## Graphique comparatif
f1 = Figure();

ax = Axis(f1[1, 1];
    xlabel="Simulation",
    ylabel="Nombre de morts",
    title="Comparaison des décès avec et sans intervention"
)

scatter!(ax, 1:length(deaths_with), deaths_with, label="Avec intervention", color=:blue)
scatter!(ax, 1:length(deaths_without), deaths_without, label="Sans intervention", color=:red)

axislegend(ax)
f1

# **Figure 1: Comparaison des décès**
# Cette figure présente, pour chacune des 50 simulations, le nombre total de décès observés 
# avec et sans intervention. Chaque point correspond à une simulation indépendante. Cette 
# figure permet de visualiser la variabilité des résultats entre les répétitions et de comparer 
# directement l’effet global de l’intervention sur la mortalité. Elle met en évidence la dispersion 
# des résultats ainsi que la tendance générale à une réduction du nombre de décès en présence 
# d’une intervention.

println("Moyenne morts avec intervention = ", mean(deaths_with))
println("Moyenne morts sans intervention = ", mean(deaths_without))

f2 = Figure();

## Graphiqe avec intervention
ax1 = Axis(f2[1, 1],
    title = "Évolution AVEC intervention",
    xlabel = "Temps",
    ylabel = "Population"
)


lines!(ax1, 1:length(S_with), S_with, label="Sains")
lines!(ax1, 1:length(I_with), I_with, label="Infectieux")
lines!(ax1, 1:length(D_with), D_with, label="Décédés")

axislegend(ax1)

## Graphique sans intervention
ax2 = Axis(f2[2, 1],
    title = "Évolution SANS intervention",
    xlabel = "Temps",
    ylabel = "Population"
)

lines!(ax2, 1:length(S_without), S_without, label="Sains")
lines!(ax2, 1:length(I_without), I_without, label="Infectieux")
lines!(ax2, 1:length(D_without), D_without, label="Décédés")

axislegend(ax2)

f2

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

## Graphique pour le budget (simulation sans intervention)
sim_with = run_simulation(L, 3750, 21000; with_intervention=true)

budget_time = sim_with[7]  # 7e sortie = budget

f_budget = Figure() ## Figure

ax = Axis(
    f_budget[1, 1],
    xlabel = "Temps (ticks)",
    ylabel = "Budget restant",
    title = "Évolution du budget pendant la simulation"
)

lines!(ax, 1:length(budget_time), budget_time)

f_budget

# **Figure 4: Épuisement du budget lors de l'épidémie avec intervention**
# Cette figure montre l’évolution du budget restant au cours du temps. Le budget diminue rapidement au début, 
# indiquant une forte utilisation des ressources lors des premières phases de l’épidémie. Ensuite, la décroissance 
# ralentit, puis atteint un plateau autour de 3 000, ce qui suggère qu’une partie du budget n’est pas utilisée
# car l’épidémie est contrôlée avant la fin de la simulation.

# Les résultats reposent sur 50 simulations indépendantes avec et sans intervention, afin de tenir compte du 
# caractère stochastique du modèle. En moyenne, le nombre de décès est de 2309,34 avec intervention contre 2910,88 
# sans intervention, soit une réduction d’environ 601,54 décès. L’intervention permet donc de diminuer significativement 
# la mortalité. La variabilité diffère toutefois entre les scénarios : la variance est plus élevée avec intervention 
# (≈ 753 051) qu’en absence d’intervention (≈ 181 697), indiquant des résultats plus hétérogènes selon les simulations.
# Le budget moyen restant avec intervention est d’environ 4377, ce qui suggère que toutes les ressources ne sont pas 
# toujours utilisées avant la fin de l’épidémie. Globalement, ces résultats montrent que l’intervention ralentit la 
# propagation de la maladie et réduit la mortalité, malgré une variabilité importante entre les simulations.

# ## Discussion
# Les résultats obtenus montrent que l’intervention basée sur le dépistage et la vaccination ciblée permet de 
# réduire de manière notable le nombre total de décès par rapport au scénario sans intervention. En moyenne, 
# environ 601,54 décès sont évités, ce qui indique que la stratégie mise en place est efficace pour limiter 
# la propagation de la maladie. Cette réduction s’explique par le fait que le dépistage permet d’identifier 
# rapidement une partie des individus infectieux et d’interrompre les chaînes de transmission, tandis que la 
# vaccination des contacts proches diminue la probabilité de transmission locale (Adam & Arduin, 2023).

# La stratégie de vaccination implémentée correspond à une approche de vaccination en anneau, dans laquelle 
# les individus en contact avec un cas détecté sont vaccinés de manière ciblée (Henao-Restrepo et al., 2017). 
# Cette stratégie vise à créer une barrière immunitaire autour des foyers d’infection afin de limiter leur 
# expansion (Krauland et al., 2026). Dans le modèle, cette approche permet de contenir localement la propagation, 
# mais son efficacité dépend fortement de la rapidité du dépistage et de la capacité à intervenir avant que 
# l’infection ne se diffuse à d’autres cellules.

# Cependant, l’intervention ne permet pas d’éliminer complètement l’épidémie. Plusieurs facteurs expliquent cette 
# efficacité partielle. D’une part, les tests antigéniques ne sont pas parfaitement sensibles, ce qui signifie 
# qu’une fraction des individus infectieux n’est pas détectée et continue de transmettre la maladie (Adam & Arduin, 2023 
# ; Krauland et al., 2026). D’autre part, le vaccin n’est pas immédiatement efficace; le délai de deux générations avant 
# l’acquisition de l’immunité laisse une période durant laquelle les individus vaccinés restent susceptibles. 
# Enfin, la limitation du budget restreint le nombre total de tests et de vaccinations pouvant être réalisés, 
# ce qui empêche une couverture complète de la population.

# L’analyse du budget apporte un éclairage important sur l’efficacité de l’intervention. En moyenne, un budget restant
# d’environ 4377,18 est observé à la fin des simulations, ce qui suggère que l’épidémie peut parfois s’éteindre avant 
# que toutes les ressources ne soient utilisées. Cela indique que le facteur limitant n’est pas uniquement la quantité 
# de ressources disponibles, mais également le moment du déclenchement de l’intervention et la dynamique de propagation 
# de la maladie. Toutefois, dans certaines simulations, une consommation plus rapide du budget limite la capacité à 
# tester et vacciner de nouveaux individus au cours de l’épidémie. Cela met en évidence l’importance d’une allocation 
# efficace et précoce des ressources pour maximiser l’impact des interventions.

# Les résultats montrent également une forte variabilité entre les simulations, tant avec que sans intervention. Cette 
# variabilité est particulièrement marquée dans le scénario avec intervention, où la variance des décès (≈ 753 051) est 
# plus élevée que sans intervention (≈ 181 697). Cela indique que, bien que l’intervention réduise le nombre moyen de
# décès, son effet est plus hétérogène selon les simulations. Dans certains cas, l’épidémie est rapidement contrôlée, 
# tandis que dans d’autres, elle se propage largement malgré les mesures mises en place. Cette variabilité est liée au 
# caractère stochastique du modèle, notamment dans les déplacements des agents et les événements de transmission 
# (Adam & Arduin, 2023 ; Shakiba et al., 2021).

# D’un point de vue épidémiologique, ces résultats sont cohérents avec les observations réelles. Ils soulignent 
# l’importance des stratégies combinant dépistage et vaccination ciblée pour contrôler une épidémie (Adam & Arduin, 2023). 
# Le modèle met également en évidence le rôle crucial du délai d’action des interventions; une réponse tardive ou 
# une immunité retardée peut limiter l’efficacité globale des mesures (Krauland et al., 2026). Par ailleurs, 
# la distribution du nombre d’infections secondaires suggère une propagation hétérogène, où certains individus 
# infectent plusieurs autres, tandis que la majorité en infecte peu. Ce phénomène est caractéristique des dynamiques 
# de type « super-propagation », soit des situations où un petit nombre d’individus infectés est responsable d’un 
# nombre disproportionné de transmissions, contribuant fortement à la diffusion globale de la maladie 
# (Karen et al., 2022 ; Shakiba et al., 2021).

# Malgré les résultats obtenus, ce modèle présente plusieurs limites qui devraient être prises en compte lors de 
# l’interprétation des résultats. Tout d’abord, la dynamique de transmission de la maladie utilisée dans cette 
# simulation est très simplifiée. Par exemple, la probabilité d’infection est fixée à 0,4 pour tous les individus 
# et toutes les interactions. Cependant, dans les situations réelles, la transmission dépend de plusieurs facteurs, 
# comme la charge virale, la durée de contact et les caractéristiques individuelles (âge, état de santé, etc.) 
# (Blaser et al., 2014). De plus, ces facteurs varient selon la maladie considérée : par exemple, la transmission 
# du VIH dépend fortement de la charge virale (Blaser et al., 2014), tandis que celle de la COVID-19 dépend davantage 
# de la proximité et de la durée des contacts (Karia et al., 2020). Cette simplification limite donc la capacité du 
# modèle à refléter fidèlement la complexité des dynamiques réelles.

# Une autre limite du modèle est que les agents se déplacent de manière aléatoire dans un environnement homogène. 
# Cette modélisation ne prend pas en compte les structures sociales réelles. En pratique, les individus interagissent 
# davantage au sein de groupes spécifiques (famille, école, travail), ce qui influence fortement la transmission des 
# maladies (Mossong et al., 2008). L’absence de ces structures peut conduire à une mauvaise estimation de la vitesse 
# et des schémas de propagation.

# De plus, dans ce modèle, la maladie est toujours fatale, ce qui constitue une simplification importante. La plupart 
# des maladies infectieuses présentent des taux de mortalité variables, avec une proportion d’individus qui guérissent
# et développent une immunité naturelle (Szomolanyi et al., 2010). L’absence d’un état « rétabli » empêche de représenter 
# ces dynamiques.

# En ce qui concerne les interventions, certaines simplifications ont également été faites. Les tests sont uniquement
# administrés aux individus ayant été en contact avec un cas confirmé, alors que, dans la réalité, certaines populations 
# sont plus à risque (Liu et al., 2021 ; Chapman et al., 2025). De plus, le vaccin est supposé parfaitement efficace 
# après deux générations, ce qui ne reflète pas les variations réelles d’efficacité vaccinale (Zachreson et al., 2023). 
# Le modèle ne prend pas non plus en compte les mutations virales ni la nécessité de doses multiples (Behl et al., 2022).

# Par ailleurs, l’intervention est déclenchée seulement après le premier décès, ce qui correspond à une réponse 
# relativement tardive. Une intervention plus précoce, dès la détection des premiers cas, pourrait améliorer 
# significativement les résultats (Kucharski et al., 2020). Il serait donc pertinent d’explorer différents scénarios 
# de déclenchement.

# Finalement, le caractère stochastique du modèle permet de représenter l’incertitude inhérente aux dynamiques 
# épidémiques, mais il introduit également une forte variabilité des résultats. Bien que cette variabilité soit 
# partiellement prise en compte par la répétition des simulations, une analyse plus approfondie permettrait 
# d’identifier les paramètres les plus influents.

# On peut aussi citer des références dans le document `references.bib`, qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# `@` suivi de la clé de citation. Par exemple: @ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.
