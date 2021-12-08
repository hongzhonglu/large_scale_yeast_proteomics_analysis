def exchange_ecYeast(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem


# function copy from strain_design repo
def ecYeastMinimalMedia(model):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """
    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if '_REV' in x and 'exchange' in y]
    # first block any uptake
    for i, x in enumerate(exchange_rxn):
        rxn0 = exchange_rxn[i]
        #print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    #Allow uptake of essential components
    model.reactions.get_by_id("r_1654_REV").upper_bound = 10000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861_REV").upper_bound = 10000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100_REV").upper_bound = 10000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992_REV").upper_bound = 10000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005_REV").upper_bound = 10000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060_REV").upper_bound = 10000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832_REV").upper_bound = 10000 #H+ exchange (reversible)
    return model


def chemostatSimulation(model0, D0):
    """
    This funcion is used to simulate the chemostat growth of yeast
    :param model0:
    :param D0:
    :return:
    """
    growth = D0
    with model0:
        model0 = ecYeastMinimalMedia(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    return solution3


def simulationCompare(model_in, growth_in, objective):

    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution00 = model0.optimize()
        GR = solution00.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution1 = model0.optimize()

    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()

    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()

    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = objective
        solution4 = model0.optimize()
    return solution1, solution2, solution3, solution4


