__version__ = "0.5.1"

algo_name_aliases = {
        # savings heuristics
        "ps":"ps",   "cw64-ps":"ps",
                     "cw":"ps",
                     "parallelsavings":"ps",
                     "clarkewright":"ps",
        "gpl":"gpl", "ga67-ps":"gpl","ga67-ps|pi":"gpl","ga67-ps|lambda":"gpl",
                     "gaskellpi":"gpl", "gaskelllambda":"gpl",
        "ss":"ss",   "we64-ss":"ss", "sequentialsavings":"ss",
        
        "gps":"gps", "pa88-ps":"gps","generalizedsavings":"gps",
                     "generalizedparallelsavings":"gps",
                     "paessens":"gps",
        "ims":"ims", "hp76-ps":"ims","hp76-ps|ims":"ims",
                     "suppression":"ims",
                     "iterativemergesuppression":"ims",
                     "mergesuppressionsavings":"ims",
        "ps2o":"ps2o", "rt79-ps":"ps2o", "cawlip":"ps2o",
                     "savingswith2opt":"ps2o", 
        
        # insertion heuristics
        "si":"si",   "ci":"si",
                     "cheapestinsertion":"si",
                     "sequentialinsertion":"si",
        "mj":"mj",   "mj76-si":"si",
                     "molejameson":"mj","molejamesoninsertion":"mj",
        "pi":"pi",   "parallelinsertion":"pi", "parallelcheapestinsertion":"pi",
        
        # maximum mathcing heuristics
        "mbsa":"mbsa",   "dv89-mbsa":"mbsa", "matching":"mbsa",
                     "maximummatching":"mbsa", "mm":"mbsa",
                     "desrochersverhoog":"mbsa",
        
        # the 2-phase heuristic
        "cmt":"cmt", "cmt79-2p":"cmt","tp":"cmt","twophase":"cmt",
                     "cmt2p":"cmt","cmt2phase":"cmt","cmttwophase":"cmt",
        
        "sn":"sn", "nn":"sn","snn":"sn", # parallel nearest neighbour
        "pn":"pn", "pnn":"pn", # parallel nearest neighbour
        "ty":"ty",   "ty68-snn":"ty","tyagi":"ty", # tyagi nearest neighbour
        
        "swp":"swp", "sweep":"swp", # plain sweep
        "wh":"wh",   "wh72-swls":"wh", "whs":"wh", "whswp":"wh",
                     "wrenholliday":"wh", "wrenhollidaysweep":"wh", 
        "gm":"gm",   "gm74-swri":"gm", "gms":"gm", "gmswp":"gm",
                     "gilletmiller":"gm", "gilletmillersweep":"gm",
        
        "rfcs":"rfcs", "be83-rfcs":"rfcs","routefirstclustersecond":"rfcs",
                       "beasley":"rfcs",
                       "newtonthomas":"rfcs",
        
        # generalized assignment problem heuristic
        "gap":"gap", "fj81-gap":"gap", "fisherjaikumar":"gap",
        
        # set covering "petal" heuristic
        "ptl":"ptl", "fr76-ptl":"ptl", "petal":"ptl", "fosterryan":"ptl",
        
        # lagrangian relaxation 3-opt* heuristic
        "lr3o":"lr3o", "sg82-lr3opt":"lr3o", "lr3opt":"lr3o",
        
        # TO ENABLE THEM ALL!
        "all":"all",
        "classical":"classical", "classic":"classical",
}

def get_algorithms(names):
    has_gurobi = True
    try:
        import gurobipy
    except:
        has_gurobi = False
    
    # get algorithms
    algos = []
    for algo_name in names:
        algo_name = algo_name.lower() # ignore case
        if algo_name in algo_name_aliases:
            # translate name to standard abbreviation
            algo_name = algo_name_aliases[algo_name]

            # import only those that are needed, when needed
            if algo_name in ["ps", "all", "classical"]:
                #"ps":"Clarke & Wright (1964) parallel savings heuristic",
                from verypy.classic_heuristics.parallel_savings import get_ps_algorithm
                algos.append( ("ps",)+get_ps_algorithm() )
            if algo_name in ["ps2o", "all"]:
                #"ps2o":"Robbins and Turner (1979) CAWLIP savings algorithm",
                from verypy.classic_heuristics.cawlip_savings import get_ps2o_algorithm
                algos.append( ("ps2o",)+get_ps2o_algorithm() )
            if algo_name in ["gpl", "all", "classical"]:
                #"gpl":"Gaskell (1967) best of pi/lambda savings heuristic",
                from verypy.classic_heuristics.gaskell_savings import get_gs_algorithm
                algos.append( ("gpl",)+get_gs_algorithm() )
            if algo_name in[ "ss", "all", "classical"]:
                #"ss":"Webb (1964) sequential savings algorithm"
                from verypy.classic_heuristics.sequential_savings import get_ss_algorithm
                algos.append( ("ss",)+get_ss_algorithm(lambda_multiplier=1.0) )
            if algo_name in ["gps", "all", "classical"]:
                #"gps":"Paessens (1988) parallel savings heuristic with generalized savings function and M4 strategy",
                from verypy.classic_heuristics.paessens_savings import get_gps_algorithm
                algos.append( ("gps",)+get_gps_algorithm() )
                pass
            if algo_name in ["ims", "all", "classical"]:
                #"ims":"Holmes & Parker (1976) iterative parallel savings merge suppression heuristic",
                from verypy.classic_heuristics.suppression_savings import get_ims_algorithm 
                algos.append( ("ims",)+get_ims_algorithm() )
            if algo_name in ["si", "all"]:
                #"si":"Mole & Jameson (1976) sequential insertion heuristic without local search"),
                from verypy.classic_heuristics.cheapest_insertion import get_si_algorithm 
                algos.append( ("si",)+get_si_algorithm() )        
            if algo_name in ["mj", "all", "classical"]:
                #"mj":"Mole & Jameson (1976) sequential insertion heuristic with local search"),                
                from verypy.classic_heuristics.mole_jameson_insertion import get_mj_algorithm 
                algos.append( ("mj",)+get_mj_algorithm() )
            if algo_name in ["pi", "all"]:
                #"pi":"van Breedam (1995) parallel insertion heuristic"),
                from verypy.classic_heuristics.cheapest_insertion import get_pi_algorithm
                algos.append( ("pi",)+get_pi_algorithm() )
            if algo_name in ["mbsa", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [mbsa/DV89-MM] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:
                    #"mbsa":"Desrochers and Verhoog (1989) matching based savings algorithm"
                    from verypy.classic_heuristics.matchingvrp import get_mm_algorithm
                    algos.append( ("mbsa",)+get_mm_algorithm())
            if algo_name in ["cmt", "all", "classical"]:
                #"cmt":"Christofides, Mingozzi & Toth (1979) two phase heuristic"
                from verypy.classic_heuristics.cmt_2phase import get_cmt2p_algorithm
                algos.append( ("cmt",)+get_cmt2p_algorithm() )
            if algo_name in ["sn", "all"]:
                #"sn":"Sequential Nearest Neighbor construction heuristic"
                from verypy.classic_heuristics.nearest_neighbor import get_snn_algorithm
                algos.append( ("sn",)+get_snn_algorithm() )
            if algo_name in ["pn", "all"]:
                #"pn":"Parallel Nearest Neighbor construction heuristic"
                from verypy.classic_heuristics.nearest_neighbor import get_pnn_algorithm
                algos.append( ("pn",)+get_pnn_algorithm() )
            if algo_name in ["ty", "all", "classical"]:
                #"ty":"Tyagi (1968) Nearest Neighbor construction heuristic"
                from verypy.classic_heuristics.tyagi_nearest_neighbor import get_ty_algorithm
                algos.append( ("ty",)+get_ty_algorithm() ) 
            if algo_name in ["swp", "all"]:
                #"swp":"Plain Sweep algorithm without route improvement heuristics"
               from verypy.classic_heuristics.sweep import get_swp_algorithm
               algos.append( ("swp",)+get_swp_algorithm() )
            if algo_name in ["wh", "all", "classical"]:
                #"wh":"Wren and Holliday (1972) Sweep + Local Search algorithm"
               from verypy.classic_heuristics.wren_holliday_sweep import get_wh_algorithm
               algos.append( ("wh",)+get_wh_algorithm() ) 
            if algo_name in ["gm", "all", "classical"]:
                #"gm":"Gillett and Miller (1974) Sweep algorithm"
                from verypy.classic_heuristics.gillet_miller_sweep import get_gm_algorithm
                algos.append( ("gm",)+get_gm_algorithm() )
            if algo_name in ["rfcs", "all", "classical"]:
                #"rfcs":"Route-first-cluster-second heuristic of Beasley (1983)"
                from verypy.classic_heuristics.rfcs import get_rfcs_algorithm
                algos.append( ("rfcs",)+get_rfcs_algorithm()) 
            if algo_name in ["gap", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [gap/FJ81-GAP] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:
                    #"gap":"Fisher & Jaikumar (1981) generalized assignment problem heuristic"
                    from verypy.classic_heuristics.gapvrp import get_gap_algorithm
                    algos.append( ("gap",)+get_gap_algorithm() )
            if algo_name in ["ptl", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [ptl/FR76-1PTL] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:                    
                    #"ptl":"Foster and Ryan (1976) Petal algorithm"
                    from verypy.classic_heuristics.petalvrp import get_ptl_algorithm
                    algos.append( ("ptl",)+get_ptl_algorithm() )  
            if algo_name in ["lr3o", "all", "classical"]:
                #"lr3o":"Stewart and Golden (1982) LR3OPT heuristic"
                from verypy.classic_heuristics.lr3opt import get_lr3opt_algorithm
                algos.append( ("lr3o",)+get_lr3opt_algorithm() )  
        else:
            print(algo_name, "is not a valid algorithm name", file=sys.stderr)
    return algos