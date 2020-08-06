# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

# MODIFIED FROM http://examples.gurobi.com/traveling-salesman-problem/

from math import sqrt
from signal import signal, SIGINT, default_int_handler

from gurobipy import Model, GRB, quicksum, GurobiError
from config import MAX_MIP_SOLVER_RUNTIME, MIP_SOLVER_THREADS

# closures that include n
def _subtourelim(model, where):
  """  Callback - use lazy constraints to eliminate sub-tours """
  if where == GRB.callback.MIPSOL:
    # make a list of edges selected in the solution
    X = model.cbGetSolution(model._vars)
    n = int(sqrt(len(X)))
    selected = [(i,j) for i in range(n) for j in range(n) if X[(i,j)]>0.5]

    # find the shortest cycle in the selected edge list
    tour = _subtour(selected,n)
    if len(tour) < n:
      # add a subtour elimination constraint
      expr = quicksum(model._vars[tour[i], tour[j]]
          for i in range(len(tour))
              for j in range(i+1, len(tour)))
      model.cbLazy(expr <= len(tour)-1)

# closures that include n
def _subtour(edges,n):
  """ Given a list of edges, finds the shortest subtour """
  visited = [False]*n
  cycles = []
  costs = []
  selected = [[] for i in range(n)]
  for x,y in edges:
    selected[x].append(y)
  while True:
    current = visited.index(False)
    thiscycle = [current]
    while True:
      visited[current] = True
      neighbors = [x for x in selected[current] if not visited[x]]
      if len(neighbors) == 0:
        break
      current = neighbors[0]
      thiscycle.append(current)
    cycles.append(thiscycle)
    costs.append(len(thiscycle))
    if sum(costs) == n:
      break
  return cycles[costs.index(min(costs))]
      
def solve_tsp_gurobi(D, selected_idxs):
    #print     
    #print "============== GUROBI ============="
    
    n = len(selected_idxs)
    if selected_idxs[0]==selected_idxs[-1]:
        n = n-1
    
    # no need to invoke Gurobi for tiny TSP cases
    if n<=3:
        sol = None
        obj_f = 0.0
        if n>1:
            sol = list(selected_idxs)
            if sol[0]!=sol[-1]:
                sol.append(sol[0])            
            obj_f = sum( D[sol[i-1],sol[i]] for i in range(1,len(sol)) ) 
        return sol, obj_f 
        
    m = Model("TSP")
    m.params.OutputFlag = 0

    # Create variables    
    edgevars = {}
    
    for i in range(n):
       for j in range(i+1):
         from_node = selected_idxs[i]
         to_node = selected_idxs[j]
         edgevars[i,j] = m.addVar(obj=D[from_node, to_node], vtype=GRB.BINARY,
                              name='e'+str(i)+'_'+str(j))
         edgevars[j,i] = edgevars[i,j]
       m.update()
    
    # Add degree-2 constraint, and forbid loops
    for i in range(n):
      m.addConstr(quicksum(edgevars[i,j] for j in range(n)) == 2)
      edgevars [i,i].ub = 0
    m.update()
    
    # Optimize model
    m._vars = edgevars 
    m.params.LazyConstraints = 1
    m.setParam('TimeLimit', MAX_MIP_SOLVER_RUNTIME)
    m.setParam('Threads', MIP_SOLVER_THREADS)
    
    m.optimize(_subtourelim)

    # restore SIGINT callback handler which is changed by gurobipy
    signal(SIGINT, default_int_handler)

    status = m.Status
    if status == GRB.TIME_LIMIT:
        raise GurobiError(10023, "Gurobi timeout reached when attempting to solve TSP")
    elif m.Status == GRB.INTERRUPTED:
        raise KeyboardInterrupt()
    
    solution = m.getAttr('x', edgevars)
    selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5]
    cycles = _subtour(selected, n)
    assert len(cycles) == n
    
    # make the route always start from the 1st index
    sol = [selected_idxs[i] for i in cycles]+[selected_idxs[0]]
    obj_f = m.objVal
    return sol, obj_f

if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("gurobi", solve_tsp_gurobi)
