# CoRL2018Constraint
Code for CoRL2018: Inferring geometric constraints in human demonstrations

## Testing
* Use the Dockerfile to build a working environment. 
* The constraint equations are generated using sympy code generation. Each constraint model autogenerates into it's own python module. 
* The Dockerfile will generate all these files for you. 
```
python constraint_equations/geometric_constraints.py
```
* ```Run test.py``` to see an example. 