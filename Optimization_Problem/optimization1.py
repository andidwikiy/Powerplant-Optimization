#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 21:30:59 2024

@author: imran
"""

import pulp
import pandas as pd

# Constants and Parameters
disc = 0.1
years = 5
alpha = 0.25
PG = {'A': 150, 'B': 250, 'C': 100}
cost_inv = {'A': 300000, 'B': 350000, 'C': 250000}
cost_fuel = {'A': 20.409, 'B': 14, 'C': 25.953}
cost_fixed = {'A': 12000, 'B': 36000, 'C': 30000}
cost_var = {'A': 1, 'B': 3, 'C': 2.5}
ReqCap = {0: 600, 1: 720, 2: 840, 3: 960, 4: 1080}

# Read load data
df_demand = pd.read_csv('load data(1).csv')

# Initialize the model
model = pulp.LpProblem("Generator_Installation_Plan", pulp.LpMinimize)

# Decision Variables
P = pulp.LpVariable.dicts("Power", ((unit, hour, year) for unit in PG for hour in range(8760) for year in range(years)), lowBound=0)
X = pulp.LpVariable.dicts("total_Units", ((unit, year) for unit in PG for year in range(years)), lowBound=0, cat='Integer')
Z = pulp.LpVariable.dicts("New_Units", ((unit, year) for unit in PG for year in range(years)), lowBound=0, cat='Integer')

# Objective Function
Cinv = pulp.lpSum([(1 + disc) ** -year * cost_inv[unit] * PG[unit] * Z[unit, year] for unit in PG for year in range(years)])
Csalv = pulp.lpSum([(1 + disc) ** -years * cost_inv[unit] * PG[unit] * Z[unit, year] * (1 - alpha) ** (years - year) for unit in PG for year in range(years)])
Cfuel = pulp.lpSum([(1 + disc) ** -year * P[unit, hour, year] * cost_fuel[unit] for unit in PG for hour in range(8760) for year in range(years)])
Com = pulp.lpSum([(1 + disc) ** -year * cost_fixed[unit] * PG[unit] * X[unit, year] for unit in PG for year in range(years)]) + pulp.lpSum([(1 + disc) ** -year * P[unit, hour, year] * cost_var[unit] for unit in PG for hour in range(8760) for year in range(years)])

model += Cinv + Cfuel + Com - Csalv

# Constraints
for year in range(years):
    model += pulp.lpSum([PG[unit] * X[unit, year] for unit in PG]) >= ReqCap[year]
    for hour in range(8760):
        model += pulp.lpSum([P[unit, hour, year] for unit in PG]) == df_demand.loc[hour, str(year)]
        for unit in PG:
            model += 0 <= P[unit, hour, year] <= X[unit, year] * PG[unit]

for unit in PG:
    for year in range(1, years):
        model += Z[unit, year] == X[unit, year] - X[unit, year - 1]
    model += Z[unit, 0] == X[unit, 0]  # For the first year, new units are the total units

# Solve the model
model.solve()

# Calculate and print the number of new units installed for each generator type every year
for year in range(years):
    print(f"Year {year}:")
    for unit in PG:
        if year == 0:
            new_units = X[unit, year].varValue
        else:
            new_units = X[unit, year].varValue - X[unit, year - 1].varValue
        print(f"Unit {unit}: {max(0, new_units)} new units installed")

# Display total cost
print("Total Plan Cost =", pulp.value(model.objective))




