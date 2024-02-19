# Script that plots the data from the csv files and saves them as png files
# Takes the first column as the x-axis and the rest as the y-axis
# The labels are inferred from the column names on the first row after comments
from re import L
import sys
import os
from os import path
from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# If no filename passed, exit
if len(sys.argv) < 2:
	print("Please pass the filename as an argument.")
	print("Only files in the output folder can be plotted.")
	exit()

# Find the path of CSV file
basepath = path.dirname(__file__)
filepath = path.abspath(path.join(basepath, "..", "output", sys.argv[1]))

# Read the CSV file using pandas
data = pd.read_csv(filepath, comment='#', index_col=0)

# Compose the title from the comments
title = ""
description = "("
with open(filepath, 'r') as f:
	title = f.readline()[1:].strip().replace('_', ' ').capitalize()
	for line in f:
		if line[0] == '#':
			description += line[1:].strip() + ", "
	description = description[:-2] + ")"


#plt.figure(figsize=(12, 10))
sns.set_style("darkgrid")
ax = None


# ====== Plot for the error_iteration.csv ======
if sys.argv[1] == "saspso_error_iterations.csv":
	fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

	# error plot
	ax = sns.lineplot(ax=axes[0], data=data[['static_err','static_p_err','dynamic_err','dynamic_p_err']], linewidth=3, legend=True)
	ax.set_yscale('symlog', linthresh=1e-16)	# Needed for having a log scale showing 0
	ax.set_ylim(bottom=-1e-16)
	ax.set_ylabel("Error", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	#ax.legend(['Townsend', 'Townsend Parallel', 'Gomez-Levy', 'Gomez-Levy Parallel'])

	# constraint violation plot
	ax = sns.lineplot(ax=axes[1], data=data[['static_viol','static_p_viol','static_viol','static_p_viol']], linewidth=3)
	ax.set_ylabel("Constraint Violation", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.set_yscale('symlog', linthresh=1e-16)	# Needed for having a log scale showing 0
	ax.set_ylim(bottom=-1e-16)
	#ax.legend(['Townsend', 'Townsend Parallel', 'Gomez-Levy', 'Gomez-Levy Parallel'])

	# set axis where to put the title
	ax = axes[0]


# ====== Plot for the time_numparticles.csv ======
# TODO: update the plot to show the speedup
elif sys.argv[1] == "time_numparticles.csv":
	ax = sns.lineplot(data=data.Serial_time, color='blue', linewidth=3)
	sns.lineplot(data=data.Parallel_time, color='red', ax=ax, linewidth=3)
	ax2 = ax.twinx()
	sns.lineplot(data=data.Speedup, ax=ax2, color='green', linewidth=3)
	ax.legend(handles=[Line2D([], [], marker='_', color='blue', label='Serial'),
					   Line2D([], [], marker='_', color='red',  label='Parallel'),
					   Line2D([], [], marker='_', color='green',  label='Speedup')],)
	ax.set_ylabel("Time (ms)", fontsize=16)
	ax2.set_ylabel("Speedup", fontsize=16)
	ax.set_xlabel("Number of particles", fontsize=16)


# ====== If the csv is not recognized ======
else:
	print("File not supported.")
	exit()

# Set description below title
ax.text(x=0.5, y=1.1, s=title, fontsize=18, weight='bold', ha='center', va='bottom', transform=ax.transAxes)
ax.text(x=0.5, y=1.05, s=description, fontsize=14, alpha=0.75, ha='center', va='bottom', transform=ax.transAxes)


# Save the plot in the output directory
plt.savefig(os.path.join(os.getcwd(), "..", "output", sys.argv[1][:-4] + ".png"))