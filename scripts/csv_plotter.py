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

# Set the font size
sns.set_theme(font_scale = 1.3)

# Compose the title from the comments
title = ""
description = "("
with open(filepath, 'r') as f:
	title = f.readline()[1:].strip().replace('_', ' ').capitalize()
	for line in f:
		if line[0] == '#':
			description += line[1:].strip() + ", "
	description = description[:-2] + ")"


sns.set_style("darkgrid")
ax = None


# ====== Plot for the saspso_static_adaptive.csv ======
if sys.argv[1] == "saspso_static_adaptive.csv":
	fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

	# error plot
	ax = sns.lineplot(ax=axes[0], data=data[['static_err','static_p_err','adaptive_err','adaptive_p_err']], linewidth=3, legend=True)
	ax.set_yscale('symlog', linthresh=1e-16)	# Needed for having a log scale showing 0
	ax.set_ylim(bottom=-1e-16)
	ax.set_xlim(right=1400)
	ax.set_ylabel("Error", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.legend(fontsize=14)

	# constraint violation plot
	ax = sns.lineplot(ax=axes[1], data=data[['static_viol','static_p_viol','adaptive_viol','adaptive_p_viol']], linewidth=3)
	ax.set_ylabel("Constraint Violation", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.set_yscale('symlog', linthresh=1e-16)	# Needed for having a log scale showing 0
	ax.set_ylim(bottom=-1e-16)
	ax.legend(fontsize=14)

	# set axis where to put the title
	ax = axes[0]


# ====== Plot for the saspso_time_numparticles.csv ======
elif sys.argv[1] == "saspso_time_numparticles.csv":
	plt.figure(figsize=(12, 10))
	ax = sns.lineplot(data=data.Serial_time, color='blue', linewidth=3)
	sns.lineplot(data=data.Parallel_time, color='red', ax=ax, linewidth=3)
	ax2 = ax.twinx()
	sns.lineplot(data=data.Speedup, ax=ax2, color='green', linewidth=3)
	ax.legend(handles=[Line2D([], [], marker='_', color='blue', label='Serial'),
					   Line2D([], [], marker='_', color='red',  label='Parallel'),
					   Line2D([], [], marker='_', color='green',  label='Speedup')],
			fontsize=14)
	ax.set_ylabel("Time (ms)", fontsize=16)
	ax2.set_ylabel("Speedup", fontsize=16)
	ax.set_xlabel("Number of particles", fontsize=16)


# ====== Plot for the saspso_optimize.csv ======
elif sys.argv[1] == "saspso_optimize.csv":
	fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

	# value plot
	ax = sns.lineplot(ax=axes[0], data=data[['value']], linewidth=3, legend=True)
	ax.set_yscale('linear')
	#ax.set_ylim(bottom=-1e-16)
	ax.set_ylabel("Value", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.axhline(y=8150)

	# constraint violation plot
	ax = sns.lineplot(ax=axes[1], data=data[['violation']], linewidth=3, legend=False)
	ax = sns.lineplot(ax=axes[1], data=data[['threshold']], linewidth=3, palette=['r',], alpha=0.5, legend=False)
	ax.set_ylabel("Constraint Violation", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.set_yscale('symlog', linthresh=1e-3)
	ax.set_ylim(bottom=-1e-3, top=10)
	ax.set_xlim(left=-100, right=14000)

	# feasible particles plot as twin of the constraint violation plot
	ax = sns.lineplot(ax=axes[1].twinx(), data=data[['feasible_particles']], linewidth=3, palette=['g',], legend=False)
	ax.set_ylabel("Feasible Particles", fontsize=16)
	ax.set_xlabel("Iterations", fontsize=16)
	ax.set_yscale('symlog', linthresh=50)

	ax.legend(handles=[Line2D([], [], marker='_', color='blue', label='Total Violation'),
					   Line2D([], [], marker='_', color='red',  label='Violation Threshold'),
					   Line2D([], [], marker='_', color='green',  label='Feasible Particles')],
					   fontsize=14,
					   loc='upper center')

	# set axis where to put the title
	ax = axes[0]


# ====== Plot for the abc_optimize.csv ======
elif sys.argv[1] == "abc_optimize.csv":
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

    # Value plot
    ax = sns.lineplot(ax=axes[0], data=data[['value']], linewidth=3, legend=True)
    ax.set_yscale('linear')
    ax.set_ylabel("Value", fontsize=16)
    ax.set_xlabel("Iterations", fontsize=16)
    ax.axhline(y=8150)

    # Constraint violation plot
    ax = sns.lineplot(ax=axes[1], data=data[['violation']], linewidth=3, legend=False)
    ax.set_ylabel("Constraint Violation", fontsize=16)
    ax.set_xlabel("Iterations", fontsize=16)
    ax.set_yscale('symlog', linthresh=1e-3)
    ax.set_ylim(bottom=-1e-3, top=10)
    ax.set_xlim(left=-100, right=6000)

    # Feasible particles plot as twin of the constraint violation plot
    ax = sns.lineplot(ax=axes[1].twinx(), data=data[['feasible_particles']], linewidth=3, palette=['g',], legend=False)
    ax.set_ylabel("Feasible Particles", fontsize=16)
    ax.set_xlabel("Iterations", fontsize=16)
    ax.set_yscale('symlog', linthresh=50)

    ax.legend(handles=[Line2D([], [], marker='_', color='blue', label='Total Violation'),
                       Line2D([], [], marker='_', color='green',  label='Feasible Particles')],
                       fontsize=14,
                       loc='upper center')

    # Set axis where to put the title
    ax = axes[0]


# ====== If the csv is not recognized ======
else:
	print("File not supported.")
	exit()

# Set description below title
ax.text(x=0.5, y=1.1, s=title, fontsize=18, weight='bold', ha='center', va='bottom', transform=ax.transAxes)
ax.text(x=0.5, y=1.05, s=description, fontsize=14, alpha=0.75, ha='center', va='bottom', transform=ax.transAxes)


# Save the plot in the output directory
plt.savefig(os.path.join(os.getcwd(), "..", "output", sys.argv[1][:-4] + ".png"))
