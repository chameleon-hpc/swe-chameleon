import os #, sys
import statistics as st
from datetime import datetime
import csv
import numpy as np
import matplotlib.pyplot as plt

NUM_EXECUTIONS = 5

F_SIZE = 16

plt.rc('font', size=F_SIZE)             # controls default text sizes
plt.rc('axes', titlesize=F_SIZE)        # fontsize of the axes title
plt.rc('axes', labelsize=F_SIZE)        # fontsize of the x and y labels
plt.rc('xtick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('ytick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('legend', fontsize=F_SIZE)       # legend fontsize
plt.rc('figure', titlesize=F_SIZE)      # fontsize of the figure title

def remove_duplicates_sorted(input):
  print(input)
  last_element = None
  output = []
  for element in input:
    if element != last_element:
      output.append(element)
    last_element = element
  print(output)
  return output

def to_float_array(input):
  return [float(val) for val in input]

def to_medians_mins_maxs(input, length):
  medians = []; mins = []; maxs = [];
  for i in range(length):
    values = input[i*NUM_EXECUTIONS:(i+1)*NUM_EXECUTIONS]
    medians.append(st.median(values))
    mins.append(abs(min(values) - st.median(values)))
    maxs.append(abs(max(values) - st.median(values)))
  return medians, mins, maxs

def plotData(target_file_path, x_values, y_names, y_valuess, text_x_axis, text_y_axis):

    # calc medians, mins and maxs
    y_valuess_medians = []
    y_valuess_mins = []
    y_valuess_maxs = []
    for y_values in y_valuess:
      medians, mins, maxs = to_medians_mins_maxs(y_values, int(len(y_values)/NUM_EXECUTIONS))
      y_valuess_medians.append(medians)
      y_valuess_mins.append(mins)
      y_valuess_maxs.append(maxs)

    xtick = list(range(len(x_values)))

    tmp_colors = ['darkorange', 'green', 'red']

    # ========== Time Plot
    path_result_img = target_file_path + ".pdf"
    fig = plt.figure(figsize=(8, 4.5),frameon=False)
    ax = fig.gca()
    labels = []
    # plot versions
    for i in range(len(y_valuess_medians)):
        #ax.plot(xtick, np.array(y_valuess[i]), 'x-', linewidth=2, color=tmp_colors[i])
        ax.errorbar(xtick, np.array(y_valuess_medians[i]), [y_valuess_mins[i], y_valuess_maxs[i]], elinewidth=1, capsize=2)
        labels.append(y_names[i])
    plt.xticks(xtick, x_values)
    ax.minorticks_on()
    ax.legend(labels, fancybox=True, shadow=False)                
    ax.grid(b=True, which='major', axis="both", linestyle='-', linewidth=1)
    ax.grid(b=True, which='minor', axis="both", linestyle='-', linewidth=0.4)
    ax.set_xlabel(text_x_axis)
    ax.set_ylabel(text_y_axis)
    ax.set_ylim(ymin=0)
    #ax.set_title(text_header + " - Execution Time" )
    fig.savefig(path_result_img, dpi=None, facecolor='w', edgecolor='w', 
        format="pdf", transparent=False, bbox_inches='tight', pad_inches=0,
        frameon=False, metadata=None)
    print("saved to "+path_result_img)
    plt.close(fig)

if __name__ == "__main__":
  source_folder = os.path.dirname(os.path.abspath(__file__))
  target_folder_plot  = os.path.join(source_folder, "../../2019_ma_convent_chameleon_task_based_models/Thesis/img")

  # read scaling data
  scaling_values = []
  with open('scaling.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)

  target_file = os.path.join(target_folder_plot, "scaling")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  x_values = remove_duplicates_sorted(scaling_values[0][1:])

  plotData(target_file, x_values, y_names, y_valuess, "#Nodes", "Walltime [s]")

  # read imbalance data
  scaling_values = []
  with open('imbalance.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)

  target_file = os.path.join(target_folder_plot, "imbalance")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  x_values = remove_duplicates_sorted(scaling_values[0][1:])

  plotData(target_file, x_values, y_names, y_valuess, "Dry Fraction", "Walltime [s]")

  # read granularity data
  scaling_values = []
  with open('granularity.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)

  target_file = os.path.join(target_folder_plot, "granularity")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  x_values = remove_duplicates_sorted(scaling_values[0][1:])

  plotData(target_file, x_values, y_names, y_valuess, "#Blocks", "Walltime [s]")
