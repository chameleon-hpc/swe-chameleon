import os #, sys
import statistics as st
from datetime import datetime
import csv
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

NUM_EXECUTIONS = 5

F_SIZE = 13

plt.rc('font', size=F_SIZE)             # controls default text sizes
plt.rc('axes', titlesize=F_SIZE)        # fontsize of the axes title
plt.rc('axes', labelsize=F_SIZE)        # fontsize of the x and y labels
plt.rc('xtick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('ytick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('legend', fontsize=F_SIZE)       # legend fontsize
plt.rc('figure', titlesize=F_SIZE)      # fontsize of the figure title

def remove_duplicates_sorted(input):
  #print(input)
  last_element = None
  output = []
  for element in input:
    if element != last_element:
      output.append(element)
    last_element = element
  #print(output)
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

    xtick = np.array(list(range(len(x_values))))

    tmp_colors = ['darkorange', 'green', 'red', 'blue']

    # ========== Time Plot
    path_result_img = target_file_path + ".pdf"
    fig = plt.figure(figsize=(8, 6.5),frameon=False)
    labels = []
    plt.subplot(2,1,1)
    plt.grid(b=True, which='major', axis="both", linestyle='-', linewidth=1, zorder=0)
    plt.grid(b=True, which='minor', axis="both", linestyle='-', linewidth=0.4, zorder=0)
    ax = fig.gca()
    # plot versions
    for i in range(len(y_valuess_medians)):
        #ax.plot(xtick, np.array(y_valuess[i]), 'x-', linewidth=2, color=tmp_colors[i])
        #plt.errorbar(xtick, np.array(y_valuess_medians[i]), [y_valuess_mins[i], y_valuess_maxs[i]], elinewidth=1, capsize=2, marker=".")
        plt.bar(xtick+(i*0.2), np.array(y_valuess_medians[i]), 0.2, zorder=3)
        labels.append(y_names[i])
    plt.xticks(xtick+0.2, x_values)
    plt.minorticks_on()
    #plt.legend(labels,  shadow=False)                
    ax.legend(labels, fancybox=True, loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2)
    ax.set_xlabel(text_x_axis)
    ax.set_ylabel(text_y_axis)
    _,current_top_ylim=ax.get_ylim()
    ax.set_ylim([0, current_top_ylim*1.1])

    #speedups
    y_valuess_medians_tmp = deepcopy(y_valuess_medians)
    for i in range(len(y_valuess_medians_tmp[0])):
      for j in range(len(y_valuess_medians_tmp)):
        y_valuess_medians_tmp[j][i] = y_valuess_medians[0][i] / y_valuess_medians[j][i]
    y_valuess_medians = y_valuess_medians_tmp
    plt.subplot(2,1,2)
    ax = fig.gca()
    # plot versions
    for i in range(len(y_valuess_medians)):
        plt.plot(xtick, np.array(y_valuess_medians[i]), '.-', linewidth=2)
        #plt.errorbar(xtick, np.array(y_valuess_medians[i]), [y_valuess_mins[i], y_valuess_maxs[i]], elinewidth=1, capsize=2, marker=".")
        labels.append(y_names[i])
    plt.xticks(xtick, x_values)
    plt.minorticks_on()
    ax.tick_params(axis='x', which='minor', bottom=False)
    #plt.legend(labels, fancybox=True, shadow=False)                
    plt.grid(b=True, which='major', axis="both", linestyle='-', linewidth=1)
    plt.grid(b=True, which='minor', axis="both", linestyle='-', linewidth=0.4)
    ax.set_xlabel(text_x_axis)
    ax.set_ylabel("Speedup compared\nto MPI+OpenMP")
    _,current_top_ylim=ax.get_ylim()
    ax.set_ylim([0, current_top_ylim*1.1])
    
    plt.subplots_adjust(hspace=0.75)

    #plt.set_title(text_header + " - Execution Time" )
    fig.savefig(path_result_img, dpi=None, facecolor='w', edgecolor='w', 
        format="pdf", transparent=False, bbox_inches='tight', pad_inches=0,
        frameon=False, metadata=None)
    print("saved to "+path_result_img)
    plt.close(fig)

def plotDataBarsOnly(target_file_path, x_values, y_names, y_valuess, text_x_axis, text_y_axis):

    # calc medians, mins and maxs
    y_valuess_medians = []
    y_valuess_mins = []
    y_valuess_maxs = []
    for y_values in y_valuess:
      medians, mins, maxs = to_medians_mins_maxs(y_values, int(len(y_values)/NUM_EXECUTIONS))
      y_valuess_medians.append(medians)
      y_valuess_mins.append(mins)
      y_valuess_maxs.append(maxs)

    xtick = np.array(list(range(len(x_values))))

    tmp_colors = ['darkorange', 'green', 'red', 'blue']

    # ========== Time Plot
    path_result_img = target_file_path + ".pdf"
    fig = plt.figure(figsize=(8, 3),frameon=False)
    labels = []
    plt.grid(b=True, which='major', axis="both", linestyle='-', linewidth=1, zorder=0)
    plt.grid(b=True, which='minor', axis="both", linestyle='-', linewidth=0.4, zorder=0)
    ax = fig.gca()
    # plot versions
    for i in range(len(y_valuess_medians)):
        #ax.plot(xtick, np.array(y_valuess[i]), 'x-', linewidth=2, color=tmp_colors[i])
        #plt.errorbar(xtick, np.array(y_valuess_medians[i]), [y_valuess_mins[i], y_valuess_maxs[i]], elinewidth=1, capsize=2, marker=".")
        plt.bar(xtick+(i*0.2), np.array(y_valuess_medians[i]), 0.2, zorder=3)
        labels.append(y_names[i])
    plt.xticks(xtick+0.2, x_values)
    plt.minorticks_on()
    #plt.legend(labels,  shadow=False)                
    ax.legend(labels, fancybox=True, loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2)
    ax.set_xlabel(text_x_axis)
    ax.set_ylabel(text_y_axis)
    _,current_top_ylim=ax.get_ylim()
    ax.set_ylim([0, current_top_ylim*1.1])

    #plt.set_title(text_header + " - Execution Time" )
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

#  target_file = os.path.join(target_folder_plot, "scaling_bar")
#  plotDataBar(target_file, x_values, y_names, y_valuess, "#Nodes", "Walltime [s]")

  # speedups
  target_file = os.path.join(target_folder_plot, "scaling_speedups")
  y_valuess_tmp = deepcopy(y_valuess)
  for i in range(len(y_valuess_tmp[0])):
    for j in range(len(y_valuess_tmp)):
      y_valuess_tmp[j][i] = y_valuess[0][i] / y_valuess[j][i]
  y_valuess = y_valuess_tmp

  plotData(target_file, x_values, y_names, y_valuess, "#Nodes", "Speedup over MPI+OpenMP")

  # read interference data
  scaling_values = []
  with open('interference.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)

  target_file = os.path.join(target_folder_plot, "interference")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  x_values = remove_duplicates_sorted(scaling_values[0][1:])

  plotData(target_file, x_values, y_names, y_valuess, "Interference Fraction", "Walltime [s]")

  # read single interference data
  scaling_values = []
  with open('single_interference.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)

  target_file = os.path.join(target_folder_plot, "single_interference")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  x_values = remove_duplicates_sorted(scaling_values[0][1:])
  x_values = [x.replace("1.", "0.") for x in x_values]

  plotData(target_file, x_values, y_names, y_valuess, "Interference Fraction (Single Thread)", "Walltime [s]")

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

  # imbalance speedups
  target_file = os.path.join(target_folder_plot, "imbalance_speedups")
  y_valuess_tmp = deepcopy(y_valuess)
  for i in range(len(y_valuess_tmp[0])):
    for j in range(len(y_valuess_tmp)):
      y_valuess_tmp[j][i] = y_valuess[0][i] / y_valuess[j][i]
  y_valuess = y_valuess_tmp

  plotData(target_file, x_values, y_names, y_valuess, "Dry Fraction", "Speedup over MPI+OpenMP")

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

  plotDataBarsOnly(target_file, x_values, y_names, y_valuess, "Blocksize", "Walltime [s]")
