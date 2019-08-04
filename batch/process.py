import os #, sys
import statistics as st
from datetime import datetime
import csv
import numpy as np
import matplotlib.pyplot as plt

F_SIZE = 16

plt.rc('font', size=F_SIZE)             # controls default text sizes
plt.rc('axes', titlesize=F_SIZE)        # fontsize of the axes title
plt.rc('axes', labelsize=F_SIZE)        # fontsize of the x and y labels
plt.rc('xtick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('ytick', labelsize=F_SIZE)       # fontsize of the tick labels
plt.rc('legend', fontsize=F_SIZE)       # legend fontsize
plt.rc('figure', titlesize=F_SIZE)      # fontsize of the figure title

def to_float_array(input):
  return [float(val) for val in input]

def plotData(target_file_path, x_values, y_names, y_valuess, text_x_axis, text_y_axis):
    xtick = list(range(len(x_values)))

    tmp_colors = ['darkorange', 'green', 'red']

    # ========== Time Plot
    path_result_img = target_file_path + ".pdf"
    fig = plt.figure(figsize=(8, 4.5),frameon=False)
    ax = fig.gca()
    labels = []
    # plot versions
    for i in range(len(y_valuess)):
        #ax.plot(xtick, np.array(y_valuess[i]), 'x-', linewidth=2, color=tmp_colors[i])
        ax.plot(xtick, np.array(y_valuess[i]))
        labels.append(y_names[i])
    plt.xticks(xtick, x_values)
    ax.minorticks_on()
    ax.legend(labels, fancybox=True, shadow=False)                
    ax.grid(b=True, which='major', axis="both", linestyle='-', linewidth=1)
    ax.grid(b=True, which='minor', axis="both", linestyle='-', linewidth=0.4)
    ax.set_xlabel(text_x_axis)
    ax.set_ylabel(text_y_axis)
    #ax.set_title(text_header + " - Execution Time" )
    fig.savefig(path_result_img, dpi=None, facecolor='w', edgecolor='w', 
        format="pdf", transparent=False, bbox_inches='tight', pad_inches=0,
        frameon=False, metadata=None)
    print("saved to "+path_result_img)
    plt.close(fig)

if __name__ == "__main__":
  source_folder = os.path.dirname(os.path.abspath(__file__))
  target_folder_plot  = os.path.join(source_folder, "../../2019_ma_convent_chameleon_task_based_models/Thesis/img")

  # read data
  scaling_values = []
  with open('scaling.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      scaling_values.append(row)
  scaling_values = [*zip(*scaling_values)]
  #print(scaling_values)

  target_file = os.path.join(target_folder_plot, "scaling")
  y_names = [row[0] for row in scaling_values[1:]]
  y_valuess = [row[1:] for row in scaling_values[1:]]
  y_valuess = [to_float_array(y_values) for y_values in y_valuess]
  plotData(target_file, scaling_values[0][1:], y_names, y_valuess, "#Nodes", "Walltime")



    # plot results
#    tmp_target_file_name = "plot_granularity_" + str(int(gran))
#    tmp_target_file_path = os.path.join(target_folder_plot, "imbalance")
#    node_counts = [1, 2, 4, 8]
#    plotData(tmp_target_file_path, unique_n_threads, arr_types, arr_time_original, y_valuess, arr_n_tasks_remote, title_prefix + " - Granularity " + str(int(gran)), "# threads per rank")
#    
#    plotDataMinMaxAvg(tmp_target_file_path + "_bytes_send", unique_n_threads, arr_types, arr_bytes_send_per_msg_min, arr_bytes_send_per_msg_max, arr_bytes_send_per_msg_avg, title_prefix + " - Granularity " + str(int(gran)) + " - Bytes Send", "# threads per rank", "Data Volume [KB]", True, divisor=1000)
#    plotDataMinMaxAvg(tmp_target_file_path + "_bytes_recv", unique_n_threads, arr_types, arr_bytes_recv_per_msg_min, arr_bytes_recv_per_msg_max, arr_bytes_recv_per_msg_avg, title_prefix + " - Granularity " + str(int(gran)) + " - Bytes Recv", "# threads per rank", "Data Volume [KB]", True, divisor=1000)
#    plotDataMinMaxAvg(tmp_target_file_path + "_throughput_send", unique_n_threads, arr_types, arr_throughput_send_min, arr_throughput_send_max, arr_throughput_send_avg, title_prefix + " - Granularity " + str(int(gran)) + " - Throughput Send", "# threads per rank", "Throughput [MB/s]")
#    plotDataMinMaxAvg(tmp_target_file_path + "_throughput_recv", unique_n_threads, arr_types, arr_throughput_recv_min, arr_throughput_recv_max, arr_throughput_recv_avg, title_prefix + " - Granularity " + str(int(gran)) + " - Throughput Recv", "# threads per rank", "Throughput [MB/s]")
