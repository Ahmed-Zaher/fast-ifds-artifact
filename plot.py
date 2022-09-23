import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rc('axes', titlesize=40)
plt.rc('axes', labelsize=40)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
xticksthing = np.linspace(0, 500000, 2000000)
def plotAvgQueryTime(df, analysis, b, anchor_x = 0, anchor_y = 1):
  fig = plt.figure(figsize=(15, 8), dpi=120)
  ax = fig.add_subplot()

  ax.scatter(df['m_GExp'], (df['preprocessing_time'] + df['all_query_time']) / df['n_G'], s=200, c='#377eb8', marker="s", label='PARAM')
  if b == "IFDS":
    ax.scatter(df['m_GExp'], df['old_avg_query_time'], s=400, c='#ff7f00', marker="^", label='IFDS')
  else:
    ax.scatter(df['m_GExp'], df['on_demand_avg_query_time'], s=200, c='#4daf4a', marker="o", label='DEM')

  ax.set_xlabel("Instance size")
  ax.set_ylabel("Average query time (s)")
  title = "Reachability"
  if (analysis == "null_ptr"):
    title = "Null-pointer"
  if (analysis == "uninit_var"):
    title = "Uninitialized variables"
  ax.set_title(title)
  ax.ticklabel_format(useOffset=False, style='plain')
  ax.legend(loc='upper left', bbox_to_anchor=(anchor_x, anchor_y), prop={'size': 35})
  #plt.grid()
  plt.savefig("final_output/" + analysis + "_" + b + "_avg_query_time.pdf")

def plotPerAnalysis(analysis, b):
  df = pd.read_csv('result.csv')
  df = df.rename(columns={c: c.strip() for c in df.columns})
  df = df[(df['analysis'] == analysis)]
  df.head()
  plotAvgQueryTime(df, analysis, b)
plotPerAnalysis("reachability", "IFDS")
plotPerAnalysis("null_ptr", "IFDS")
plotPerAnalysis("uninit_var", "IFDS")
plotPerAnalysis("reachability", "DEM")
plotPerAnalysis("null_ptr", "DEM")
plotPerAnalysis("uninit_var", "DEM")
