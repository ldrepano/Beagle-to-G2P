
'''
Beagle to G2P Pipeline
October 2024 
Laura Drepanos and Ganna Reint


Purpose: 

Example usage: python3 BE_G2P_output_visualization.py G2P_PPM1D_O15297_protein_features_CBE.csv


Outputs: html file to visualize G2P information layered with BE screen results 


'''

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import stats
import plotly.io as pio
import plotly.graph_objs as go
from itertools import cycle
import matplotlib.colors

import plotly.express as px



def main():
	if len(sys.argv)<2:
		print("\nPlease provide filepath for G2P output to visualize as a command line argument \n")
		exit()

	inputfile=sys.argv[1]
	g2p_output = pd.read_csv(inputfile, delimiter=',') 
	#remove residues not targeted by any z-score
	g2p_output = g2p_output.dropna(subset=['Strongest Z-score'])

	traces=[]
	potential_feature_cols=[col for col in g2p_output.columns if col not in ['residueId','Strongest Z-score','Strongest Z-Score mutation','Strongest Z-Score sgRNA','Ref AA','AA','Strongest Z-Score Alt AA','Sequence','Structure', 'Mean Z-Score', 'Confidence','Most Severe Possible Mutation',"Chain (UniProt)"]]
	
	#mutation type determines point shape
	mutation_shape_dict={"Missense":"square","Nonsense":"x","Silent":"circle"}
	g2p_output["Mutation Shape"]=g2p_output['Strongest Z-Score mutation'].map(mutation_shape_dict)
	

	#Just adding the shapes for different mutations to the legend 
	traces.append(go.Scatter(x=[None], y=[None], mode='markers',
	                       marker=dict(size=10, symbol="square",color="white",line=dict(color='grey',width=2)),
	                       legendgroup='Marker Shape', legendgrouptitle_text="Marker Shape",showlegend=True, name='Missense Mutations'))
	traces.append(go.Scatter(x=[None], y=[None], mode='markers',
	                       marker=dict(size=10, symbol="circle",color="white",line=dict(color='grey',width=2)),
	                       legendgroup='Marker Shape', showlegend=True, name='Silent Mutations'))
	traces.append(go.Scatter(x=[None], y=[None], mode='markers',
	                       marker=dict(size=10, symbol="x",color="white",line=dict(color='grey',width=2)),
	                       legendgroup='Marker Shape', showlegend=True, name='Nonsense Mutations'))

	color_list = ["pink","red","blue","yellow",'orange',"purple","green","tan","grey","black"]
	for feature,feature_color in zip(potential_feature_cols,color_list):
		g2p_output_feature= g2p_output.dropna(subset=[feature])
		new_trace=go.Scatter(x=g2p_output_feature['residueId'], y=g2p_output_feature['Strongest Z-score'], 
			mode='markers',name=feature, marker=dict(size=20, color=feature_color, opacity=0.9,symbol=g2p_output_feature['Mutation Shape']),hovertemplate=g2p_output_feature[feature])
		new_trace.hoverinfo="text"
		new_trace.text=feature 
		traces.append(new_trace)
	layout = go.Layout(
	    xaxis=dict(title='Residue', zeroline=True, zerolinewidth=1, zerolinecolor='black', color='black', showgrid=True, gridcolor='lightgray'),
	    yaxis=dict(title='Z-Score', zeroline=True, zerolinewidth=1, zerolinecolor='black', color='black', showgrid=True, gridcolor='lightgray'),
	    legend=dict(x=1.05, y=1, bgcolor='white', borderwidth=2, bordercolor='black',title="Target Info (select features below)"),
	    width=1600,  # Set width of the plot
	    height=600,  # Set height of the plot,
	    plot_bgcolor='white'  # Set background color to white

	)
	# Create figure
	fig = go.Figure(data=traces, layout=layout)

	pio.write_html(fig, file=f"{inputfile[:-4]}.html")



if __name__ == '__main__':
    main()


