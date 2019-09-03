#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 15:56:15 2019

@author: bh
"""
import pandas as pd
import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool 
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Button
from bokeh.layouts import column
from bokeh.io import curdoc

output_file('histograma_massa_invariante.html', title="Invariant mass histogram")

#leitura do arquivo de dados
df = pd.DataFrame()
df = pd.read_csv('/home/bh/Documentos/DoubleMuRun2011A.csv', engine='python') #Lê o arquivo com os dados.


btn_bin = Button(label="REBIN", button_type="success")

#add hover
hover = HoverTool(tooltips=[('Counts', '@inv_mass'),('Invariant mass (GeV)', '$x')])
hover.point_policy = "follow_mouse"
bins = 500

def rebin(newbin):
    

    weights = []
    for i in df["M"]:
        weights.append(newbin/np.log(10)/i)
  
    
    #create histogram    
    hist, edges = np.histogram(df["M"], bins=np.logspace(-0.05, 2.50, num=newbin+1), weights=weights)

    # Put the information in a dataframe
    inv_mass = pd.DataFrame({'inv_mass': hist, 
                             'left': edges[:-1], 
                             'right': edges[1:]})

    # Convert dataframe to column data source
    inv_mass_plotable = ColumnDataSource(inv_mass)
    

    p = figure(title='Dimuon invariant mass spectrum \n', tools=[hover, 'reset','zoom_in', 'zoom_out', 'pan', 'save', 'box_zoom'],
               background_fill_color="white", plot_height = 600, plot_width = 1350, y_range=[10**0, max(hist)],x_axis_type="log",
               y_axis_type="log")
    p.quad(source=inv_mass_plotable, top='inv_mass', bottom=1, left='left', right='right',
           fill_color="blue", line_color="black", hover_fill_color='yellow', hover_line_color="red")
    #p.y_range.start = 0
    p.xaxis.axis_label = 'Invariant mass (GeV)'
    p.yaxis.axis_label = 'Counts'
    p.grid.grid_line_color="grey"
    
    return p

dimuon = rebin(bins)

btn_bin.on_click(rebin(100))

curdoc().add_root(column(dimuon, btn_bin))