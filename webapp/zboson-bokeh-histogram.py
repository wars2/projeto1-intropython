#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 15:56:15 2019

@author: bh
"""
import pandas as pd
import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, ColumnDataSource, Panel, Tabs
from bokeh.models.widgets import Button, Dropdown
from bokeh.layouts import column

output_file('histograma_massa_invariante_nonserver.html', title="Invariant mass histogram")

#leitura do arquivo de dados
df = pd.DataFrame()
df = pd.read_csv('/home/bh/Documentos/DoubleMuRun2011A.csv', engine='python') #Lê o arquivo com os dados.


btn_bin = Button(label="calcular fit")

options = [("Breit-Wigner", "Breit-Wigner"), ("Duas Gaussianas + Exponencial", "Duas Gaussianas + Exponencial"), ("CrystalBall + Exponencial", "CrystalBall + Exponencial")]
dropdown = Dropdown(label="fit selection", button_type="warning", menu=options)

#add hover
hover = HoverTool(tooltips=[('Counts', '@inv_mass'),('Invariant mass (GeV)', '$x')])
hover.point_policy = "follow_mouse"
bins = 500

def plot(newbin, begin, end, y_type = 'linear'):
    

    weights = []
    for i in df["M"]:
        weights.append(newbin/np.log(10)/i)
  
    
    #create histogram    
    hist, edges = np.histogram(df["M"], bins=np.logspace(begin, end, num=newbin+1), weights=weights)

    # Put the information in a dataframe
    inv_mass = pd.DataFrame({'inv_mass': hist, 
                             'left': edges[:-1], 
                             'right': edges[1:]})

    # Convert dataframe to column data source
    inv_mass_plotable = ColumnDataSource(inv_mass)
    

    p = figure(tools=[hover, 'reset','zoom_in', 'zoom_out', 'pan', 'save', 'box_zoom'],background_fill_color="white",
               plot_height = 600, plot_width = 1350, y_range=[10**0, max(hist)],x_axis_type="log",y_axis_type=y_type)
    
    p.quad(source=inv_mass_plotable, top='inv_mass', bottom=1, left='left', right='right',
           fill_color="blue", line_color="black", hover_fill_color='yellow', hover_line_color="red")
    #p.y_range.start = 0
    p.xaxis.axis_label = 'Invariant mass (GeV)'
    p.yaxis.axis_label = 'Counts'
    p.grid.grid_line_color="grey"
    
    return p

dimuon = plot(bins, -0.5, 2.5, 'log')
dimuon_tab = Panel(child=dimuon, title='Dimuon invariant mass spectrum')

rho = plot(bins/10, -0.15, -0.05 )
rho_tab = Panel(child=rho, title=u"\u03c1 peak")

phi = plot(bins/10, -0.09, 0.09)
phi_tab = Panel(child=phi, title=u"\u03a6 peak")

jpsi = plot(bins/10, 0.47, 0.52)
jpsi_tab = Panel(child=jpsi, title=u"J/\u03c8 peak")

psi = plot(bins/10, 0.51, 0.62)
psi_tab = Panel(child=psi, title=u"\u03c8' peak")

upsilon = plot(bins/10, 0.92, 1.02)
upsilon_tab = Panel(child=upsilon, title=u"\u03a5 peak")

z = plot(bins/5, 1.85, 2.05)
z_tab = Panel(child=z, title="Z peak")

#btn_bin.on_click(plot(100))

tabs = Tabs(tabs=[ dimuon_tab, rho_tab, phi_tab,jpsi_tab, psi_tab,upsilon_tab, z_tab])

show(column(tabs, btn_bin,dropdown))