from shiny import App, render, ui
import numpy as np
import matplotlib.pyplot as plt

app_ui = ui.page_fluid(
    ui.input_slider("bins", "Number of bins:", min=5, max=50, value=20),
    ui.output_plot("histogram")
)

def server(input, output, session):
    @output
    @render.plot
    def histogram():
        data = np.random.randn(1000)
        plt.hist(data, bins=input.bins(), color='blue', edgecolor='black')
        plt.title("Histogram")

app = App(app_ui, server)