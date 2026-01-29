import vtk
from vtk import (vtkUnstructuredGridReader, vtkDataSetMapper, vtkActor,
                 vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor)
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

#def plot_density_contour(vtk_file, scalar_name, contour_values):
#    # Read the VTK file
#    reader = vtk.vtkUnstructuredGridReader()
#    reader.SetFileName(vtk_file)
#    reader.ReadAllScalarsOn()
#    reader.Update()
#
#    # Create a contour filter and set its input
#    contour_filter = vtk.vtkContourFilter()
#    contour_filter.SetInputConnection(reader.GetOutputPort())
#    # Set specific contour values
#    for i, value in enumerate(contour_values):
#        contour_filter.SetValue(i, value)
#
#    # Create a mapper and set the contour filter as its input
#    contour_mapper = vtk.vtkPolyDataMapper()
#    contour_mapper.SetInputConnection(contour_filter.GetOutputPort())
#
#    # Create an actor to display the contours
#    contour_actor = vtk.vtkActor()
#    contour_actor.SetMapper(contour_mapper)
#
#    # Create a renderer, render window, and interactor
#    renderer = vtk.vtkRenderer()
#    render_window = vtk.vtkRenderWindow()
#    render_window.AddRenderer(renderer)
#    render_window_interactor = vtk.vtkRenderWindowInteractor()
#    render_window_interactor.SetRenderWindow(render_window)
#
#    # Add the actor to the scene
#    renderer.AddActor(contour_actor)
#    renderer.SetBackground(1, 1, 1)  # Background color white
#
#    # Render and start interaction
#    render_window.Render()
#    render_window_interactor.Start()
#    # Create the RendererWindowInteractor and display the vtk_file
#    interactor = vtkRenderWindowInteractor()
#    interactor.SetRenderWindow(render_window)
#    interactor.Initialize()
#    interactor.Start()
#
## Replace with the path to your VTK file
#vtk_file = 'soln_c00200.vtk'
#
## Define specific density values for contours
#contour_values = [0.1, 0.2, 0.3, 0.4, 0.5]  # Modify as needed
#
## Plot the density contour
#plot_density_contour(vtk_file, 'density', contour_values)



# read the data
grid = pv.read('soln_c00200.vtk')

# plot the data with an automatically created Plotter
grid.plot(show_scalar_bar=False, show_axes=False)
