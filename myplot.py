''' Class for convenient and pretty plotting '''

import matplotlib
import matplotlib.pyplot as plt
import sys
# import numpy as np
import pandas as pd

save = 'save' in sys.argv

class Plot:
    def __init__(self,figsize,n=1,sharex=False):
        '''
        figsize [(w,h)]: tuple with figure width x height
        n [int|string]: number of subplots (int) or type of plot ('paramspace')

        '''
        self.n = n # need later for tight layout

        if n == 'paramspace':
            # definitions for the axes
            left, width = 0.13, 0.65
            bottom, height = 0.1, 0.65
            spacing = 0.005

            rect_scatter = [left, bottom, width, height]
            rect_histx = [left, bottom + height + spacing, width, 0.2]
            rect_histy = [left + width + spacing, bottom, 0.2, height]

            # start with a rectangular Figure
            self.fig = plt.figure(figsize=(8, 8))
            # main axis (square plot in the middle)
            self.axes = [plt.axes(rect_scatter)]
            self.axes[0].tick_params(direction='in', top=True, right=True)
            # plot on top of the main one
            self.axes.append(plt.axes(rect_histx))
            self.axes[1].tick_params(direction='in', labelbottom=False)
            # plot at the side
            self.axes.append(plt.axes(rect_histy))
            self.axes[2].tick_params(direction='in', labelleft=False)
        else:
            fig, ax = plt.subplots(n,figsize=figsize, sharex=sharex)
            self.fig = fig
            # convenient if only have one subplot (not ty type axes[0] all the time)
            self.ax = ax
            # array of all our axes; new axes get added here too (e.g. two subplots, one has 2 y axes, means in total 5 axes in the list)
            self.axes = list(ax) if n > 1 else [ax]


    def legend(self, out=False, ncol=1, title=None, pos=None):
        '''
        Create a customized legend.

        out [True|False]: legend outside of the plot in the bottom (default False)
        ncol [int]: number of columns in the legend (default 1)
        pos [list of int|string]: position of the legend (default "best").
            Important for plotting with a second Y axis (see below)
        '''

        if out:
            # draw legend at the bottom
            self.axes[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),\
                fancybox=True, shadow=False, ncol=ncol,fontsize=15)
            # remove legends from all other subplots
            for ax in self.axes[:-1]:
                legend = ax.legend()
                legend.remove()
        else:
            # draw legend on each axis
            # in case these are not subplots, but one plot with two y axes, legends might crash
            # -> choose given location
            for i in range(len(self.axes)):
                # loc = pos[i] if pos else None
                if pos != None: loc = pos if type(pos) == int else pos[i]
                else: loc = None

                # loc = pos if type(pos) pos[i] if pos else None
                self.axes[i].legend(ncol=ncol, loc=loc, title=title)


    def pretty(self, large=3, stretch = None, grid='major'):
        '''
        Make the plot pretty.

        large [int]: controls fontsizes (default 3)
        stretch [None|"float"|"year"] change x axis range (default None).
            Useful for cases when the automatic ranges makes the leftmost and
        rightmost point fall right at the frame and get "eaten up"
            None: do not do anything
            "float": increase the range by 10%
            "year": increase the range by 1

        grid ["major"|"minor"]: only major grid, or also minor;
            same as the grid argument of the function ax.grid()
        '''

        for ax in self.axes:
            ## make better x and y limits
            if stretch:
                lims = {'x': list(ax.get_xlim()), 'y': list(ax.get_ylim())}
                # xlim = list(ax.get_xlim())
                # ylim = list(ax.get_ylim())

                strength = 0.005

                if stretch == 'float':
                    for axx in lims:
                        for i in range(2):
                            if lims[axx][i] > 0: lims[axx][i] = lims[axx][i] * (1 + (-1)**(i+1)*strength)
                            elif lims[axx][i] < 0: lims[axx][i] = lims[axx][i] * (1 + (-1)**i*strength)
                            else: lims[axx][i] = (-1)**(i+1) * 0.01
                    # xlim_orig = xlim[0]
                    # xlim[0] = - xlim[1] * strength if xlim[0] == 0 \
                    #     else xlim[0] *(1 - strength)
                    # xlim[1] = xlim[1] * (1 - strength) if xlim_orig == 0  else xlim[1] + xlim_orig * strength

                # elif stretch == 'year':
                #     xlim[0] = xlim[0] - 1
                #     xlim[1] = xlim[1] + 1


                # ylim[0] = - ylim[1]*0.02 if ylim[0] == 0 else ylim[0]*0.98

                # ax.set_xlim(xlim[0], xlim[1])
                # ax.set_ylim(ylim[0], ylim[1]*1.02)

                ax.set_xlim(lims['x'][0], lims['x'][1])
                ax.set_ylim(lims['y'][0], lims['y'][1])

            ## increase tick sizes (numbers)
            for t in ax.get_xaxis().get_ticklabels():
                t.set_fontsize(15 + large)
            for t in ax.get_yaxis().get_ticklabels():
                t.set_fontsize(15 + large)

            ## increase title sizes
            ax.get_xaxis().get_label().set_fontsize(17 + large)
            ax.get_yaxis().get_label().set_fontsize(17 + large)
            ax.title.set_fontsize(20)

            ## increase legend font
            legend = ax.get_legend()
            if legend:
                if legend.get_title(): legend.get_title().set_fontsize(19 + large)
                for t in legend.get_texts():
                    t.set_fontsize(17 + large)

            ## add gridlines
            if grid: ax.grid(linestyle='--',zorder=0, which=grid)

        ## increase fontsize of texts on the figure
        for t in self.fig.get_children():
            if type(t) == matplotlib.text.Text:
                t.set_fontsize(20)



        ## remove white borders around the plot
        # self.fig.tight_layout()
        if type(self.n) == int:
            self.fig.tight_layout(rect=[0,0,1,0.97])


    def figure(self, name=None):
        '''
        Show the plot or save an image.

        name [string]: name of the figure to be saved (with extension)
        '''
        #self.pretty(2)
        print ('Image:', name)

        if save:
            if name:
                self.fig.savefig(name)
                print ('(saved)')
            else:
                print ("You didn't provide an image name!")
        else:
            print ('(NOT SAVED)')
            plt.show()


    def add_axis(self, col=None, ax=0):
        '''
        Add a second Y axis.

        col [string]: colour of the axis. Can be useful to differentiate
            the curves and see which one belongs to which Y axis
        ax [int]: a second axis is created as a twin of the given
            existing axis (default 0).
        '''
        ax2 = self.axes[ax].twinx()
        if col: ax2.tick_params(axis='y', colors=col)
        self.axes.append(ax2)

    def __del__(self):
        plt.close(self.fig)
