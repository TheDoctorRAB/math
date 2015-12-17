########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.16.December.2015
########################################################################
# 
# Binomial distribution
#
########################################################################
#
#
#
########################################################################
#
# imports
#
import numpy
import scipy.stats
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.ticker import MultipleLocator
from win32api import GetSystemMetrics
#
########################################################################
#
#
#
#######
#
# diagnostics
#
matplotlib.rcParams.update({'font.size': 16}) #set plot font
width=GetSystemMetrics (0)
height=GetSystemMetrics (1)
#
# input data
#
domain_x=10 #success domain
trials=10 #total trials
x=numpy.linspace(1,domain_x,domain_x)
probability=0.5
#
#######
#
# compute distribution
#
binomial_pdf=scipy.stats.binom.pmf(x,trials,probability)
bimax=numpy.amax(binomial_pdf)
#
#######
#
# plot
#
fig,left_axis=plot.subplots()
title='Binomial PDF example'
xtitle='x'
ytitle='f(x)'
#
plot.title(title)
left_axis.set_xlabel(xtitle)
left_axis.set_ylabel(ytitle)
#
xmin=0.9
xmax=domain_x+0.1
ymin=0
ymax=bimax+0.05*bimax
#
xmajortick=1
xminortick=0.5
ymajortick=0.5*bimax
yminortick=.025*bimax
#
plot.xlim(xmin,xmax)
left_axis.axis(ymin=ymin,ymax=ymax)
#
left_axis.xaxis.set_major_locator(MultipleLocator(xmajortick))
left_axis.yaxis.set_major_locator(MultipleLocator(ymajortick))
left_axis.xaxis.set_minor_locator(MultipleLocator(xminortick))
left_axis.yaxis.set_minor_locator(MultipleLocator(yminortick))
#
left_axis.tick_params(axis='both',which='major',direction='inout',length=7)
#
left_axis.grid(which='major',axis='both',linewidth='1.1')
#
#left_axis.plot(x,binomial_pdf)
left_axis.vlines(x,[0],binomial_pdf,linewidth=7)
#
plot.get_current_fig_manager().resize(width,height)
plot.gcf().set_size_inches((0.01*width),(0.01*height))
plot.show()
#
########################################################################
#
# EOF
#
########################################################################
