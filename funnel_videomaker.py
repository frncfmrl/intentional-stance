
# -*- coding: utf-8 -*-

"""
Created on Fri Jul 17 15:01:05 2020
"""


############ Translation from matlab in progress ###########

import funnel_manifold
import numpy as np

nDots=13

sig1=0
sig2=0
dt=.001

manifold=funnel()
sf=manifold.suggestedFrame

r=np.empty(2,nDots)
r[0,:]=sf('xMin')+(sf('xMax')-sf('xMin'))*np.random.rand(nDots)
r[1,:]=sf('yMin')+(sf('yMax')-sf('yMin'))*np.random.rand(nDots)

###################### PLOTTING ###########################################################

funnelLinewidth=5
tMax=250
showManifold=0

figure

axis tight manual
set(gca,'nextplot','replacechildren')

v = VideoWriter('funnel.avi')
open(v)


for tInd = 1:tMax
    set(gca,'visible','off');
    set(gca,'xtick',[]);
    disp(tInd)
    clf
    sgtitle([' t = ' num2str(tInd)])
    hold on

    ys=[yMin:.05:yMax];
    if showManifold==1
        ys=[yMin:.05:yMax];
        plot(xFunnel(ys),ys,'lineWidth',5,'Color','b');
        plot(-xFunnel(ys),ys,'lineWidth',5,'Color','b');
    end
    
    for n=1:nDots
        r(:,n)=keepInFunnel(evolve(r(:,n),dt,sig1,sig2));
        plot(r(1,n),r(2,n),'or','MarkerFaceColor','r','MarkerSize',14);
    end
    xlim([-1,1]); 
    ylim([yMin,yMax]);
    pause(.05);
    set(gca,'visible','off');
    set(gca,'xtick',[]);
    drawnow

    frame = getframe(gcf);
    writeVideo(v,frame);
