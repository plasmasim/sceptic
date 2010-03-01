function postproc(filename)

%Postproc(filename). Plots angular distribution of current for file
%filename='... .dat'

    tcc=0;
    rcc=0;
    readoutput;
    nastep
    flux0=sqrt(2*Ti)/(2*sqrt(pi));

    fluxofangle=double(ninth)*double(nthused-1)/(4*pi*rhoinf*dt*double(nastep))/flux0;
    
    figure
    plot(tcc,fluxofangle)
    xlabel('cos\theta','FontSize',18);
    ylabel('\Gamma_i /(v_{ti} / 2\pi^{1/2})','FontSize',18)
    title('Angular current distribution','FontSize',18)
    


end