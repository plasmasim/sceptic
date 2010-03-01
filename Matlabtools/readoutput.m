%function readoutput(filename)

% Open the file
    filename
	fid=fopen(filename,'r');
     
% Check if the run has been performed with collisions or not
    if(or(strcmp(filename(end-7),'c')==1,strcmp(filename(end-7),'C')==1))
        ncoll=0;
    else
        ncoll=1;
    end
    % Read the header  

    if(ncoll)
	TS=textscan(fid,'%f',12,'headerlines',1);
    TS=TS{1};
    
    dt=TS(1);vd=TS(2);Ti=TS(3);steps=TS(4);rhoinf=TS(5);
    phiinf=TS(6);fave=TS(7);dbl=TS(8);Vp=TS(9);damplen=TS(10);
    Bz=TS(11);nrused=TS(12);%coll=TS(13);nrused=TS(14);

    else
    TS=textscan(fid,'%f',14,'headerlines',1);
    TS=TS{1};
    
    dt=TS(1);vd=TS(2);Ti=TS(3);steps=TS(4);rhoinf=TS(5);
    phiinf=TS(6);fave=TS(7);dbl=TS(8);Vp=TS(9);damplen=TS(10);
    Bz=TS(11);nrused=TS(12);coll=TS(13);nrused=TS(14);
    
    end

% Read mesh diag informations
    TS=textscan(fid,'%f%f%f',nrused,'headerlines',1);
    rcc=TS{1};diagphi=TS{2};diagrho=TS{3};
 
% Read flux to probe at each step
    TS=textscan(fid,'%d',steps,'headerlines',3);
    fluxprobe=TS{1};
    
% Read particle flux data
    TS=textscan(fid,'%f',2,'headerlines',2);
    TS=TS{1};
    nthused=TS(1);
    format='';
%    for k=1:nthused
%        format=strcat(format,'%d');
%    end
%    TS=textscan(fid,format,steps,'headerlines',1);
%    ninthstep=cell2mat(TS)
    TS=textscan(fid,'%d','headerlines',1);
    
% Read averaged particle flux data
    TS =textscan(fid,'%d',1,'headerlines',1);
    nastep=TS{1};
    TS=textscan(fid,'%d');
    ninth=TS{1};
% fix theta boundaries
    ninth(1)=2*ninth(1);
    ninth(nthused)=2*ninth(nthused);
    
    
        
% Read Mesh potential
    TS=textscan(fid,'%s',1,'delimiter','\n');
    TS{1};
    for k=1:nrused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1);
            phi=TS;
        else
            phi=horzcat(phi,TS);
        end
    end
    phi=phi';
    
% Read Mesh density (normalised to rhoinf)
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:nrused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1)
            rho=TS;
        else
            rho=horzcat(rho,TS);
        end
    end
    rho=rho';

% Read inverse volume of cells
    TS=textscan(fid,'%f','headerlines',2);
    volinv=TS{1};
    
 
	TS=textscan(fid,'%f',8,'headerlines',1);
    TS=TS{1};
    
    dt=TS(1);vd=TS(2);Ti=TS(3);steps=TS(4);rmax=TS(5);rhoinf=TS(6);
    dbl=TS(7);Vp=TS(8);
   
    
% Read psum
    TS=textscan(fid,'%s',2,'delimiter','\n','headerlines',1);
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            psum=TS;
        else
            psum=horzcat(psum,TS);
        end
    end
    psum=psum';
    
% Read vrsum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vrsum=TS;
        else
            vrsum=horzcat(vrsum,TS);
        end
    end
    vrsum=vrsum';
    
% Read vtsum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vtsum=TS;
        else
            vtsum=horzcat(vtsum,TS);
        end
    end
    vtsum=vtsum';
    
% Read vpsum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vpsum=TS;
        else
            vpsum=horzcat(vpsum,TS);
        end
    end
    vpsum=vpsum';
 
    
% Read v2sum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            v2sum=TS;
        else
            v2sum=horzcat(v2sum,TS);
        end
    end
    v2sum=v2sum';
    
% Read vr2sum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vr2sum=TS;
        else
            vr2sum=horzcat(vr2sum,TS);
        end
    end
    vr2sum=vr2sum';
    
% Read vtp2sum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vtp2sum=TS;
        else
            vtp2sum=horzcat(vtp2sum,TS);
        end
    end
    vtp2sum=vtp2sum';
    
    % Read vzsum
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            vzsum=TS;
        else
            vzsum=horzcat(vzsum,TS);
        end
    end
    vzsum=vzsum';
    
% Read diagvr
    TS=textscan(fid,'%s',2,'delimiter','\n');
    
    for k=1:nthused
        TS=textscan(fid,'%f',nrused);
        TS=TS{1};
        if (k==1)
            diagvr=TS;
        else
            diagvr=horzcat(diagvr,TS);
        end
    end
    diagvr=diagvr';
    
% Read crap
    TS=textscan(fid,'%f',nrused,'headerlines',2);
    TS=textscan(fid,'%f','headerlines',2);
    
% Read tcc
    TS=textscan(fid,'%[^t]',1);
    TS=textscan(fid,'%f',nthused,'headerlines',1);
    tcc=TS{1};

% Read various forces ...
   
    
    TS=textscan(fid,'%f','headerlines',2);
    TS=TS{1};
    
    charge1=TS(1);charge2=TS(7);ffield(1)=TS(2);ffield(2)=TS(8);
    felec1=TS(3);felec2=TS(9);fion1=TS(4);fion2=TS(10);
    fcol1=TS(5);fcol2=TS(11);ftot1=TS(6);ftot2=TS(12);
    
    TS=textscan(fid,'%f','headerlines',1+ncoll);
    TS=TS{1};
    fion3=TS(1);
    TS=textscan(fid,'%f','headerlines',1);
    TS=TS{1};
    encoll=TS(1);

% close file

    fclose(fid);
        
%end