function readlive(filename)

% Open the file
    file=filename;
    short=false;new=false;
    filename=strcat(filename,'.dat');
    readoutput();
	fid=fopen(strcat(file,'.liv'),'r');
    
    t(1)=0;
    for i=1:300
        
    
% Read the header     

        TS=textscan(fid,'%s%s',1);
        TS=textscan(fid,'%f',2);
        TS=TS{1};
        if (numel(TS)==0)
            break
        end
        step(i)=TS(1);dt(i)=TS(2);
        t(i+1)=t(i)+dt(i);

 % Read Mesh potential
   TS=textscan(fid,'%s',1);

    for k=1:nrused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1);
            p=TS;
        else
            p=horzcat(p,TS);
        end
    end
    p';

    Phi(i,:,:)=p';
    
    end
% close file
    fclose(fid);
    figure;
    set(gca,'nextplot','replacechildren');
    for k=1:10:i-1
        plot(rcc,Phi(k,:,1));
        axis([1,rcc(nrused),-1.5,0.5])
        text(3/4*rcc(nrused),-1,num2str(t(k)),'FontSize',16)
        %pause;
       F(k)=getframe;
    end
    waitforbuttonpress();
    movie(F,5,1);
      
    end