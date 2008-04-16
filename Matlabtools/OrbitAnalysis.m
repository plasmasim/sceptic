function OrbitAnalysis(filename)

% Open the file
	fid=fopen(filename,'r');
    
    TS=textscan(fid,'%d',1,'headerlines',1);
    norb=TS{1};
    TS=textscan(fid,'%d',1,'headerlines',1);   
    for (k=1:norb)
       
       TS=textscan(fid,'%d',1,'headerlines',1);
       steps=TS{1};
       Delta=0;
       l=1;
       TS=textscan(fid,'%f%f%f%f%f%f',1);
       Xorbit{k}(l)=TS{1};   
       Yorbit{k}(l)=TS{2};
       Zorbit{k}(l)=TS{3};
       VXorbit{k}(l)=TS{4};
       VYorbit{k}(l)=TS{5};
       VZorbit{k}(l)=TS{5}; 
       Rorbit{k}(l)=sqrt(TS{1}*TS{1}+TS{2}*TS{2});
       eps=sign(TS{1});
       
       RR=5;
       
% Hope this allows to eliminate the risk of plotting two orbits
        while(l<steps)
            
            l=l+1;
            TS=textscan(fid,'%f%f%f%f%f%f',1);
            Delta=abs(TS{3}-Zorbit{k}(l-1))+abs(TS{2}-Yorbit{k}(l-1))+abs(TS{1}-Xorbit{k}(l-1));
            
            if (Delta>=0.4)
                    if (RR<1.5)
                    ctheta(k)=Zorbit{k}(l-1)/RR;
                    h(k)= k*k;
                end
                break
            end
            Xorbit{k}(l)=TS{1};   
            Yorbit{k}(l)=TS{2};
            Zorbit{k}(l)=TS{3};
            VXorbit{k}(l)=TS{4};
            VYorbit{k}(l)=TS{5};
            VZorbit{k}(l)=TS{5}; 
           % Rorbit{k}(l)=eps*sign(TS{1})*sqrt(TS{1}*TS{1}+TS{2}*TS{2});
            Rorbit{k}(l)=sqrt(TS{1}*TS{1}+TS{2}*TS{2});
            RR=sqrt(TS{1}*TS{1}+TS{2}*TS{2}+TS{3}*TS{3});
        end

        TS=textscan(fid,'%f%f%f%f%f%f');

    end
    
    for k=1:100
        Xprobe(k)=cos(pi*k/100);
        Yprobe(k)=sin(pi*k/100);
        Xbound(k)=8*cos(pi*k/100);
        Ybound(k)=8*sin(pi*k/100);
    end 
% close file
   fclose(fid);
    figure;
    for k=1:1:norb
        plot(Zorbit{k},Rorbit{k});
        hold all
    end
    plot(Xprobe,Yprobe,'k','Linewidth',1);
    plot(Xbound,Ybound,'k','Linewidth',1);
    
    axis equal;axis tight;
%    figure;
%    plot(h,ctheta);

%end
     
        
end