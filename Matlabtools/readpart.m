%function readpart(filename)

% Open the file
    filename
	fid=fopen(filename,'r');
     

    % Read the header  
    
    dVZ=0.05;
    VZ=-5:dVZ:5;
    Nf=prod(size(VZ));
    fz1=zeros(Nf,1);
    fth1=zeros(Nf,1);
    fz2=zeros(Nf,1);
    fth2=zeros(Nf,1);
    fz3=zeros(Nf,1);
    fth3=zeros(Nf,1);
    fz4=zeros(Nf,1);
    fth4=zeros(Nf,1);
    Ncount1=0;Ncount2=0;Ncount3=0;Ncount4=0;
    
    imax
    for ii=0:imax;%35;
        
  % strcat(filename,'00',num2str(ii),'.dat')
  if(ii<10)
        fid=fopen(strcat(filename,'00',num2str(ii),'.dat'),'r');   
  else
        fid=fopen(strcat(filename,'0',num2str(ii),'.dat'),'r'); 
  end
    TS=textscan(fid,'%f',6);
    TS=TS{1};
    
    npart=TS(2);
    
 
    for k=1:npart
       if(mod(k,10000)==0)

           disp(strcat(num2str((ii*npart+k)*100/(npart*(1+imax))),' %'))
       end
        TS=textscan(fid,'%f',6);
        TS=TS{1};
        x=TS(1);
        y=TS(2);
        z=TS(3);
        vx=TS(4);
        vy=TS(5);
        vz=TS(6);
        
        
        rho=sqrt(x^2+y^2);
        r=sqrt(rho^2+z^2);
        phi=asin(y/(rho+0.0001));
        theta=acos(z/(rho+0.0001));
        vth=vy*sin(phi)+vx*cos(phi);
        vr=[vx vy vz]*[x y z]'/r;
     %   if(and(z>20,and(z<28,rho<0.5)));
     
       % if(z^2+rho^2>20^2)
        %if(and(z>0,and(z<1,rho^2<0.5)))
        if(and(rho^2+z^2>4.95^2,vr<0))
            Ncount1=Ncount1+1;
            if(vz<VZ(1))
                fz1(1)=fz1(1)+1;
            elseif(vz>VZ(Nf))
                fz1(Nf)=fz1(Nf)+1;
            else
                fz1(floor((vz-VZ(1))/dVZ+1))=fz1(floor((vz-VZ(1))/dVZ)+1)+1;
            end
            if(vth<VZ(1))
                fth1(1)=fth1(1)+1;
            elseif(vth>VZ(Nf))
                fth1(Nf)=fth1(Nf)+1;
            else
                fth1(floor((vth-VZ(1))/dVZ+1))=fth1(floor((vth-VZ(1))/dVZ)+1)+1;
            end
        elseif(and(rho^2+z^2>4.95^2,vr>0))
            Ncount2=Ncount2+1;
            if(vz<VZ(1))
                fz2(1)=fz2(1)+1;
            elseif(vz>VZ(Nf))
                fz2(Nf)=fz2(Nf)+1;
            else
                fz2(floor((vz-VZ(1))/dVZ+1))=fz2(floor((vz-VZ(1))/dVZ)+1)+1;
            end
            if(vth<VZ(1))
                fth2(1)=fth2(1)+1;
            elseif(vth>VZ(Nf))
                fth2(Nf)=fth2(Nf)+1;
            else
                fth2(floor((vth-VZ(1))/dVZ+1))=fth2(floor((vth-VZ(1))/dVZ)+1)+1;
            end
         elseif(and(z>3,and(z<4,rho^2<1)))
            Ncount3=Ncount3+1;
            if(vz<VZ(1))
                fz3(1)=fz3(1)+1;
            elseif(vz>VZ(Nf))
                fz3(Nf)=fz3(Nf)+1;
            else
                fz3(floor((vz-VZ(1))/dVZ+1))=fz3(floor((vz-VZ(1))/dVZ)+1)+1;
            end
            if(vth<VZ(1))
                fth3(1)=fth3(1)+1;
            elseif(vth>VZ(Nf))
                fth3(Nf)=fth3(Nf)+1;
            else
                fth3(floor((vth-VZ(1))/dVZ+1))=fth3(floor((vth-VZ(1))/dVZ)+1)+1;
            end
         elseif(and(z>4,and(z<5,rho^2<1)))
            Ncount4=Ncount4+1;
            if(vz<VZ(1))
                fz4(1)=fz4(1)+1;
            elseif(vz>VZ(Nf))
                fz4(Nf)=fz4(Nf)+1;
            else
                fz4(floor((vz-VZ(1))/dVZ+1))=fz4(floor((vz-VZ(1))/dVZ)+1)+1;
            end
            if(vth<VZ(1))
                fth4(1)=fth4(1)+1;
            elseif(vth>VZ(Nf))
                fth4(Nf)=fth4(Nf)+1;
            else
                fth4(floor((vth-VZ(1))/dVZ+1))=fth4(floor((vth-VZ(1))/dVZ)+1)+1;
            end
        end
    end
    
    fclose(fid);
    
    end

 %   end
 
 
         
         
         