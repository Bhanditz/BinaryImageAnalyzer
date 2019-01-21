%Name: Kehong You 
%Project: Binary Image Analyzer / Connected Components Labeler 
%Please change the line (11) to change input file!
%Please change line (127) to change output file name!
% Line (23)-(127) is where the connected component are labeled, and
% output is generated. It is displayed and saved as a new tif file.
% Line (130) - (279) is where the various features are calculated,
% and displayed on a formatted chart in the command window when ran.


T = imread('pic1bin.tif'); %***INPUT IMAGE, CAN CHANGE FILE NAME HERE******
T = mat2gray(T);
[M,N] = size(T);
f = zeros(1,100); %equivalence tree
B = zeros(M,N); %binary image
l = 2; %initial label value
a=0;
b=0;
c=0;

%Below is the looping code that will calculate the connect components
% First pass, then followed by the 2nd pass
for i = 2:(M-1)
    for j=2: (N-1)
        if (T(i,j)==1)
           a=B(i,j); 
           b=B(i,j-1); 
           c =B(i-1,j);
  
          if((b==0)&&(c==0)) 
              
            
             B(i,j) = l; 
             l=l+1;
             
              
         
          elseif ((b==0)&&(c>1)) 
              B(i,j)=c;
          elseif ((b>1)&&(c==0))
              B(i,j) =b;
          elseif ((b>1)&&(c>1))
              if(b==c)
                  B(i,j)=b;
              else
                  while (f(b)>0)
                      
                       b=f(b); %find root
                     
                       
                  end
                  while  (f(c)>0)
                      c=f(c); %find root
                  end
                  
                  
                  if (b==c)
                      B(i,j)=b;
                  elseif (b<c)
                      B(i,j) = b;
                      f(c) =b;
                  elseif(b>c)
                     B(i,j) = c; 
                     f(b) = c;
                  end
                  
                  
                  
              end
             
          end
       
        end
    end

    %pause;
  %  imshow(B);
end


%second pass on the image
for i=1:M
    for j=1:N
        if (B(i,j)>0)
            k=B(i,j);
            while(f(k)>0)
                k = f(k);
            end
            B(i,j) = k;
        end
    end
end

CI = zeros(M,N,3); %Created to store the colored version of the output
 
%loop to generate color for each unique label number
for z=2:l
    b = zeros(M,N,3);
red = rand;  % random number in the 0 to 1.0 range
green = rand; % random number in the 0 to 1.0 range
blue = rand;  % random number in the 0 to 1.0 range
b(:,:,1) = red;
b(:,:,2)=green;
b(:,:,3)=blue;


for i=1:M
    for j=1:N
     if B(i,j)==z
         CI(i,j,1) = red;
          CI(i,j,2) =green;
            CI(i,j,3) =blue;
         
     end
     
    end
end

 
    
end

 imshow(CI);
     
 
 
imwrite(CI,'pic1.tif', 'tiff' );   %*************SAVES FILE HERE********


comp = 0;

%counting total number of objects in image
for i=1:100
if (f(i)>0)
comp = comp +1;    
end
end

BB = reshape(B,1,4800); %flatten array
uniq = unique(BB); %getting all unique number for each component
comp = size(uniq);
comp1 = comp(2)-1;%number of components
listdata = zeros(10,comp1); %holds data for part 2
uniq = uniq(2:end); %excludes 0
listdata(1,:) = uniq; %sets first row of list as the label numbers


tempx = 0; %holds temporary values for X during loop
tempy = 0; %holds temporary values for Y during loop
    r0 = M/2; %finding center (center is origin)
    c0 = N/2;



%area, perimeter, position calculations in the following loop:
for z=1:comp1
    for i=1:M-1
        for j=1:N-1    
           
            
            if (B(i,j)==listdata(1,z))          %To find matching label values
            listdata(2,z) =  listdata(2,z)+1;   % If it matches, add to area counter
            
            tempx = tempx + T(i,j)*i;        %sum of T(i,j)*i
            tempy = tempy + T(i,j)*j;        %sum of T(i,j)*j
            
            end
            
            if (B(i,j)==listdata(1,z))    %perimeter
            a = B(i,j); 
            b = B(i,j-1); 
            c = B(i-1,j);
            d = B(i+1,j+1);
            e = B(i+1,j);
            f = B(i,j+1);
            
            g = B(i-1,j-1);
            h = B(i-1,j+1);
            l = B(i+1,j-1);
            
            zz = b+c+d+e+f+g+h+l;
            if((a*8)~=zz)   %when at least one of the 8 adjacent pts is 0
            listdata(5,z) =  listdata(5,z)+1;   
            end
            
            
            end
        end
 
    end
   

   
      %converting i,j yo x,y
   tempx = tempx/listdata(2,z);         
   tempy = tempy/listdata(2,z);        
   
     listdata(3,z) = tempy-c0;        %Position X
      listdata(4,z) = r0-tempx;      %Position Y
   
   tempx = 0; tempy=0;
end

%Theta, 2nd Moments, and Elongation are calculated in the following loop:
tA = 0; %representing A,B,C in eq 2.17~19
tB = 0;
tC = 0;
for z=1:comp1
    for i=1:M-1
        for j=1:N-1    
            if (B(i,j)==listdata(1,z)) 
                
                xcord = j-c0;
                ycord = r0-i;
                
                
                xprime = xcord-listdata(3,z); %eq. 2.15
                yprime = ycord-listdata(4,z); %eq. 2.15
                
                tA = xprime*xprime*T(i,j) + tA; %eq 2.17
                tB = xprime*yprime*T(i,j) +2* tB;  %eq 2.18
                tC = yprime*yprime*T(i,j) + tC; %eq 2.19
            end
            
        end
    end
    
 
    
    deno = tB*tB+(tA-tC)*(tA-tC); 
    deno = sqrt(deno);                  %denominator of eq 2.22
    listdata(6,z) = (asind(tB/deno))/2; %eq 2.22
    theta1 = listdata(6,z);
    theta2 = (acosd((tA-tC)/deno))/2;  %eq. 2.23
    listdata(7,z) = listdata(3,z)*cosd(listdata(6,z))+listdata(4,z)*sind(listdata(6,z)); %rho eq.2.14 
    
    
    Xsquare1 = (0.5*(tA+tC)+.5*(tA-tC)*cosd(2*theta1)+0.5*(tB*sind(2*theta1))); %eq 2.20
    Xsquare2 = (0.5*(tA+tC)+.5*(tA-tC)*cosd(-2*theta2)+0.5*(tB*sind(-2*theta2))); %eq 2.20

    listdata(8,z) = Xsquare1;
     listdata(9,z) = Xsquare2;
     
     if (listdata(8,z) < listdata(9,z) ) %makes sure max 2nd moment is bigger 
         temp = listdata(8,z); %if not, then swap the two values
         listdata(8,z) = listdata(9,z);
         listdata(9,z) = temp;
  
     end
     
      listdata(10,z) = sqrt( listdata(9,z))/sqrt( listdata(8,z)); %elongation eq. 2.24
     
    tA = 0; %reset parameter/temp values for each unique label
    tB = 0;
    tC = 0;
    xprime = 0;
    yprime = 0;

end



%formatting table:
if (comp1 == 9)
 rowNames = {'Unique Label','Area','Position: x','Position: y','Perimeter','Theta °','Rho °','Max 2nd Moment','Min 2nd Moment','Elongation'};
colNames = {'s1','s2','s3','s4','s5','s6','s7','s8','s9'};
sTable = array2table(listdata,'RowNames',rowNames,'VariableNames',colNames);
     disp(sTable); %displays formatted table for part 2  
elseif (comp1 == 7)
rowNames = {'Unique Label','Area','Position: x','Position: y','Perimeter','Theta °','Rho °','Max 2nd Moment','Min 2nd Moment','Elongation'};
colNames = {'s1','s2','s3','s4','s5','s6','s7'};
sTable = array2table(listdata,'RowNames',rowNames,'VariableNames',colNames);
     disp(sTable); %displays formatted table for part 2  

else
     disp(listdata); %displays unformatted table if image does not have 7 or 9 objects for part 2    
end

