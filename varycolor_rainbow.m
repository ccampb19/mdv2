function ColorSet=varycolor_rainbow(NumberOfPlots)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.  
%
%     Yellow and White colors were not used because of their poor
%     translation to presentations.
% 
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
% 
%         figure
%         hold on;
% 
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

%Created by Daniel Helmick 8/12/2008

narginchk(1,1)%correct number of input arguements??
nargoutchk(0, 1)%correct number of output arguements??

%Take care of the anomolies
if NumberOfPlots<1
    ColorSet=[];
elseif NumberOfPlots==1
    ColorSet=[0 0 0];
elseif NumberOfPlots==2
    ColorSet=[0 1 0; 0 1 1];
elseif NumberOfPlots==3
    ColorSet=[0 1 0; 0 1 1; 0 0 1];
elseif NumberOfPlots==4
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1];
elseif NumberOfPlots==5
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0];
elseif NumberOfPlots==6
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0; 0 0 0];

else %default and where this function has an actual advantage

    vigbyor = [
        148 0 211;
        75 0 130;
        0 0 255;
        0 255 0;
        255 255 0;
        255 127 0;
        255 0 0;
        ];
    
    ColorSet = interp1(1:7, vigbyor, linspace(1,7,NumberOfPlots))./255;
%     %we have 5 segments to distribute the plots
%     EachSec=floor(NumberOfPlots/5); 
%     
%     %how many extra lines are there? 
%     ExtraPlots=mod(NumberOfPlots,5); 
%     
%     %initialize our vector
%     ColorSet=zeros(NumberOfPlots,3);
%     
%     %This is to deal with the extra plots that don't fit nicely into the
%     %segments
%     Adjust=zeros(1,5);
%     for m=1:ExtraPlots
%         Adjust(m)=1;
%     end
%     
%     SecOne   =EachSec+Adjust(1);
%     SecTwo   =EachSec+Adjust(2);
%     SecThree =EachSec+Adjust(3);
%     SecFour  =EachSec+Adjust(4);
%     SecFive  =EachSec;
% 
%     for m=1:SecOne
%         ColorSet(m,:)=[0 1 (m-1)/(SecOne-1)];
%     end
% 
%     for m=1:SecTwo
%         ColorSet(m+SecOne,:)=[0 (SecTwo-m)/(SecTwo) 1];
%     end
%     
%     for m=1:SecThree
%         ColorSet(m+SecOne+SecTwo,:)=[(m)/(SecThree) 0 1];
%     end
%     
%     for m=1:SecFour
%         ColorSet(m+SecOne+SecTwo+SecThree,:)=[1 0 (SecFour-m)/(SecFour)];
%     end
% 
%     for m=1:SecFive
%         ColorSet(m+SecOne+SecTwo+SecThree+SecFour,:)=[(SecFive-m)/(SecFive) 0 0];
%     end
    
end