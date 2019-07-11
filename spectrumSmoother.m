function final=spectrumSmoother(data,type,wind,thresh)
%SPECTRUM_SMOOTHER(DATA,OPTION,WINDOW,THRESHOLD)
%   simple spectral filter algorithm
%
%   DATA is a column vector of the data (Nx1);
%   OPTION selects the type of filter
%       1 = boxcar average, must specify window (default is 3)
%               WINDOW defined as how much on either side to average (ie, 3 means +/-3)
%       2 = spike killer, average method
%               deviation of spike in % counts replaced by average of
%               neighbors (see WINDOW, THRESHOLD)
%       3 = positive spectrum filter, inserts lowest non-zero value
%       4 = spike killer, nearest neighbor method
%               deviation of spike by THRESHOLD (in % counts) by average of WINDOW nearest
%               neigbbor(s). default window is 3;
%

%   WINDOW              
%       number of nearest neighbors to average
%       NOTE:   IF using option 2, a zero here results in spike replaced by
%       spectrum average   
%   THRESHOLD
%       used in options 2; % deviation from average is replaced with
%       average over WINDOW
%       either option will pad boundaries to return vector of original size & length
%
%   NOTES
%       AEC 03/03


if nargin<3
    wind=3;thresh=0;
elseif nargin<4
    thresh=0;
end


pts = length(data);
final=zeros(1,pts);
temp=final;

if type==1   %boxcar
    %check window appropriate
    if 2*wind>=pts
        disp('WARNING: BOXCAR algorithm has too few points for the selected window')
      %  break
    end
    
    %perform average
    for i=1:pts
        if(i<=wind || i>=pts-wind )
            final(i)=data(i);       %ensures that function returns same number of points (pad here)
        else
            final(i)=sum(data(i-wind:i+wind))./(2*wind+1);    
        end
    end
    
elseif type==2  %spike kill, average method
    avg = sum(data)/pts;
    
    for i=1:pts
        if  i<=wind || i>=pts-wind
            final(i)=data(i);       %pad array so size is the same
        elseif data(i)>=avg*(1+thresh/100)
            if wind    %if not zero, replace with window average (but skip spike!)
                final(i)=(sum(data(i-wind:i-1))+sum(data(i+1:i+wind)))/(2*wind);   %use points around spike
            else        %if zero, replace with spectrum average
                final(i)=avg;
            end
        else    
            final(i)=data(i);    
        end
    end
    
elseif    type==3   %negative spike 
    %find minimum nonzero value
    min_val=min(abs(data));  %find true min that is >=0
    if min_val==0;   %if 0, then find next highest
        for i = 1:pts
            if data(i)<=0
                temp(i)=-1; %make negative so we can filter out easily
            else
                temp(i)=data(i);    %keep value
            end
        end
        %end
        min_val=min(abs(temp));     %use next highest value for min
    end

    for i=1:pts
        if data(i)<min_val
            final(i)=min_val;
        else
            final(i)=data(i);
        end
    end
 
elseif   type==4        %spike killer, nearest neighbor method
    
    for i=1:pts
        if  i<=wind || i>=pts-wind
            final(i)=data(i);
        else
            neighbor_avg=(sum(data(i-wind:i-1))+sum(data(i+1:i+wind)))/(2*wind);
            if  data(i)>(1+thresh/100)*neighbor_avg
                final(i)=neighbor_avg;
            else
                final(i)=data(i);
            end
        end
    end
    
end %end of all functions
final=final';


