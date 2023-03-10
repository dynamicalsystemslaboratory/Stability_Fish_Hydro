function maxlines = plot_lines(X,al,type,dis)
% X is the 2D array of all scatter points (:,ialph)
% each column is for a particular alpha
% al is the 1D array of alpha values
% type= 1: Stable, -1:Unstable
% dis: tolerance distance 
% Sample case ---
%     X = [1 2 3 3.1 1.45; 
%          0 1.1 2.1 1.3 2.2; 
%          0 0 1.2 2.15 3.21;
%          0 0 0 0 1.2]
%     al= [1 2 3 4 5]; dis=0.2

    maxn = size(X,1); Nal = size(X,2);
    Xlines = zeros(maxn,Nal);

    ial = 1;
    n = 0; %n=find(X(:,ial),1,'last');

    for ial = 1:Nal
        nprev = n;
        n=find(X(:,ial),1,'last');
        if (isempty(n))
            n=0;
        end

        Xtmp = X(1:n,ial); % Xtmp== zi1S, Xlines==ziS
        if (nprev==0 && n>0) % for first set of points to plot as separate lines
            Xlines(1:n,ial) = sort(Xtmp);
        else
            flag_assignedj = zeros(maxn,1); % Xlineprev(jj) has been assigned to
            for ii=1:n
                flag_assigned = 0; % Xtmp(ii) has been assigned to a line
                for jj=1:maxn
                    if ( (abs(Xtmp(ii)-Xlines(jj,ial-1)) < dis) && flag_assignedj(jj)==0) % if close to not-assigned line
                        Xlines(jj,ial) = Xtmp(ii);
                        flag_assigned = 1; flag_assignedj(jj)=1;
%                         disp('Here1')
                        break;
                    elseif ( (abs(Xtmp(ii)-Xlines(jj,ial-1)) < dis)) % if close to already assigned line
                        jjtmp = find(Xlines(:,ial-1),1,'last')+1;
                        Xlines(jjtmp,ial) = Xtmp(ii); % create new branch
                        Xlines(jjtmp,ial-1) = Xlines(jj,ial-1); % connect new branch to the main line
                        flag_assigned = 1; flag_assignedj(jjtmp)=1;
%                         disp('Here2')
                        break;
                    end
                end
                if (flag_assigned == 0)
                    Xlines(find(Xlines(:,ial-1),1,'last')+1,ial) = Xtmp(ii);
                end
            end
        end
        
    end

    % max number of lines we need to plot
    maxlines = find(sum(abs(Xlines),2),1,'last');

    % Set bottom zeros to nan before plotting
    for ial = 1:Nal
        n=find(Xlines(:,ial),1,'last');
        if (isempty(n))
            n=0;
        end
        if (n<maxn)
            Xlines(n+1:maxn,ial) = nan;
        end
    end
    Xlines(Xlines==0)=nan;

    if (type>0)
        for i=1:maxlines
            plot(al,Xlines(i,:),'-g'); hold on;
        end
    else
        for i=1:maxlines
            plot(al,Xlines(i,:),'-r'); hold on;
        end
    end

end