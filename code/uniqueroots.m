function PossRoots = uniqueroots(AllRoots,tol) 
    tol2 = tol; tol3 = 0.02;
    % Remove imag roots
    iRoots = abs(imag(AllRoots));
    if (max(max(iRoots))>tol)
        disp('Imaginary roots! Ignoring ... ');
        disp(sum(iRoots>tol))
        idx = iRoots<tol; % 2D matrix with 1 where real, 0 where imag
        idxrow = min(idx,[],2); % find min in each row, i.e. set 0 if any component is imag
        AllRoots = real(AllRoots(idxrow,:));
    else
        AllRoots = real(AllRoots);
    end
    % Remove roots with zi1 or zi2 outside the domain (-0.5,0.5)
    idxrow1 = abs(AllRoots(:,1))<=0.5;
    idxrow2 = abs(AllRoots(:,2))<=0.5;
    idxrow = min([idxrow1 idxrow2],[],2);
    AllRoots = AllRoots(idxrow,:);

    % Wrap angles to (-pi,pi) and eliminate repeated roots
    AllRoots(:,3) = wrapToPi(AllRoots(:,3));
    AllRoots(:,4) = wrapToPi(AllRoots(:,4));

    AllRoots(abs(AllRoots(:,3)+pi)<tol2,3)=pi; % set -pi equal to pi to avoid repetition
    AllRoots(abs(AllRoots(:,4)+pi)<tol2,4)=pi;
    
    PossRoots = uniquetol(AllRoots, tol3, 'DataScale', 1,'Byrows',true);
end