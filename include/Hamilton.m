%%
% Let us create the following graph
%       (1)--(2)--(3)-------(4)
%        |   / \   |		 |
%        |  /   \  |		 |
%        | /     \ |		 |
%       (5)-------(6)		 |
%		 |		 			 |
%		 |					 |
%		 |					 |
%		(7)-----------------(8)
%   

 % (1,2) = 1 -> certify node 1 and 2 is connected.
% Graph_existance_check=[0 1 0 0 1 0 0 0;
%                        1 0 1 0 1 1 0 0;
%                        0 1 0 1 0 1 0 0;
%                        0 0 1 0 0 0 0 1;
%                        1 1 0 0 0 1 1 0;
%                         0 1 1 0 1 0 0 0;
%                         0 0 0 0 1 0 0 1;
%                         0 0 0 1 0 0 1 0];
% Start = 5; % Source
% Goal = 1; % Destination
%
% P =test_hamilton(Graph_existance_check,Start,Goal);
%
% P will be an array mentioning the path/cycle, if path/cycle found; or a
% string: 'No Path Found', if path/cycle not found
%
% #Note: This code can be used for finding Hamiltonian cycle also. For
% that, make sure Source and Destination are same.

%%
%{
    Main Function
%}
function hamPath = Hamilton(Graph, Source, Destination)

% Input Checking
if ~isreal(Graph)
    error('Graph must be in real form');
elseif ~isnumeric(Graph)
    error('Matrix must be numeric');
elseif ~ismatrix(Graph)
    error('Check Matrix Dimensions');
else
    [r, c] = size (Graph);
    if r~=c
        error('Matrix must be square matrix');
    end
end

if ~(isreal(Source)||isreal(Destination)||(Source>0 && Source<=r) || (Destination>0 && Destination<=r))
    error('improper Source/Destination');
end

deg = 0;
for i = 1:length(Graph(1,:))
    
    if Graph(Source,i) == 1
        deg = deg + 1;
    end
    if Graph(Destination,i) == 1
        deg = deg + 1;
    end 
    
end
N = length(Graph(1,:));
if deg < N
    fprintf("probablly not Hamilton Graph\n");
else
    fprintf("This is Hamilton Graph\n");
end
clear c;

% Function call
hamPath = find_Hamilton(Graph, Source, Destination, r);

end

%%
%{
    This functions sets some initial parameters, and calls the actual
    function.
%}
function hamPath = find_Hamilton(Graph, Source, Destination, totalNodes)

hamPath = zeros(size(Graph(1,:)));

hamPath(1) = Source;

[Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, 1);

if Status == 0
    if Source ~= Destination
        hamPath = 'No Path Found';
    else
        hamPath = 'No Cycle Found';
    end
    return;
end

end

%%
%{
    This function recursively call itself, hence finding the solution
%}
function [Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, nodesFound)

% Ending Condition check
if ( (nodesFound == totalNodes-1 && Source~=Destination) || (nodesFound == totalNodes && Source==Destination) )
    if ( Graph(hamPath(nodesFound), Destination) ~= 0)
        hamPath(nodesFound+1) = Destination;
        Status = 1;
        return;
    else
        Status = 0;
        return;
    end
end

for i=1:totalNodes
    if i==Destination
        continue;
    end
    
    if isSafe(Graph, hamPath, nodesFound, i)
        hamPath(nodesFound+1) = i;
        
        [Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, nodesFound+1);
        if Status
            return;
        end
        
        hamPath(nodesFound+1) = 0;
    end
end

Status = 0;

end

%%
%{
    This function is used to check whether the current node can be added
    or not for making the path/cycle.
%}
function Flag = isSafe(Graph, hamPath, nodesFound, i)

if Graph(hamPath(nodesFound),i) == 0
    Flag = 0;
    return;
end

for ii=1:nodesFound
    if hamPath(ii) == i
        Flag = 0;
        return;
    end
end

Flag = 1;

end