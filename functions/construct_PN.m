function [Pre,Post] = construct_PN(adj)
% CONSTRUCT_PN   Construye las matrices Pre y Post de la red de Petri
%   [Pre,Post] = construct_PN(adj)
%   adj : matriz de adyacencia (NxN), con 1 en (i,j) si hay transición i→j
%   Pre : matriz (NxE) con Pre(i,e)=1 si la transición e sale del lugar i
%   Post: matriz (NxE) con Post(j,e)=1 si la transición e entra al lugar j

    n = size(adj,1);

    % 1) Extraer pares (i,j) de todas las aristas, sin bucle i→i
    [I,J] = find(adj);
    keep = I~=J;
    I = I(keep);
    J = J(keep);
    
    % 2) Número de transiciones = número de aristas
    E = numel(I);

    % 3) Montar Pre y Post como sparse
    Pre  = sparse(I, 1:E, 1, n, E);
    Post = sparse(J, 1:E, 1, n, E);
    
    % Si realmente necesitas matrices densas, descomenta:
    %Pre  = full(Pre);
    %Post = full(Post);
end