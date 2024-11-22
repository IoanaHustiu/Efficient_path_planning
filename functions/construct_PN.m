function [Pre,Post,m0] = construct_PN(T)
%construct Pre, Post and initial marking based on partition and initial cells

Post= zeros(size(T.adj,1),[]);
Pre = zeros(size(T.adj,1),[]);

for i = 1 : size(T.adj,1)
    temp = setdiff(find(T.adj(i,:)),i);
    for j = 1 : length(temp)
        % add a transition from place i to temp(j)
%         fprintf(1,'Add a transition from place p%d to p%d\n',i, temp(j));
        Pre = [Pre zeros(size(Pre,1),1)];
        Post = [Post zeros(size(Post,1),1)];
        Pre(i,size(Pre,2)) = 1;
        Post(temp(j),size(Post,2))=1;
    end
end

m0=zeros(length(T.Q),1);    %initial marking (for PN)
for i=T.Q
    m0(i)=sum(ismember(T.R0,i));  %number of robots initially in state i
end

end
