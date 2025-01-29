function [Pre,Post] = construct_PN(adj)
%construct Pre, Post

Post= zeros(size(adj,1),[]);
Pre = zeros(size(adj,1),[]);

for i = 1 : size(adj,1)
    temp = setdiff(find(adj(i,:)),i);
    for j = 1 : length(temp)
        Pre = [Pre zeros(size(Pre,1),1)];
        Post = [Post zeros(size(Post,1),1)];
        Pre(i,size(Pre,2)) = 1;
        Post(temp(j),size(Post,2))=1;
    end
end

end
