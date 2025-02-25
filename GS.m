%% Gram-Schmidt orthogonalization
function [V] = GS(A)
V(:,1)=A(:,1)/norm(A(:,1));
Num=0;
[Ahang,Alie]=size(A); 
for j=2:Alie
    sum(:,1)=zeros(1,Ahang); 
    for i=1:j-1
        h(i,j)=A(:,j)'*V(:,i);
         sum(:,1) =sum(:,1)+h(i,j)*V(:,i);
            Num=Num+1;
    end
    V(:,j)=A(:,j)-sum(:,1);
    V(:,j)=V(:,j)/norm(V(:,j));
end