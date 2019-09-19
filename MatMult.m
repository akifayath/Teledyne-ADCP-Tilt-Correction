function MM=MatMult(A,B)  % Multiply matrices along third dimension
nrows=size(A,1);
ncols=size(B,2);

indim=size(A,2);
if indim~=size(B,1)
    error('MatMult:InDim','Inner matrix dimensions should agree')
end
n3=max(size(A,3),size(B,3));
if (size(A,3)~=size(B,3)) && size(A,3)~=1 && size(B,3)~=1
    error('MatMult:Wrong3dim','Third dimension size must be equal or \n at least one variable must have third dimension size equal to one')
end
MM=zeros(nrows,ncols,n3);
for cntrow=1:nrows
    for cntcol=1:ncols
        for cntdim=1:indim
               MM(cntrow,cntcol,:)=MM(cntrow,cntcol,:)+A(cntrow,cntdim,:).*B(cntdim,cntcol,:);
        end
    end
end
end