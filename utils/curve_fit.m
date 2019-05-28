function [coeff,parameter]=curve_fit(knots,num_ctrl,points,outer_normal)

parameter_original = arclength_parameter(points);
num_points = size(points,2);
 while num_points <2*num_ctrl
     % double the points
     pp = zeros(2,num_points*2 -1);
     p_inter = (points(:,1:end-1)+points(:,2:end))/2;
     pp(:,1:2:end)=points;
     pp(:,2:2:end) = p_inter;
     points = pp;
     
     num_points = size(points,2);
 end

 parameter = arclength_parameter(points);
pt = interparc(parameter,points(1,:),points(2,:),'linear');
points = pt';
dim = size(points,1);
num_points = size(points,2);
coeff = zeros(size(points,1),num_ctrl);
coeff(:,1)=points(:,1);
coeff(:,end)=points(:,end);
basis = zeros(num_points,num_ctrl);
for i = 1:num_ctrl
    c = zeros(1,num_ctrl);
    c(i)=1;
    sp = spmak(knots,c);
    basis(:,i)=fnval(sp,parameter);
end


if nargin == 3
    N = basis(:,2:end-1);
    P = points;
    b = P*N-points(:,[1,end])*basis([1,end],2:end-1);
    b = (P-coeff(:,[1,end])*(basis(:,[1,end]))')*N;
    M = N'*N;
    coeff(:,2:end-1)=(M\b')';
    return;
else
    n1 = outer_normal(:,1);
    n2 = outer_normal(:,2);
    x0 = coeff';
    x0 = x0(:);
    p = points';
    p = p(:);
    H = sparse(numel(x0),numel(x0));
    B = basis'*basis;
    ii=1;
    while ii*num_ctrl<=numel(x0)
        H(((ii-1)*num_ctrl+1):ii*num_ctrl,...
            ((ii-1)*num_ctrl+1):ii*num_ctrl)=B;
        ii = ii+1;
    end
    f = (- points*basis)';
    f = f(:);
    begining_indices = ((1:dim)-1)*num_ctrl+1;
    ending_indices = (1:dim)*num_ctrl;
    eq_indices=[begining_indices,ending_indices];
    beq = x0(eq_indices);
    Aeq= sparse(numel(beq),numel(x0));
    for ii = 1:numel(beq)
        Aeq(ii,eq_indices(ii))=1;
    end
    next_begining_indices = ((1:dim)-1)*num_ctrl+2;
    before_ending_indices = ending_indices-1;
    A_ieq = sparse(2,numel(x0));
    A_ieq(1,next_begining_indices)=n1;
    A_ieq(1,begining_indices)= -n1;
    A_ieq(2,before_ending_indices)=n2;
    A_ieq(2,ending_indices)= -n2;
    b_ieq = zeros(2,1);
    op = optimset('Display','off'); 
%     x = quadprog(H,f,[],[],Aeq,beq,[],[],[],op);
        x = quadprog(H,f,A_ieq,b_ieq,Aeq,beq,[],[],[],op);

    coeff = reshape(x,num_ctrl,dim)';
end
parameter = parameter_original;
end

function parameter = arclength_parameter(points)
parameter = zeros(1,size(points,2));
diff = (points(:,1:end-1)-points(:,2:end));
length = sqrt(sum(diff.^2,1));
parameter(2:end)=cumsum(length);
parameter = parameter/parameter(end);
parameter(end)=1.0;
end
