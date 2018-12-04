function RMAT = mrhobb(A,B,MUA,RHO,LAMBDA,NIND,REFF)

musp = A.*LAMBDA.^(-B);
tempr = zeros(length(RHO),length(LAMBDA));
for i = 1:length(RHO)
    tempr(i,:) = Rtheory(MUA,musp,RHO(i),NIND,0,REFF);
end
RMAT = tempr(2:end,:)./tempr(1:end-1,:);
        