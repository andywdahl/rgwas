lpmat2pmat	<- function( pmat ){
	pmat		<- exp(t( apply( pmat, 1, function(x) x - max(x) ) ))
	pmat / matrix( rowSums(pmat), nrow(pmat), ncol(pmat), byrow=FALSE )
}

Lam2Sig	<- function( P, Cmat, Shat_qi ){
	B					<- ncol( Cmat )
	Schur_D		<- solve( Shat_qi )
	Schur_D_C	<- Schur_D %*% Cmat 
	Sigma			<- matrix( NA, P, P )
	Sigma[1:B		,1:B]		<- diag(B) + t(Cmat) %*% Schur_D_C
	Sigma[1:B		,-(1:B)]<- -t(Schur_D_C)
	Sigma[-(1:B),1:B]		<- -Schur_D_C
	Sigma[-(1:B),-(1:B)]<- Schur_D
	Sigma
}

binarize	<- function( y, na.rm=T ){
	b	<- rep( NA, length(y) )
	b[ y == min(y, na.rm=na.rm) ]	<- 0
	b[ y == max(y, na.rm=na.rm) ]	<- 1
	b
}
