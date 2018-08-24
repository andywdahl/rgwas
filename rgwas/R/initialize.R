initialize_fxn	<- function(
	Yb, Y, G, X, K,
	init_sd, seed,
	alpha0, beta0, p0, Sigma0,
	...
){
	set.seed( seed )

	B	<- ifelse( is.null( Yb ), 0, ncol(Yb) )
	P	<- ncol(Y) + B
	N	<- nrow(Y)
	Q	<- ifelse( is.null( X ), 0, ncol(X) )
	S	<- ncol(G)

	if( !(	missing(alpha0) & missing(beta0) & missing(p0) & missing(Sigma0) ) ){
		if(		missing(alpha0) | missing(beta0) | missing(p0) | missing(Sigma0) )
			stop( 'Init all or nothing' )
		return( list( al=alpha0, be=beta0, Si=Sigma0, p=p0, B=B, P=P, N=N, Q=Q, S=S, K=K ) )
	}

	pars	<- list( al=NA, be=array(0,dim=c(K,S,P)), Si=NA, p=rep(1/K,K) )
	if( Q > 0 )
		pars$al	<- cbind( matrix( 0, Q, B ), matrix( solve( t(X) %*% X ) %*% ( t(X) %*% Y ), Q, P-B ) )
	pars$be[,1,]		<- init_sd*rnorm( K*P )
	if( B == 0 ){
		pars$Si	<- var( Y )
	} else {
		Lambda	<- diag(P)
		Lambda[-(1:B),-(1:B)]	<- solve( var(Y) )
		pars$Si	<- solve( Lambda )
	}
	c( pars, list( B=B, P=P, N=N, Q=Q, S=S, K=K ) )
}
