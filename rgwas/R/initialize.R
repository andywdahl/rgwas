initialize_fxn	<- function(
	Yb, Yq, G, X, K,
	init_sd, seed,
	alpha0, beta0, p0, Sigma0,
	...
){
	set.seed( seed )

	B	<- ifelse( is.null( Yb ), 0, ncol(Yb) )
	P	<- ifelse( is.null( Yq ), 0, ncol(Yq) ) + B
	N	<- ifelse( is.null( Yq ), nrow(Yb), nrow(Yq) )
	Q	<- ifelse( is.null( X ), 0, ncol(X) )
	S	<- ncol(G)

	if( !(	missing(alpha0) & missing(beta0) & missing(p0) & missing(Sigma0) ) ){
		if(		missing(alpha0) | missing(beta0) | missing(p0) | missing(Sigma0) )
			stop( 'Init all or nothing' )
		return( list( al=alpha0, be=beta0, Si=Sigma0, p=p0, B=B, P=P, N=N, Q=Q, S=S, K=K ) )
	}

	pars	<- list( al=NA, be=array(0,dim=c(K,S,P)), Si=NA, p=rep(1/K,K) )
	if( Q > 0 ){
		al_b	<- matrix( 0, Q, B )

		if( is.null( Yq ) ){
			al_q	<- NULL
		} else {
			al_q	<- solve( t(X) %*% X ) %*% ( t(X) %*% Yq )
		}
		pars$al	<- cbind( al_b, al_q )
	}
	pars$be[,1,]		<- init_sd*rnorm( K*P )

	if( is.null( Yb ) ){
		pars$Si	<- var( Yq )
	} else if( is.null( Yq ) ){
		pars$Si	<- NA #diag(P)
	} else {
		Lambda	<- diag(P)
		Lambda[-(1:B),-(1:B)]	<- solve( var( Yq ) )
		pars$Si	<- solve( Lambda )
	}
	c( pars, list( B=B, P=P, N=N, Q=Q, S=S, K=K ) )
}
