mfmr_em_alg_hom_noquant	<- function( Yb, X, trace, maxit=1e3, tol=1e-3 ){

	N	<- nrow(Yb)
	Q	<- ncol(X)
	B	<- ncol(Yb)
	alpha	<- matrix( 0, Q, B )

	Xalpha	<- X %*% alpha
	XtXiXt	<- solve(t(X)%*%X) %*% t(X)

	for( it in 1:maxit ){

		olda	<- alpha

		mu0		<- Xalpha

		muY		<- SigY	<- matrix( NA, N, B )
		for( b in 1:B ){
			binout		<- bin_fxn( mu0[,b], sd0=1, zeros=which(Yb[,b]==0) )
			muY	[,b]	<- binout$muY
			SigY[,b]	<- binout$SigY
			rm( binout )
		}
 		Sigb	<- colMeans( SigY )
 		rm( SigY )

		alpha		<- XtXiXt %*% muY
		Xalpha	<- X %*% alpha
		al_eps	<- sqrt( mean((alpha-olda)^2)/mean(alpha^2) )
		if( Q == 1 ) alpha		<- matrix( alpha,		Q, B )

		sqrt.err	<- muY - Xalpha
		Shat		<- t(sqrt.err) %*% sqrt.err / N

		diag( Shat )<- diag( Shat ) + Sigb

		if( it > 10 & al_eps < tol ) break

	}
	ll	<- ll_fxn_hom_noquant( Yb, Xalpha, N, B )

	list( pvec=NA, alpha=alpha, beta=NA, Sigma=diag(B), ll=ll, niter=it, conv=it<maxit )
}

ll_fxn_hom_noquant	<- function( Yb, Xalpha, N, B ){

	mu0			<- Xalpha[,1:B]

	lls	<- rowSums( sapply( 1:B, function(b){
		l.binp_fxn( z = - mu0[,b], zeros=which(Yb[,b]==0) )
	}))

	-N*B/2*log(2*pi) + sum( lls )
}
