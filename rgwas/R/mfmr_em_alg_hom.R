mfmr_em_alg_hom	<- function( Yb, Yq, X, trace, maxit=1e3, tol=1e-3, nrun=1 ){
	if( nrun != 1 ){
		warning( 'Hom MFMR only supported for nrun=1' )
		nrun	<- 1
	}

	if( is.null( Yq ) )
		return( mfmr_em_alg_hom_noquant( Yb, X, trace, maxit=1e3, tol=1e-3 ) )

	N	<- nrow(cbind(Yb,Yq))
	Q	<- ncol(X)
	B	<- ncol(Yb)
	P	<- ncol(cbind(Yb,Yq))
	alpha	<- matrix( 0, Q, P )
	Sigma	<- diag(P)

	Xalpha	<- X %*% alpha
	XtXiXt	<- solve(t(X)%*%X) %*% t(X)

	Lambda_qq		<- solve( Sigma[-(1:B),-(1:B)] )
	Lambda_fast	<- solve( Sigma[1:B,1:B] - Sigma[1:B,-(1:B)] %*% Lambda_qq %*% Sigma[-(1:B),1:B] )

	for( it in 1:maxit ){

		olda	<- alpha

		Lam_qerr	<- ( Yq - Xalpha[,-(1:B)] ) %*% Lambda_qq
		mu0				<- Xalpha[,1:B]	+ Lam_qerr %*% Sigma[-(1:B),1:B]

		muY		<- SigY	<- matrix( NA, N, B )
		for( b in 1:B ){
			binout		<- bin_fxn( mu0[,b], sd0=1/sqrt(Lambda_fast[b,b]), zeros=which(Yb[,b]==0) )
			muY	[,b]	<- binout$muY
			SigY[,b]	<- binout$SigY
			rm( binout )
		}
 		Sigb	<- colMeans( SigY )
 		rm( SigY )

		alpha		<- XtXiXt %*% cbind( muY, Yq )
		Xalpha	<- X %*% alpha
		al_eps	<- sqrt( mean((alpha-olda)^2)/mean(alpha^2) )
		if( Q == 1 ) alpha		<- matrix( alpha,		Q, P )

		sqrt.err	<- cbind(muY,Yq) - Xalpha
		Shat			<- t(sqrt.err) %*% sqrt.err / N

		diag( Shat )[1:B]<- diag( Shat )[1:B] + Sigb
		Shat_qi	<- solve( Shat[-(1:B),-(1:B)] )

		Cmat	<- -Shat_qi %*% Shat[-(1:B),1:B]
		Dmat	<- Shat_qi + Cmat %*% t(Cmat)

		Lam.tmp									<- diag( P )
		Lam.tmp[-(1:B),1:B]			<- Cmat
		Lam.tmp[1:B		,-(1:B)]	<- t(Cmat)
		Lam.tmp[-(1:B),-(1:B)]	<- Dmat
		Sigma	<- solve( Lam.tmp  )
		#Sigma	<- Lam2Sig( P, Cmat, Lambda_qq )
		rm( Lam.tmp )

		Lambda_qq		<- solve( Sigma[-(1:B),-(1:B)] )
		Lambda_fast	<- solve( Sigma[1:B,1:B] - Sigma[1:B,-(1:B)] %*% Lambda_qq %*% Sigma[-(1:B),1:B] )

		if( it > 10 & al_eps < tol ) break

	}
	ll	<- ll_fxn_hom( Yb, Yq, Xalpha, Sigma, N, P, B )

	list( pvec=NA, alpha=alpha, beta=NA, Sigma=Sigma, ll=ll, niter=it, conv=it<maxit )
}

ll_fxn_hom	<- function( Yb, Yq, Xalpha, Sigma, N, P, B ){
	if( is.null(Yq) )
		return( ll_fxn_hom_noquant( Yb, Xalpha, N, B ) )

	Lambda_qq		<- solve( Sigma[-(1:B),-(1:B)] )

	qerr		<- Yq - Xalpha[,-(1:B)]
	Lam_qerr<- qerr %*% Lambda_qq
	mu0			<- Xalpha[,1:B]	+ Lam_qerr %*% Sigma[-(1:B),1:B]

	lls	<- -1/2*rowSums( Lam_qerr * qerr ) +
	rowSums( sapply( 1:B, function(b){
		l.binp_fxn( z = - mu0[,b], zeros=which(Yb[,b]==0) )
	}))

	-N*P/2*log(2*pi) + N/2*as.numeric(determinant(Lambda_qq,logarithm=TRUE)$mod) + sum( lls )
}
