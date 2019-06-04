mfmr_em_alg	<- function(
	Yb, Y, G, X, K, trace,
	N, Q, S, P, B,
	pvec, alpha, beta, Sigma,
	maxit=1e3, tol=rep(1e-3,2)
){

	if( is.null( Yb ) )
		return( mfmr_em_alg_nobin( Y, G, X, K, trace, N, Q, S, P, pvec, alpha, beta, Sigma, maxit, tol ) )

	al_eps	<- 0
	Xalpha	<- matrix( 0, N, P )
	if( Q > 0 ){
		Xalpha	<- X %*% alpha
		XtXiXt	<- solve(t(X)%*%X) %*% t(X)
	}

	Lambda_qq	<- solve( Sigma[-(1:B),-(1:B)] )
	Cmat			<- -Lambda_qq %*% Sigma[-(1:B),1:B]
	rm( Sigma )

	if( trace )
		cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", 0, NA, NA, NA, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
	for( it in 1:maxit ){
		olda		<- alpha
		oldb		<- beta

		pmat	<- matrix( NA, N, K )
		muY		<- array( NA, dim=c(N,K,B) )
		for( k in 1:K ){

 			qerr			<- Y - Xalpha[,-(1:B)] - G %*% beta[k,,-(1:B)]
			Lam_qerr	<- qerr %*% Lambda_qq
			mu0				<- Xalpha[,1:B]	+ G %*% beta[k,,1:B] - qerr %*% Cmat

			l.binp	<- matrix( NA, N, B )
			for( b in 1:B ){
				binout			<- bin_fxn( mu0[,b], sd0=1, zeros=which(Yb[,b]==0) )
				l.binp[,b]	<- binout$l.binp
				muY	[,k,b]	<- binout$muY
				rm( binout )
			}
			pmat[,k]	<- log(pvec[k]) + rowSums(l.binp) - 1/2*rowSums( Lam_qerr * qerr )
		}
 		pmat	<- lpmat2pmat( pmat )
		pvec	<- colMeans( pmat )

		if( !is.null(X) ){
			meanmat	<- apply( muY, 3, function(muY_b) rowSums( muY_b * pmat ) )
			delta		<- cbind( meanmat, Y ) - sapply( 1:P, function(p) rowSums( G * (pmat %*% beta[,,p]) ))
			alpha		<- XtXiXt %*% delta
			Xalpha	<- X %*% alpha
			al_eps	<- sqrt( mean((alpha-olda)^2)/mean(alpha^2) )
			if( Q == 1 ) alpha		<- matrix( alpha,		Q, P )
		}

		for( k in 1:K ){
			xG				<- apply( G, 2, function(g) g*pmat[,k] )
			beta[k,,]	<- solve( t(G)%*%xG ) %*% t(xG)%*%as.matrix( cbind(muY[,k,],Y) - Xalpha )
		}
		be_eps	<- sqrt( mean((beta-oldb)^2)/ mean( sapply( 1:K, function(k) c( beta[k,,] - apply( beta, 2:3, mean ) )^2 ) ) )

		errmat	<- 0
		for( k in 1:K ){
			sqrt.err.k	<- apply( cbind(muY[,k,],Y) - Xalpha - G%*%beta[k,,], 2, function(e) e*sqrt(pmat[,k]) )
			errmat			<- errmat + t(sqrt.err.k) %*% sqrt.err.k
		}
		Shat	<- errmat/N
		#diag( Shat )[1:B]<- diag( Shat )[1:B] + Sigb
		Lambda_qq	<- solve( Shat[-(1:B),-(1:B)] )
		Cmat			<- -Lambda_qq %*% Shat[-(1:B),1:B]
		rm( Shat )

		if( trace )
			cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", it, al_eps, be_eps, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
		if( it > 10 & all( c( al_eps, be_eps ) < tol, na.rm=T ) ) break

	}
	Sigma	<- Lam2Sig( P, Cmat, Lambda_qq )
	ll	<- ll_fxn( Yb, Y, Xalpha, G, beta, Sigma, pvec, N, P, K, B )
	list( pvec=pvec, alpha=alpha, beta=beta, Sigma=Sigma, ll=ll, niter=it, conv=it<maxit )
}

ll_fxn	<- function( Yb, Y, Xalpha, G, beta, Sigma, pvec, N, P, K, B ){

	Lambda_qq		<- solve( Sigma[-(1:B),-(1:B)] )

	all_lls_NxK	<- sapply( 1:K, function(k){

		qerr		<- Y - Xalpha[,-(1:B)] - G %*% beta[k,,-(1:B)]
		Lam_qerr<- qerr %*% Lambda_qq
		mu0			<- Xalpha[,1:B]	+ G %*% beta[k,,1:B] + Lam_qerr %*% Sigma[-(1:B),1:B]

		log(pvec[k]) -1/2*rowSums( Lam_qerr * qerr ) + rowSums( sapply( 1:B, function(b){
			l.binp_fxn( z = - mu0[,b], zeros=which(Yb[,b]==0) )
		}))

	})
	offsets		<- apply( all_lls_NxK, 1, max )
	all_lls_N	<- sum(log(rowSums(exp( all_lls_NxK - (offsets %o% rep(1,K)) )))) + sum(offsets)

	-N*P/2*log(2*pi) + N/2*as.numeric(determinant(Lambda_qq,logarithm=TRUE)$mod) + all_lls_N
}
