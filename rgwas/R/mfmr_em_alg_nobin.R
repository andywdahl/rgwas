mfmr_em_alg_nobin	<- function(
	Y, G, X, K, trace,
	N, Q, S, P,
	pvec, alpha, beta, Sigma,
	maxit=1e3, tol=rep(1e-3,2)
){

	al_eps	<- 0
	Xalpha	<- matrix( 0, N, P )
	if( Q > 0 ){
		Xalpha	<- X %*% alpha
		XtXiXt	<- solve(t(X)%*%X) %*% t(X)
	}

	Lambda	<- solve( Sigma )
	rm( Sigma )

	if( trace )
		cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", 0, NA, NA, NA, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
	for( it in 1:maxit ){
		olda		<- alpha
		oldb		<- beta

		pmat	<- matrix( NA, N, K )
		for( k in 1:K ){
 			qerr			<- Y - Xalpha - G %*% beta[k,,]
			Lam_qerr	<- qerr %*% Lambda
			pmat[,k]	<- log(pvec[k])  - 1/2*rowSums( Lam_qerr * qerr )
		}
 		pmat	<- lpmat2pmat( pmat )
		pvec	<- colMeans( pmat )

		if( !is.null(X) ){
			delta		<- Y - sapply( 1:P, function(p) rowSums( G * (pmat %*% beta[,,p]) ))
			alpha		<- XtXiXt %*% delta
			Xalpha	<- X %*% alpha
			al_eps	<- sqrt( mean((alpha-olda)^2)/mean(alpha^2) )
			if( Q == 1 ) alpha		<- matrix( alpha,		Q, P )
		}

		for( k in 1:K ){
			xG				<- apply( G, 2, function(g) g*pmat[,k] )
			beta[k,,]	<- solve( t(G)%*%xG ) %*% t(xG)%*%as.matrix( Y - Xalpha )
		}
		be_eps	<- sqrt( mean((beta-oldb)^2)/ mean( sapply( 1:K, function(k) c( beta[k,,] - apply( beta, 2:3, mean ) )^2 ) ) )

		errmat	<- 0
		for( k in 1:K ){
			sqrt.err.k	<- apply( Y - Xalpha - G%*%beta[k,,], 2, function(e) e*sqrt(pmat[,k]) )
			errmat			<- errmat + t(sqrt.err.k) %*% sqrt.err.k
		}
		Lambda	<- solve( errmat/N )

		if( trace )
			cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", it, al_eps, be_eps, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
		if( it > 10 & all( c( al_eps, be_eps ) < tol, na.rm=T ) ) break

	}
	ll	<- ll_fxn_nobin( Y, Xalpha, G, beta, Lambda, pvec, N, P, K )
	list( pvec=pvec, alpha=alpha, beta=beta, Sigma=solve(Lambda), ll=ll, niter=it, conv=it<maxit )
}

ll_fxn_nobin	<- function( Y, Xalpha, G, beta, Lambda, pvec, N, P, K ){

	all_lls_NxK	<- sapply( 1:K, function(k){
		qerr		<- Y - Xalpha - G %*% beta[k,,]
		Lam_qerr<- qerr %*% Lambda
		log(pvec[k]) -1/2*rowSums( Lam_qerr * qerr )
	})
	offsets		<- apply( all_lls_NxK, 1, max )
	all_lls_N	<- sum(log(rowSums(exp( all_lls_NxK - (offsets %o% rep(1,K)) )))) + sum(offsets)

	-N*P/2*log(2*pi) + N/2*as.numeric(determinant(Lambda,logarithm=TRUE)$mod) + all_lls_N
}
