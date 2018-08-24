interxn_test	<- function( X, y, g, pmat, bin=FALSE ){

	family	<- ifelse( bin, 'binomial', 'gaussian' )

	gxe			<- apply( pmat, 2, function(x) x*g )

	nullfit	<- glm( y ~ X						, family=family )
	gfit		<- glm( y ~ X + g				, family=family )
	gxefit	<- glm( y ~ X + g + gxe	, family=family )
	pvals	<- c(
		anova( nullfit, gfit	, test='Chisq')$Pr[2],
		anova( gfit		, gxefit, test='Chisq')$Pr[2],
		anova( nullfit, gxefit, test='Chisq')$Pr[2]
	)
	names(pvals)	<- c( 'Hom', 'Het', 'Global' )

	gxe	<- cbind( gxe, (1-rowSums(pmat))*g )
	fit	<- glm( y ~ X + gxe, family=family )

	return(list( pvals=pvals, fit=fit, pvec=colMeans( pmat, na.rm=T ) ))

}
