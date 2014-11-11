rm(list=ls())

atf <- function(reflectance, condition)
    {
    ## assigns luminance to povray reflectance values for specified condition
    ## here practically reflectance and x are always identical
    atf = read.table(paste('luminance_measurements.txt', sep=''), head=T)
    atf.means = tapply(atf$luminance, data.frame(atf$reflectance, atf$medium), mean)
    x   = as.numeric(row.names(atf.means))
    atf.fit = lm(atf.means[,condition] ~ x)
    return(reflectance * atf.fit$coefficients[2] + atf.fit$coefficients[1])
    }


## 1. compute global fit 'pretending' all data came from the same population
make_single_data=function(agg_matches)
    {
    common_data = as.data.frame(matrix(NA, nrow=60, ncol=4, dimnames=list(c(1:60), c('medium', 'x_reflect', 'x_lum', 'match'))))
    x_values=as.numeric(row.names(agg_matches))
    count = 1
    for (my_medium in colnames(agg_matches))
        {
        idx = (count-1)*10
        common_data$medium[   (idx+1):(idx+10)] = rep(my_medium, 10)
        common_data$x_reflect[(idx+1):(idx+10)] = x_values
        common_data$match[    (idx+1):(idx+10)] = agg_matches[, my_medium]
        count = count + 1
        }

    for (k in c(1:60))
        {
        common_data$x_lum[k] = atf(common_data$x_reflect[k], common_data$medium[k])
        }
    return(common_data)
    }

my_orange    = '#FDAE61'
my_red       = '#D7191C'
my_darkblue  = '#2C7BB6'
my_lightblue = '#ABD9E9'

all_cols = list(plain = 'black', shadow='grey', transp_med=my_red, transp_med_steep=my_orange, transp_dark=my_darkblue, transp_dark_steep=my_lightblue)

observers= c('if','js','sk','wd')

global_rsq = matrix(NA, nrow=4, ncol=1, dimnames=list(observers, 'R_sq') )
slopes     = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
intercepts = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
rsq        = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))

svg('figs/ltf.svg', w=18,h=5, pointsize=8)

par(mfrow=c(1,4), lwd=1.5,tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=3, font=1, font.axis=3, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

for (obs in observers)
	
	{
	## combine log and design data
	log_data<-read.table(paste('data/log_', obs, '.txt', sep=''),header=TRUE)
	
	## calculating matches as a function of target reflectance, target location and medium
	m.matches <-tapply(log_data$match_lum, data.frame(log_data$target_lum, log_data$location, log_data$medium), mean)
	sd.matches<-tapply(log_data$match_lum, data.frame(log_data$target_lum, log_data$location, log_data$medium), sd)
	n.matches <-tapply(log_data$match_lum, data.frame(log_data$target_lum, log_data$location, log_data$medium), length)
	
	## aggregate across locations
	agg.m.matches = matrix(0, nrow=10, ncol=6)
	row.names(agg.m.matches) = sort(unique(log_data$target_lum))
	colnames(agg.m.matches)  = sort(unique(log_data$medium))
	
	agg.sd.matches = agg.m.matches
	agg.n.matches  = agg.m.matches
	
	for (curr_med in colnames(agg.m.matches))
		{
		agg.m.matches[ ,curr_med] = apply(m.matches[,,curr_med], 1, mean, na.rm=T)
		agg.sd.matches[,curr_med] = apply(sd.matches[,,curr_med],1, mean, na.rm=T)
		agg.n.matches[ ,curr_med] = apply(n.matches[,,curr_med], 1, mean, na.rm=T)
		}
	
	agg.sem.matches = agg.sd.matches/sqrt(agg.n.matches)
	
    ## rearrange data for global fit
    single_data     = make_single_data(agg.m.matches)
    single.fit      = summary(lm(single_data$match ~ single_data$x_lum))
    global_rsq[obs,]= single.fit$r.squared

	x_values=as.numeric(row.names(agg.m.matches))
	
	## PLAIN data points
	plot(atf(x_values,  'plain'),  agg.m.matches[,'plain'],  pch=21, bg='white', col='black', xlab='local luminance', ylab='', yaxt='n', xlim=c(0,432), ylim=c(0,420))
	axis(2, at=c(0,100,200,300,400),  cex=1.5)

	if (obs =='if')
		{
 		mtext(2, line=2, text='match luminance', cex=2.2)
		}
	
	## loop over linear fits and data points for all conditions
	for (curr_med in colnames(agg.m.matches))
		{
		curr.fit = lm(agg.m.matches[,curr_med] ~ atf(x_values,  curr_med))
        abline(coef(curr.fit), col=unlist(all_cols[curr_med]))
		arrows(atf(x_values, curr_med), agg.m.matches[,curr_med]-agg.sem.matches[,curr_med], atf(x_values,   curr_med), agg.m.matches[,curr_med]+agg.sem.matches[,curr_med], length=0.01, angle=90, code=3)
		points(atf(x_values, curr_med), agg.m.matches[,curr_med], pch=21, bg=unlist(all_cols[curr_med]), col='black')
        
        slopes[obs, curr_med]     = summary(curr.fit)$coefficients[2]
        intercepts[obs, curr_med] = summary(curr.fit)$coefficients[1]
        rsq[obs, curr_med]        = summary(curr.fit)$r.squared
        }

	legend('topleft', obs, bty='n', cex=2)
}

legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=.9, cex=1.2)

dev.off()

write.table(format(global_rsq, digits=4), 'stat_out/ltf_global_rsq.txt',quote=FALSE)
write.table(format(slopes, digits=4),     'stat_out/ltf_slopes.txt',    quote=FALSE)
write.table(format(intercepts, digits=4), 'stat_out/ltf_intercepts.txt',quote=FALSE)
write.table(format(rsq, digits=4),        'stat_out/ltf_rsq.txt',       quote=FALSE)



