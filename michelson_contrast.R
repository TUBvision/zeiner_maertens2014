rm(list=ls())

## CONTRAST
## compute surround intensities for eight adjacent checkerboards at target positions
compute_surround<-function(data, my_medium, my_location)
    {
	## create list containing target and surround intensities, and match settings 
	## ==========================================================================
    ## select data for medium eg. plain and check location eg. f2
    curr_data = subset(data, medium==my_medium & location==my_location)
    
    ## contrast based on EIGHT adjacent checks
    surr_checks = list(e3=c('e2','d3','e4','f3', 'd2','d4','f2','f4'), f2=c('f1','e2','f3','g2', 'e1','e3','g1','g3'), f3=c('f2','e3','f4','g3', 'e2','e4','g2','g4') )

    ## select surround coordinates for current check location
    curr_surr = unlist(surr_checks[my_location])
    
    ## initialize matrices for surround and target grayvalues [0,1]
    surr_int           = matrix(0, nrow=dim(curr_data)[1], length(curr_surr))
    colnames(surr_int) = curr_surr
    targ_int           = rep(0, dim(curr_data)[1])
    
    ## read out target and surround intensities from corresponding columns
    for (idx in 1:dim(curr_data)[1])
        {
        line = curr_data[idx,]
        targ_int[idx] = as.numeric(line[paste(as.character(line$location), '_v', sep='')])
        surr_count = 1
        for (s in curr_surr)
            {
            surr_int[idx, surr_count]=as.numeric(line[paste(s, '_v', sep='')])
            surr_count = surr_count+1
            }
        }
    return(list(targ=targ_int, surr=surr_int, match=curr_data$match_lum))
    }

compute_contrast<-function(surr_list)
    {
    ## calculate variant of Rayleigh contrast with checks weighted equally
	weighted_surr = apply(surr_list$surr, 1, mean)
    
    ## initialize output variables
    surr_list$contrasts = (surr_list$targ-weighted_surr)/(surr_list$targ+weighted_surr)
    return(surr_list)
    }

bin_eq_contrast<-function(surr_f2, surr_f3, surr_e3, n_per_bin=12)
	{
	##
	matches    = c(surr_f2$match    , surr_f3$match,     surr_e3$match)
	contrasts  = c(surr_f2$contrasts, surr_f3$contrasts, surr_e3$contrasts)

	cont.sorted    = sort(contrasts, index.return=T)
	matches.sorted = matches[cont.sorted$ix]
  
	lower_bound = seq(1,         length(cont.sorted$x), by=n_per_bin)
	upper_bound = seq(n_per_bin, length(cont.sorted$x), by=n_per_bin)
  
	con_binned   = rep(0, length(lower_bound))
	match_binned = rep(0, length(lower_bound))
	sd_binned    = rep(0, length(lower_bound))
	n_binned     = rep(0, length(lower_bound))
  
	for (k in 1:length(lower_bound))
		{
		con_binned[k]   = mean(cont.sorted$x[lower_bound[k]:upper_bound[k]])
		match_binned[k] = mean(matches.sorted[lower_bound[k]:upper_bound[k]])
		n_binned[k]     = length(matches[lower_bound[k]:upper_bound[k]])
		sd_binned[k]    = sd(matches.sorted[lower_bound[k]:upper_bound[k]])
		}
	return(list(contr=con_binned, match=match_binned, sd=sd_binned, n=n_binned, sem=sd_binned/sqrt(n_binned), or_matches = matches, or_contrasts=contrasts))
	}


## prepare dat for computing global fit 'pretending' all data came from the same population
make_single_data=function(agg_matches)
    {
    common_data = as.data.frame(matrix(NA, nrow=60, ncol=3, dimnames=list(c(1:60), c('medium', 'contrast', 'match'))))

    count = 1
    for (my_medium in contexts)
        {
        idx = (count-1)*10
        common_data$medium[  (idx+1):(idx+10)] = rep(my_medium, 10)
        common_data$contrast[(idx+1):(idx+10)] = agg_matches[[my_medium]]$contr
        common_data$match[   (idx+1):(idx+10)] = agg_matches[[my_medium]]$match
        count = count + 1
        }
    return(common_data)
    }

my_orange    = '#FDAE61'
my_red       = '#D7191C'
my_darkblue  = '#2C7BB6'
my_lightblue = '#ABD9E9'

all_cols = list(plain='black', shadow='grey', transp_med=my_red, transp_med_steep=my_orange, transp_dark=my_darkblue, transp_dark_steep=my_lightblue)

contexts = c('plain', 'shadow', 'transp_med', 'transp_med_steep', 'transp_dark', 'transp_dark_steep')
observers= c('if','js','sk','wd')

global_rsq = matrix(NA, nrow=4, ncol=1, dimnames=list(observers, 'R_sq'))
slopes     = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
intercepts = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
rsq        = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))


svg('figs/michelson.svg', w=18,h=5, pointsize=8)
par(mfrow=c(1,4), lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=3, font=1, font.axis=3, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))


for (obs in observers)
    {
    ## combine log and design data
    log_data = read.table(paste('data/log_', obs, '.txt', sep=''),         header=TRUE)
    int_data = read.table(paste('data/intensities_', obs, '.txt', sep=''), header=TRUE)
    data = cbind(log_data, int_data)

    bin_matches = list()

    for (my_medium in contexts)
        {
        pl.f2=compute_contrast(compute_surround(data, my_medium,  'f2'))
        pl.f3=compute_contrast(compute_surround(data, my_medium,  'f3'))
        pl.e3=compute_contrast(compute_surround(data, my_medium,  'e3'))
        
        bin_matches[[my_medium]] = bin_eq_contrast(pl.f2, pl.f3, pl.e3, 12)
        }

    single_data = make_single_data(bin_matches)
    single.fit  = summary(lm(single_data$match ~ single_data$contrast))
    global_rsq[obs,] = single.fit$r.squared

    plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-1.05,1.05), ylim=c(0,420), xlab='local contrast', ylab='', xaxt='n', yaxt='n')

    axis(1, at=c(-1,-0.5,0,0.5,1),    cex=1.5)
    axis(2, at=c(0,100,200,300,400),  cex=1.5)
    
    if (obs =='if') {mtext(2, text='match luminance', line=2, cex=2.2)}

    for (curr_med in contexts)
        {
        curr.fit = lm(bin_matches[[curr_med]]$match ~ bin_matches[[curr_med]]$contr)
        slopes[obs, curr_med]     = summary(curr.fit)$coefficients[2]
        intercepts[obs, curr_med] = summary(curr.fit)$coefficients[1]
        rsq[obs, curr_med]        = summary(curr.fit)$r.squared

        abline(curr.fit, col= unlist(all_cols[curr_med]), lwd=2)

        arrows(bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match - bin_matches[[curr_med]]$sem), bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match + bin_matches[[curr_med]]$sem), length=0.01, angle=90, code=3)

        points(bin_matches[[curr_med]]$contr, bin_matches[[curr_med]]$match, pch=21, bg=unlist(all_cols[curr_med]))
        }
    legend('topleft', c(obs), bty='n', cex=2)
    }

legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=0.9, cex=1.2)

dev.off()

write.table(format(global_rsq, digits=4), 'stat_out/michelson_global_rsq.txt',quote=FALSE)
write.table(format(slopes, digits=4),     'stat_out/michelson_slopes.txt',    quote=FALSE)
write.table(format(intercepts, digits=4), 'stat_out/michelson_intercepts.txt',quote=FALSE)
write.table(format(rsq, digits=4),        'stat_out/michelson_rsq.txt',       quote=FALSE)



