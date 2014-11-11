rm(list=ls())


compute_mean_surround <- function(data, n_surround=4)
    {
    ## contrast based on FOUR adjacent checks
    surr_checks = list(e3=c('e2','d3','e4','f3'), f2=c('f1','e2','f3','g2'), f3=c('f2','e3','f4','g3'))
    
    if (n_surround == 8)
        {
        surr_checks = list(e3=c('e2','d3','e4','f3', 'd2','d4','f2','f4'), f2=c('f1','e2','f3','g2', 'e1','e3','g1','g3'), f3=c('f2','e3','f4','g3', 'e2','e4','g2','g4'))
        }

    surr_int = matrix(0, nrow=dim(data)[1], length(surr_checks$e3))
    targ_int = rep(0, dim(data)[1])
    
    ## read out target and surround intensities from corresponding columns
    for (idx in 1:dim(data)[1])
        {
        line          = data[idx,]
        curr_surr     = unlist(surr_checks[line$location])
        targ_int[idx] = as.numeric(line[paste(as.character(line$location), '_v', sep='')])

        surr_count = 1
        for (s in curr_surr)
            {
            surr_int[idx, surr_count]=as.numeric(line[paste(s, '_v', sep='')])
            surr_count = surr_count+1
            }
        }
    data$targ_int = targ_int
    data$surr_int = apply(surr_int, 1, mean)
    return(data)
    }


## EDGE INTEGRATION
edge_integration <- function(data, w_near=1)
    {
    ## add edge column to data frame, default weights for incremental targets
    ## w_n * log(D/A) + w_f * log(A/B), whereby background intensity is the average plain check luminance
    bg_int    = mean(data$targ_int[data$medium =='plain'])
    data$edge = rep(0, dim(data)[1])
    
    ## edge integration computation after Rudd, 2013 p. 9, eq. 3a and differential weighting for near and far weights, ratios p.11
    for (k in 1:length(data$targ_int))
        {
        w_far = 0.3
        if (data$targ_int[k]/data$surr_int[k] > 1)  {w_far = 0.79}
        
        data$edge[k] = w_near * log(data$targ_int[k]/data$surr_int[k]) + w_far * log(data$surr_int[k]/bg_int)
        }
    
    return(data)
    }


# edge_integration <- function(targ_int, surr_int)
# ## OLD function used for first submission to JOV, errors: bg_int specified in luminance instead of HRL units as the targ_int and surr_int variables, all targets were considered as increments here
#     {
#     w_far  = 1
#     bg_int = 120
#     w_near = rep(1.26, length(targ_int))
#     ## edge integration computation after Rudd, 2013 p. 9, eq. 3a and differential weighting for near and far weights, ratios p.11
#     w_near * log(targ_int/surr_int) + w_far * log(surr_int/bg_int)
# #     w_near * (targ_int/surr_int) + w_far * (surr_int/bg_int)
#     }

bin_eq_n<-function(x,y, n_per_bin=12)
    {
    ## binning of scaled contrast values, input format is list that contains scaled data from all contexts (medium)
    x.sorted       = sort(x, index.return=T)
    matches.sorted = y[x.sorted$ix]
        
    lower_bound = seq(1,         length(x.sorted$x), by=n_per_bin)
    upper_bound = seq(n_per_bin, length(x.sorted$x), by=n_per_bin)
        
    con_binned   = rep(0, length(lower_bound))
    match_binned = rep(0, length(lower_bound))
    sd_binned    = rep(0, length(lower_bound))
    n_binned     = rep(0, length(lower_bound))
        
    for (k in 1:length(lower_bound))
        {
        con_binned[k]   = mean(x.sorted$x[lower_bound[k]:upper_bound[k]])
        match_binned[k] = mean(matches.sorted[lower_bound[k]:upper_bound[k]])
        n_binned[k]     = length(matches.sorted[lower_bound[k]:upper_bound[k]])
        sd_binned[k]    = sd(matches.sorted[lower_bound[k]:upper_bound[k]])
        }
    return(list(contr=con_binned, matches=match_binned, sd=sd_binned, n=n_binned, sem=sd_binned/sqrt(n_binned)))
    }


log_bin_eq_n<-function(x,y, n_per_bin=12)
    {
    ## binning of scaled contrast values, input format is list that contains scaled data from all contexts (medium)
    y = log(y)
    y[is.infinite(y)]=NA
    x.sorted       = sort(x, index.return=T)
    matches.sorted = y[x.sorted$ix]
        
    lower_bound = seq(1,         length(x.sorted$x), by=n_per_bin)
    upper_bound = seq(n_per_bin, length(x.sorted$x), by=n_per_bin)
        
    con_binned   = rep(0, length(lower_bound))
    match_binned = rep(0, length(lower_bound))
    sd_binned    = rep(0, length(lower_bound))
    n_binned     = rep(0, length(lower_bound))
        
    for (k in 1:length(lower_bound))
        {
        con_binned[k]   = mean(x.sorted$x[lower_bound[k]:upper_bound[k]])
        match_binned[k] = mean(matches.sorted[lower_bound[k]:upper_bound[k]])
        n_binned[k]     = length(matches.sorted[lower_bound[k]:upper_bound[k]])
        sd_binned[k]    = sd(matches.sorted[lower_bound[k]:upper_bound[k]])
        }
    return(list(contr=con_binned, matches=match_binned, sd=sd_binned, n=n_binned, sem=sd_binned/sqrt(n_binned)))
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

obs = 'js'
mean_bg = 120

## combine log and design data
log_data = read.table(paste('../log_', obs, '.txt', sep=''),         header=TRUE)
int_data = read.table(paste('../intensities_', obs, '.txt', sep=''), header=TRUE)
data = cbind(log_data, int_data)

data = compute_mean_surround(data)
data = edge_integration(data)


bin_matches = list()

for (curr_med in  unique(data$medium))
    {
    x = data$edge[data$medium==curr_med]
    y = data$match_lum[data$medium==curr_med]

    bin_matches[[curr_med]] = bin_eq_n(x,y)
    }

## plotting
X11()
par(lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-3.4,1.4), ylim=c(0,420), xlab='integrated edge', ylab='', xaxt='n', yaxt='n')

axis(1, at=c(-10:-4),    cex=1.5)
axis(2, at=c(0,100,200,300,400),  cex=1.5)

for (curr_med in unique(data$medium) )
    {
    arrows(bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match - bin_matches[[curr_med]]$sem), bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match + bin_matches[[curr_med]]$sem), length=0.01, angle=90, code=3)

    points(bin_matches[[curr_med]]$contr, bin_matches[[curr_med]]$match, pch=21, bg=unlist(all_cols[curr_med]))
    }


## plotting in log coordinates
global_rsq = matrix(NA, nrow=4, ncol=1, dimnames=list(observers, 'R_sq'))
slopes     = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
intercepts = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
rsq        = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))

# png('figs/rudd.png', w=1680,h=450)
svg('figs/rudd.svg', w=18,h=5, pointsize=8)
par(mfrow=c(1,4), lwd=1.5, tck=.01, bty='n', cex=1.5, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

for (obs in observers)
    {
    ## combine log and design data
    log_data = read.table(paste('../log_', obs, '.txt', sep=''),         header=TRUE)
    int_data = read.table(paste('../intensities_', obs, '.txt', sep=''), header=TRUE)
    data = cbind(log_data, int_data)

    data = compute_mean_surround(data)
#     data$edge = edge_integration(data$targ_int, data$surr_int) old edge_integration function
    data = edge_integration(data)

    bin_matches = list()

    for (curr_med in  unique(data$medium))
        {
        x = data$edge[data$medium==curr_med]
        y = data$match_lum[data$medium==curr_med]

        bin_matches[[curr_med]] = log_bin_eq_n(x,y)
        }

    single_data = make_single_data(bin_matches)
    single.fit  = summary(lm(single_data$match ~ single_data$contrast))
    global_rsq[obs,] = single.fit$r.squared

    plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-3.4,1.4), ylim=c(2,6), xlab='integrated edge', ylab='', xaxt='n', yaxt='n')
#     plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-10,-4), ylim=c(2,6), xlab='integrated edge', ylab='', xaxt='n', yaxt='n')

#     axis(1, at=c(-10:-4), cex=1.5)
    axis(1, at=seq(-3.4,1.4, length.out=5),    cex=1.5)
    axis(2, at=c(2:6),    cex=1.5)

    if (obs =='if') {mtext(2, text='log match luminance', line=2, cex=2.2)}

    for (curr_med in unique(data$medium) )
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

# write.table(format(global_rsq, digits=4), 'stat_out/rudd_global_rsq.txt',quote=FALSE)
# write.table(format(slopes, digits=4),     'stat_out/rudd_slopes.txt',    quote=FALSE)
# write.table(format(intercepts, digits=4), 'stat_out/rudd_intercepts.txt',quote=FALSE)
# write.table(format(rsq, digits=4),        'stat_out/rudd_rsq.txt',       quote=FALSE)


## single observer plot
obs = 'js'

png(paste('figs/', obs, '_rudd.png', sep=''), w=600,h=550)
par(lwd=2, tck=.01, bty='n', cex.axis=1.5, cex.lab=1.5, cex.main=1.5, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,1,1), mgp=c(2,0.5,0),  cex=1.5)

## combine log and design data
log_data = read.table(paste('../log_', obs, '.txt', sep=''),         header=TRUE)
int_data = read.table(paste('../intensities_', obs, '.txt', sep=''), header=TRUE)
data = cbind(log_data, int_data)

data = compute_mean_surround(data)
#     data$edge = edge_integration(data$targ_int, data$surr_int) old edge_integration function
data = edge_integration(data)

bin_matches = list()

for (curr_med in  unique(data$medium))
    {
    x = data$edge[data$medium==curr_med]
    y = data$match_lum[data$medium==curr_med]

    bin_matches[[curr_med]] = log_bin_eq_n(x,y)
    }

single_data = make_single_data(bin_matches)
single.fit  = summary(lm(single_data$match ~ single_data$contrast))
global_rsq[obs,] = single.fit$r.squared

plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-3.4,1.4), ylim=c(2,6), xlab='integrated edge', ylab='', xaxt='n', yaxt='n')
#     plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-10,-4), ylim=c(2,6), xlab='integrated edge', ylab='', xaxt='n', yaxt='n')

#     axis(1, at=c(-10:-4), cex=1.5)
axis(1, at=seq(-3.4,1.4, length.out=5),    cex=1.5)
axis(2, at=c(2:6),    cex=1.5)

mtext(2, text='log match luminance', line=2, cex=2.2)

for (curr_med in unique(data$medium) )
    {
    curr.fit = lm(bin_matches[[curr_med]]$match ~ bin_matches[[curr_med]]$contr)
    slopes[obs, curr_med]     = summary(curr.fit)$coefficients[2]
    intercepts[obs, curr_med] = summary(curr.fit)$coefficients[1]
    rsq[obs, curr_med]        = summary(curr.fit)$r.squared

    abline(curr.fit, col= unlist(all_cols[curr_med]), lwd=2)

    arrows(bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match - bin_matches[[curr_med]]$sem), bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match + bin_matches[[curr_med]]$sem), length=0.01, angle=90, code=3)

    points(bin_matches[[curr_med]]$contr, bin_matches[[curr_med]]$match, pch=21, bg=unlist(all_cols[curr_med]))
    }
legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=0.9, cex=1.2)

dev.off()


## NON log units
# png('figs/rudd.png', w=1680,h=450)
# X11()
# par(mfrow=c(1,4), lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))
# 
# for (obs in observers)
#     {
#     ## combine log and design data
#     log_data = read.table(paste('../log_', obs, '.txt', sep=''),         header=TRUE)
#     int_data = read.table(paste('../intensities_', obs, '.txt', sep=''), header=TRUE)
#     data = cbind(log_data, int_data)
# 
#     data = compute_mean_surround(data)
#     data = edge_integration(data)
# 
#     bin_matches = list()
# 
#     for (curr_med in  unique(data$medium))
#         {
#         x = data$edge[data$medium==curr_med]
#         y = data$match_lum[data$medium==curr_med]
# 
#         bin_matches[[curr_med]] = bin_eq_n(x,y)
#         }
# 
#     single_data = make_single_data(bin_matches)
#     single.fit  = summary(lm(single_data$match ~ single_data$contrast))
#     global_rsq[obs,] = single.fit$r.squared
# 
#     plot(bin_matches$plain$contr, bin_matches$plain$match, pch=21, xlim=c(-3,1.4), ylim=c(0,420), xlab='local contrast', ylab='', xaxt='n', yaxt='n')
# 
#     axis(1, at=c(-10:-4), cex=1.5)
#     axis(2, at=c(0,100,200,300,400),  cex=1.5)
# 
#     if (obs =='if') {mtext(2, text='log match luminance', line=2, cex=2.2)}
# 
#     for (curr_med in unique(data$medium) )
#         {
#         curr.fit = lm(bin_matches[[curr_med]]$match ~ bin_matches[[curr_med]]$contr)
#         slopes[obs, curr_med]     = summary(curr.fit)$coefficients[2]
#         intercepts[obs, curr_med] = summary(curr.fit)$coefficients[1]
#         rsq[obs, curr_med]        = summary(curr.fit)$r.squared
# 
#         abline(curr.fit, col= unlist(all_cols[curr_med]), lwd=2)
# 
#         arrows(bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match - bin_matches[[curr_med]]$sem), bin_matches[[curr_med]]$contr, (bin_matches[[curr_med]]$match + bin_matches[[curr_med]]$sem), length=0.01, angle=90, code=3)
# 
#         points(bin_matches[[curr_med]]$contr, bin_matches[[curr_med]]$match, pch=21, bg=unlist(all_cols[curr_med]))
#         }
#     legend('topleft', c(obs), bty='n', cex=1.5)
#     }
# 
# legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=0.9, cex=1.2)

# dev.off()
# 
# write.table(format(global_rsq, digits=4), 'stat_out/rudd_global_rsq.txt',quote=FALSE)
# write.table(format(slopes, digits=4),     'stat_out/rudd_slopes.txt',    quote=FALSE)
# write.table(format(intercepts, digits=4), 'stat_out/rudd_intercepts.txt',quote=FALSE)
# write.table(format(rsq, digits=4),        'stat_out/rudd_rsq.txt',       quote=FALSE)
# 
