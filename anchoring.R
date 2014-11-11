rm(list=ls())

## anchoring in local framework
anchor_to_white_simple = function(data)
    {
    ## compute reflectances within local framework by dividing all luminances by the highest luminance in framework
    ## r_i = l_i/l_max * 0.9
    data$R_framework = rep(NA, dim(data)[1])

    for (my_medium in unique(data$medium))
        {
        cdata = subset(data, medium == my_medium)

        ## highest luminance rule
        data$R_framework[data$medium==my_medium] = cdata$targ_int/max(cdata$targ_int) * 0.9
        }
    return(data)
    }


anchor_to_white_framework = function(data, weight_local, max_local)
    {
    ## compute reflectances within local framework by dividing all luminances by the highest luminance in framework
    ## r_i = l_i/l_max * 0.9
    data$R_framework = rep(NA, dim(data)[1])

    for (my_medium in unique(data$medium))
        {
        cdata = subset(data, medium == my_medium)
        curr_max_local = max_local[[my_medium]]
        max_global     = max_local[['plain']]

        ## highest luminance rule
        data$R_framework[data$medium==my_medium] = weight_local * (cdata$targ_int/curr_max_local * 0.9) + (1-weight_local) * (cdata$targ_int/max_global * 0.9)
        }
    return(data)
    }

## scale normalization after Gilchrist
perceived_lum = function(l_uncorrected)
    {
    l_range = max(l_uncorrected)/min(l_uncorrected)
    k = 1 / (1+(.56 * ( log(30) - log(l_range) ))) * l_uncorrected
    return(k)
    }


normalize_scale = function(data)
    {
    ## scale normalization applied to anchored reflectances (relative luminances)
    data$R_framework_scaled = data$R_framework

    for (my_medium in unique(data$medium))
        {
        cdata = subset(data, medium == my_medium)
        ## check whether scale is smaller than 30:1
        if (max(cdata$R_framework)/min(cdata$R_framework) < 30)
            {
            data$R_framework_scaled[data$medium == my_medium] = perceived_lum(cdata$R_framework)
            }
        }
    return(data)
    }


bin_equal_n = function(x, y, n_per_bin=12)
    {
    ## sorts values in variable x and groups them into bins of n, values in y are grouped correspondingly, returns summary statistics for the bins
    x.sorted = sort(x, index.return=T)
    y.sorted = y[x.sorted$ix]
  
    lower_bound = seq(1,         length(x.sorted$x), by=n_per_bin)
    upper_bound = seq(n_per_bin, length(x.sorted$x), by=n_per_bin)
  
    x.binned = rep(0, length(lower_bound))
    y.binned = rep(0, length(lower_bound))
    y.binned.sd = rep(0, length(lower_bound))
    y.binned.n = rep(0, length(lower_bound))
  
    for (k in 1:length(lower_bound))
        {
        x.binned[k]   = mean(x.sorted$x[lower_bound[k]:upper_bound[k]])
        y.binned[k]   = mean(y.sorted[lower_bound[k]:upper_bound[k]])
        y.binned.n[k] = length(y[lower_bound[k]:upper_bound[k]])
        y.binned.sd[k]= sd(y.sorted[lower_bound[k]:upper_bound[k]])
        }

    list(contr=x.binned, 
         match=y.binned, 
         sd=y.binned.sd, 
         n=y.binned.sd, 
         sem=y.binned.sd/sqrt(y.binned.n))
    }

bin_given_x = function(x, y)
    {
    grouped_data = tapply(y, x, mean)

    list(contr=as.numeric(names(grouped_data)), 
         match=grouped_data)
    }


## rearrange data to compute global fit 'pretending' all data came from the same population
make_single_data=function(agg_matches)
    {
    common_data = as.data.frame(matrix(NA, nrow=60, ncol=3, dimnames=list(c(1:60), c('medium', 'contrast', 'match'))))

    count = 1
    for (my_medium in names(agg_matches))
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

all_cols = list(plain = 'black', shadow='grey', transp_med=my_red, transp_med_steep=my_orange, transp_dark=my_darkblue, transp_dark_steep=my_lightblue)


local_weight = 0.5

contexts = c('plain', 'shadow', 'transp_med', 'transp_med_steep', 'transp_dark', 'transp_dark_steep')
observers= c('if','js','sk','wd')

global_rsq = matrix(NA, nrow=4, ncol=1, dimnames=list(observers, 'R_sq') )
slopes     = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
intercepts = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
rsq        = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))


svg('figs/anchoring.svg', w=18,h=5, pointsize=8)
par(mfrow=c(1,4), lwd=1.5, tck=.01, bty='n', cex=1.5, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

for (obs in observers)
    {
    ## combine log data and framework information
    log_data = read.table(paste('data/log_', obs, '.txt', sep=''),header=TRUE)
    int_data = read.table(paste('data/intensities_', obs, '.txt', sep=''), header=TRUE)
    fr_data  = read.table(paste('data/framework_intensity_', obs, '.txt', sep=''), header=TRUE)

    data = cbind(log_data, int_data, fr_data[3:8])

    ## 1. determine gray value, and luminance in HRL units at target position
    data$targ_int = rep(0,dim(data)[1])
    data$targ_gv  = rep(0,dim(data)[1])

    for (idx in 1:dim(data)[1])
        {
        line = data[idx,]
        ## read out gray value at target location, intensity [0,1], gray value [0, 255]
        data$targ_int[idx] = as.numeric(line[paste(as.character(line$location), '_v', sep='')])
        }

    ## determine max luminances (in units of graphic card intensity) in each framework
    max_local = list()
    for (my_medium in unique(data$medium))
        {
        max_local[my_medium] = max(unique(data$targ_int[data$medium==my_medium]))
        }

    ## 2. compute relative luminances == framework reflectance
    data = anchor_to_white_framework(data, local_weight, max_local)

    ## 3. normalize resulting reflectances when range < 30:1
    data = normalize_scale(data)

    anchor_matches = list()

    ## aggregate values by summarizing the luminance matches according to the newly computed framework lightness variable
    for (my_medium in unique(data$medium))
        {
        ## binning with 12 observation per group
        anchor_matches[[my_medium]] = bin_equal_n(data$R_framework_scaled[data$medium==my_medium], data$match_lum[data$medium==my_medium], 12)
        }

    ## COMPUTE GLOBAL FIT 'pretending' all data came from one function
    single_data = make_single_data(anchor_matches)
    single.fit  = summary(lm(single_data$match ~ single_data$contrast))
    global_rsq[obs,] = single.fit$r.squared


    plot(anchor_matches$plain$contr, anchor_matches$plain$match, pch=21, xlim=c(-.05,1.05), ylim=c(0,420), xlab='anchoring lightness', ylab='', xaxt='n', yaxt='n')

    axis(1, at=seq(0,1, length.out=5), cex=1.5)
    axis(2, at=c(0,100,200,300,400),   cex=1.5)

    if (obs =='if') {mtext(2, text='match luminance', line=2, cex=2.2)}

    for (my_medium in unique(data$medium))
        {
        curr.fit = lm(anchor_matches[[my_medium]]$match ~ anchor_matches[[my_medium]]$contr)
        slopes[obs, my_medium]     = summary(curr.fit)$coefficients[2]
        intercepts[obs, my_medium] = summary(curr.fit)$coefficients[1]
        rsq[obs, my_medium]        = summary(curr.fit)$r.squared

        abline(curr.fit, col= unlist(all_cols[my_medium]), lwd=2)

        arrows(anchor_matches[[my_medium]]$contr, (anchor_matches[[my_medium]]$match-anchor_matches[[my_medium]]$sem), anchor_matches[[my_medium]]$contr, (anchor_matches[[my_medium]]$match+anchor_matches[[my_medium]]$sem), length=0.01, angle=90, code=3)

        points(anchor_matches[[my_medium]]$contr, anchor_matches[[my_medium]]$match, pch=21, bg=unlist(all_cols[my_medium]))
        }
    legend('topleft', c(obs), bty='n', cex=2)
    }

legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=0.9, cex=1.2)

dev.off()

write.table(format(global_rsq, dig=4), paste('stat_out/anchoring_global_rsq', local_weight,'.txt', sep=""), quote=F)
write.table(format(slopes, dig=4),     paste('stat_out/anchoring_slopes', local_weight,'.txt', sep=""), quote=F)
write.table(format(intercepts, dig=4), paste('stat_out/anchoring_intercepts', local_weight,'.txt', sep=""), quote=F)
write.table(format(rsq, dig=4),        paste('stat_out/anchoring_rsq', local_weight,'.txt', sep=""), quote=F)
