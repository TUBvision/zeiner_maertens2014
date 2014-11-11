rm(list=ls())

## read out target intensity and surround check intensity in HRL units
compute_surround <- function(data, n_surround=8)
    {
    ## FOUR adjacent checks
    surr_checks = list(e3=c('e2','d3','e4','f3'), f2=c('f1','e2','f3','g2'), f3=c('f2','e3','f4','g3'))
    
    ## EIGHT adjacent checks
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


## MICHELSON CONTRASTS
michelson <- function(data)
    {
    ## adds contrast column to data frame
    data$contrast = (data$targ_int - data$surr_int)/(data$targ_int + data$surr_int)
    return(data)
    }


## SCALED CONTRASTS
## scale input range to range specified by min_out and max_out
scaling=function(input, min_out, max_out)
    {
    min_in = min(input)
    max_in = max(input)
    range_in = max_in  - min_in
    range_out= max_out - min_out
    return( (input-min_in)/range_in * range_out + min_out )
    }


scaled_contrast <- function(data)
    {
    ## adds scaled_contrast column to data frame
    ## min and max of PLAIN contrasts for normalization purpose
    min_pl = min(data$contrast[data$medium=='plain'])
    max_pl = max(data$contrast[data$medium=='plain'])
    
    data$scaled_contrast = rep(0, dim(data)[1])
    
    ## SCALE contrasts in each framework with those in PLAIN
    for (curr_med in unique(data$medium))
        {
        data$scaled_contrast[data$medium==curr_med] = scaling(data$contrast[data$medium==curr_med], min_pl, max_pl)
        }
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

# ## EDGE INTEGRATION
# edge_integration <- function(data, w_near=1, w_far=0.79)
#     {
#     ## add edge column to data frame, default weights for incremental targets
#     ## w_n * log(D/A) + w_f * log(A/B), whereby background intensity is the average plain check luminance
#     bg_int = mean(data$targ_int[data$medium =='plain'])
#     
#     ## edge integration computation after Rudd, 2013 p. 9, eq. 3a and differential weighting for near and far weights, ratios p.11
#     data$edge = w_near * log(data$targ_int/data$surr_int) + w_far * log(data$surr_int/bg_int)
#     return(data)
#     }


## ANCHORING
scale_normalization=function(l_uncorrected, version='new')
    {
    ## apply to predicted reflectances if their resulting range (max:min) is smaller than 30
    l_max = max(l_uncorrected)
    l_min = min(l_uncorrected)

    if (version=='old')
        {
        l_range = l_max/l_min
        l_corrected = 1 / (1+(.56 * ( log(30) - log(l_range) ))) * l_uncorrected
        }
    if (version=='new')
        {
        l_corrected = 0.44 * l_uncorrected - 0.44 * l_max + (0.82 * (l_uncorrected-l_min)/(l_max-l_min) - 0.87)
        }
    return(l_corrected)
    }


anchor_to_white=function(data, local_weight=1, model_version='old')
    {
    ## 2. determine max luminances (in units of graphic card intensity) in each framework
    data$R_framework = rep(NA, dim(data)[1])

    for (my_medium in unique(data$medium))
        {
        cdata     = subset(data, medium == my_medium)
        max_local = max(unique(cdata$targ_int))
        max_global= max(unique(data$targ_int[data$medium=='plain']))

    ## 3. apply highest luminance rule
        data$R_framework[data$medium==my_medium] = local_weight * (cdata$targ_int/max_local * 0.9) + (1-local_weight) * (cdata$targ_int/max_global * 0.9)
    
    ## 4. apply scale normalization if necessary
        cdata     = subset(data, medium == my_medium)
        if (max(cdata$R_framework)/min(cdata$R_framework) < 30)
            {
            data$R_framework[data$medium == my_medium] = scale_normalization(cdata$R_framework, version=model_version)
            }
        }
    return(data)
    }


## summarize respective dependent variable into bins
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
    return(list(contr=con_binned, match=match_binned, sd=sd_binned, n=n_binned, sem=sd_binned/sqrt(n_binned)))
    }


aggregate_subset<-function(subdata, model)
    {
    ## aggregate data for specified dependent variable (model), and rows in subdata
    bin_matches = list()
    
    for (curr_med in  unique(subdata$medium))
        {
        x = subdata[, model][subdata$medium==curr_med]
        y = subdata$match_lum[subdata$medium==curr_med]
        bin_matches[[curr_med]] = bin_eq_n(x,y)
        if (model == 'target_lum')
            {
            bin_matches[[curr_med]]$contr = atf(bin_matches[[curr_med]]$contr, curr_med)
            }
        }
    return(bin_matches)
    }


## compute global fit 'pretending' all data came from the same population
global_fit <- function(agg_matches)
    {
    ## rearrange the list of aggregated matches into data frame in order to perform the global fit
    single_data = as.data.frame(matrix(NA, nrow=60, ncol=3, dimnames=list(c(1:60), c('medium', 'contrast', 'match'))))

    count = 1
    for (my_medium in names(agg_matches))
        {
        idx = (count-1)*10
        single_data$medium[  (idx+1):(idx+10)] = rep(my_medium, 10)
        single_data$contrast[(idx+1):(idx+10)] = agg_matches[[my_medium]]$contr
        single_data$match[   (idx+1):(idx+10)] = agg_matches[[my_medium]]$match
        count = count + 1
        }
    single.fit = summary(lm(single_data$match ~ single_data$contrast))
    return(single.fit)
    }


## compute deviation from global fit, separately for each condition
deviation_from_global_fit <- function(agg_matches, single.fit)
    {
    cond.deviation = list()

    for (my_medium in names(agg_matches))
        {
        new.x   = agg_matches[[my_medium]]$contr
        new.y   = agg_matches[[my_medium]]$match
        
        ## 1. predict y-values for condition x-values using coefficients from global fit
        pred.y  = coef(single.fit)[2] * new.x + coef(single.fit)[1]
        
        ## 2. compute variance of y and residual variance
        res.var = sum((new.y-pred.y)^2)
        y.var   = sum((new.y - mean(new.y))^2)

        ## 3. compute fraction of variance in condition explained by global fit to total variance in condition
        cond.deviation[[my_medium]] = res.var/y.var
        }
    return(cond.deviation)
    }



## compute prediction benefit
prediction_benefit <- function(agg_matches)
    {
    ## 1. linear fit to plain condition as reference condition
    plain.lm = lm(agg_matches[['plain']]$match ~ agg_matches[['plain']]$contr)

    pred.benefit = list()
    
    for (my_medium in names(agg_matches))
        {
        new.x   = agg_matches[[my_medium]]$contr
        new.y   = agg_matches[[my_medium]]$match
        
        ## 2. predict y-values for new x-values using coefficients from plain fit
        pred.y  = coef(plain.lm)[2] * new.x + coef(plain.lm)[1]
        
        ## 3. compute deviance residuals between actual and predicted y
        res.var = sum((new.y-pred.y)^2)
        
        ## 4. computed variance of y
        y.var   = sum((new.y - mean(new.y))^2)
        pred.benefit[[my_medium]] = res.var/y.var
        }
    return(pred.benefit)
    }


## color definition for plotting
my_orange    = '#FDAE61'
my_red       = '#D7191C'
my_darkblue  = '#2C7BB6'
my_lightblue = '#ABD9E9'

all_cols = list(plain='black', shadow='grey', transp_dark=my_darkblue, transp_dark_steep=my_lightblue, transp_med=my_red, transp_med_steep=my_orange)

mean_bg = 150 ## mean of plain checkerboard luminances

observers = c('if', 'js', 'sk', 'wd')
models    = c('contrast', 'edge', 'R_framework', 'scaled_contrast')
model_name= c('contrast', 'edge integration', 'anchoring', 'normalized contrast')

## output variables
global_rsq = array(NA, c(length(models),4), dimnames=list(models=models, obs=observers))

pred_benefit = list()
export_pred  = array(NA, c(length(models),6,4), dimnames=list(models=models, medium=names(all_cols), obs=observers))

for (obs in observers)
    {
    ## combine log and design data
    log_data = read.table(paste('data/log_', obs, '.txt', sep=''),         header=TRUE)
    int_data = read.table(paste('data/intensities_', obs, '.txt', sep=''), header=TRUE)
    data = cbind(log_data, int_data)

    ## add target and surround intensity columns
    data = compute_surround(data)
    ## compute michelson contrast
    data = michelson(data)
    ## compute scaled contrasts
    data = scaled_contrast(data)
    ## edge integration
    data = edge_integration(data)
    ## anchoring
    data = anchor_to_white(data, local_weight=0.5)

    dev_from_global = list()

    for (curr_model in models)
        {
        bin_matches = aggregate_subset(data, curr_model)
        pred_benefit[[curr_model]] = prediction_benefit(bin_matches)

        global.fit = global_fit(bin_matches)
        global_rsq[curr_model, obs] = global.fit$r.squared
        dev_from_global[[curr_model]] = deviation_from_global_fit(bin_matches, global.fit)
        }

    for (my_medium in colnames(export_pred))
        {
        for (my_model in rownames(export_pred))
            {
            export_pred[my_model, my_medium, obs] = dev_from_global[[my_model]][[my_medium]]
            }
        }
    }

m.dev_from_global   = apply(export_pred, c(1,2), mean)
sem.dev_from_global = apply(export_pred, c(1,2), sd)/sqrt(2)

svg('figs/model_comparison_anch_0.5.svg', w=7,h=5, pointsize=8)

par(lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.2, cex.lab=1.5, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,4,2,1), mgp=c(2,1,0))

barplot(m.dev_from_global, beside=T, names.arg=c('plain','shadow','dark transp\nlow', 'dark transp\nhigh', 'light transp\nlow', 'light transp\nhigh'), col=c(my_orange, my_darkblue, my_lightblue, my_red), ylim=c(0,1), ylab="residual variance/total variance")
legend("topleft", model_name, bty="n", fill=c(my_orange, my_darkblue, my_lightblue, my_red))
dev.off()





















