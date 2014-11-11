rm(list=ls())


## CONTRAST
## compute contrast directly on raw data
compute_rayleigh<-function(data)
    {
    ## contrast based on EIGHT adjacent checks
    surr_checks = list(e3=c('e2','d3','e4','f3', 'd2','d4','f2','f4'), f2=c('f1','e2','f3','g2', 'e1','e3','g1','g3'), f3=c('f2','e3','f4','g3', 'e2','e4','g2','g4') )
    
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
    data$rayleigh = (targ_int-apply(surr_int,1,mean))/(targ_int+apply(surr_int,1,mean))
    return(data)
    }

## concatenate lists of contrast and match values from different checks
concatenate_checks=function(surr_f2, surr_f3, surr_e3)
    {
    matches    = c(surr_f2$match    , surr_f3$match,     surr_e3$match)
    contrasts  = c(surr_f2$contrasts, surr_f3$contrasts, surr_e3$contrasts)

    return(list(or_matches = matches, or_contrasts=contrasts))
    }

## scale input range to range specified by min_plain and max_plain
scaling=function(input, min_plain, max_plain)
    {
    min_in = min(input)
    max_in = max(input)
    range_in = max_in  - min_in
    range_plain= max_plain - min_plain
    return( (input-min_in)/range_in * range_plain + min_plain )
    }

bin_eq_scaled_contrast<-function(all_contrasts, n_per_bin=12)
    {
    ## binning of scaled contrast values, input format is list that contains scaled data from all contexts (medium)
    for (curr_med in names(all_contrasts))
        {
        curr_matches    = all_contrasts[[curr_med]]$or_matches
        curr_contrasts  = all_contrasts[[curr_med]]$scaled_contrasts

        cont.sorted    = sort(curr_contrasts, index.return=T)
        matches.sorted = curr_matches[cont.sorted$ix]
        
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
            n_binned[k]     = length(matches.sorted[lower_bound[k]:upper_bound[k]])
            sd_binned[k]    = sd(matches.sorted[lower_bound[k]:upper_bound[k]])
            }
        all_contrasts[[curr_med]]$scaled_contr=con_binned
        all_contrasts[[curr_med]]$scaled_match=match_binned
        all_contrasts[[curr_med]]$scaled_sem=sd_binned/sqrt(n_binned)
        }
    return(all_contrasts)
    }


## prepare dat for computing global fit 'pretending' all data came from the same population
make_single_data_scaled=function(agg_matches)
    {
    common_data = as.data.frame(matrix(NA, nrow=60, ncol=3, dimnames=list(c(1:60), c('medium', 'contrast', 'match'))))

    count = 1
    for (my_medium in names(agg_matches))
        {
        idx = (count-1)*10
        common_data$medium[  (idx+1):(idx+10)] = rep(my_medium, 10)
        common_data$contrast[(idx+1):(idx+10)] = agg_matches[[my_medium]]$scaled_contr
        common_data$match[   (idx+1):(idx+10)] = agg_matches[[my_medium]]$scaled_match
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


global_rsq = matrix(NA,nr=4, nc=1, dimnames=list(observers, 'R_sq'))
slopes     = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
intercepts = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))
rsq        = matrix(0, nr=4, nc=6, dimn=list(observers, names(all_cols)))

svg('figs/normalized_contrast.svg', w=18,h=5, pointsize=8)
par(mfrow=c(1,4), lwd=1.5, tck=.01, bty='n', cex=1.5, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=1, font=1, font.axis=1, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

for (obs in observers)
    {
    ## combine log and design data
    log_data= read.table(paste('data/log_', obs, '.txt', sep=''),         header=TRUE)
    int_data= read.table(paste('data/intensities_', obs, '.txt', sep=''), header=TRUE)
    data_in = cbind(log_data, int_data)

    ## 1. COMPUTE RAYLEIGH contrast based on 8 adjacent checks (edge + corner)
    data    = compute_rayleigh(data_in)
    
    contrasts = list()

    for (curr_med in unique(data$medium))
        {
        f2 = list(contrasts = data$rayleigh[data$medium == curr_med & data$location=='f2'], match = data$match_lum[data$medium == curr_med & data$location=='f2'])
        f3 = list(contrasts = data$rayleigh[data$medium == curr_med & data$location=='f3'], match = data$match_lum[data$medium == curr_med & data$location=='f3'])
        e3 = list(contrasts = data$rayleigh[data$medium == curr_med & data$location=='e3'], match = data$match_lum[data$medium == curr_med & data$location=='e3'])
        
        ## 2. summarize lists of matches and contrast from different checks
        contrasts[[curr_med]] = concatenate_checks(f2, f3, e3)
        }
        
        ## min and max of PLAIN contrasts for normalization purpose
        min_pl = min(contrasts[['plain']]$or_contrasts)
        max_pl = max(contrasts[['plain']]$or_contrasts)
    
    ## 3. SCALE contrasts in each framework with those in PLAIN
    for (curr_med in unique(data$medium))
        {
        contrasts[[curr_med]]$scaled_contrasts = scaling(contrasts[[curr_med]]$or_contrasts, min_pl, max_pl)
        }
    ## 4. SUMMARIZE scaled contrast values into 10 bins of 12 values each
    bin_matches = bin_eq_scaled_contrast(contrasts)

    ## COMPUTE GLOBAL FIT 'pretending' all data came from one function
    single_data = make_single_data_scaled(bin_matches)
    single.fit  = summary(lm(single_data$match ~ single_data$contrast))
    global_rsq[obs,] = single.fit$r.squared

    plot(bin_matches$plain$scaled_contr, bin_matches$plain$scaled_match, pch=21, xlim=c(-1.05,1.05), ylim=c(0,420), xlab='normalized contrast', ylab='', xaxt='n', yaxt='n')

    axis(1, at=c(-1,-0.5,0,0.5,1),    cex=1.5)
    axis(2, at=c(0,100,200,300,400),  cex=1.5)
    
    if (obs =='if') {mtext(2, text='match luminance', line=2, cex=2.2)}

    for (curr_med in contexts)
        {
        curr.fit = lm(bin_matches[[curr_med]]$scaled_match ~ bin_matches[[curr_med]]$scaled_contr)
        slopes[obs, curr_med]     = summary(curr.fit)$coefficients[2]
        intercepts[obs, curr_med] = summary(curr.fit)$coefficients[1]
        rsq[obs, curr_med]        = summary(curr.fit)$r.squared

        abline(curr.fit, col= unlist(all_cols[curr_med]), lwd=2)

        arrows(bin_matches[[curr_med]]$scaled_contr, (bin_matches[[curr_med]]$scaled_match - bin_matches[[curr_med]]$scaled_sem), bin_matches[[curr_med]]$scaled_contr, (bin_matches[[curr_med]]$scaled_match + bin_matches[[curr_med]]$scaled_sem), length=0.01, angle=90, code=3)

        points(bin_matches[[curr_med]]$scaled_contr, bin_matches[[curr_med]]$scaled_match, pch=21, bg=unlist(all_cols[curr_med]))
        }
    legend('topleft', c(obs), bty='n', cex=2)
    }

legend('bottomright', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red, my_orange), bty='n', pt.cex=0.9, cex=1.2)

dev.off()

write.table(format(global_rsq, digits=4), 'stat_out/norm_contrast_global_rsq.txt',       quote=FALSE)
write.table(format(slopes, digits=4),     'stat_out/norm_contrast_slopes.txt',    quote=FALSE)
write.table(format(intercepts, digits=4), 'stat_out/norm_contrast_intercepts.txt',quote=FALSE)
write.table(format(rsq, digits=4),        'stat_out/norm_contrast_rsq.txt',quote=FALSE)
