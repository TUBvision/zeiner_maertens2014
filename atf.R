atf<-function(reflectance, condition)
    {
    ## assigns luminance to povray reflectance values for specified condition
    ## here practically reflectance and x are always identical
    atf = read.table(paste('luminance_measurements.txt', sep=''), head=T)
    atf.means = tapply(atf$luminance, data.frame(atf$reflectance, atf$medium), mean)
    x   = as.numeric(row.names(atf.means))
    atf.fit = lm(atf.means[,condition] ~ x)
    return(reflectance * atf.fit$coefficients[2] + atf.fit$coefficients[1])
    }

## color definition
my_orange    = '#FDAE61'
my_red       = '#D7191C'
my_darkblue  = '#2C7BB6'
my_lightblue = '#ABD9E9'
my_xlightblue= '#D3EBF4'

all_cols = list(plain = 'black', shadow='grey30', transp_med=my_red, transp_med_steep=my_orange, transp_dark=my_darkblue, transp_dark_steep=my_lightblue)


## combine log and design data
log_data = read.table(paste('data/log_if.txt', sep=''),header=TRUE)
x_values = sort(unique(log_data$target_lum))

## empty matrix for fit output
params  = matrix(0, nr=6, nc=3, dimn=list( names(all_cols), c('slope','intercept','rsq')))

## plot
svg('figs/atf.svg', w=6,h=6)
par(lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=3, font=1, font.axis=3, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

plot(x_values, atf(x_values,  'plain'),  pch=20, col='black', xlab='reflectance (rendered)', ylab='luminance', xlim=c(0,2.5), ylim=c(0,425))

## loop over linear fits and data points for all conditions
for (curr_med in unique(log_data$medium))
    {
    curr_lm = lm(atf(x_values, curr_med) ~ x_values)
    abline(curr_lm, col=unlist(all_cols[curr_med]))
    points(x_values, atf(x_values, curr_med), pch=20, col=unlist(all_cols[curr_med]))
    params[curr_med, 'slope']     = summary(curr_lm)$coefficients[2]
    params[curr_med, 'intercept'] = summary(curr_lm)$coefficients[1]
    params[curr_med, 'rsq']       = summary(curr_lm)$r.squared
    }

legend('topleft', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp low', 'light transp high'), pch=rep(20, 6), col=c('black', 'grey30', my_darkblue, my_lightblue, my_red, my_orange), bty='n')

dev.off()

write.table(format(params, digits=4), 'stat_out/atf_regress.txt', quote=FALSE)


## PLOT FOR REAL CHECKERBOARDS
## combine log and design data
## transp_dark_1 = N.6, 25%
## transp_dark_2 = N.3, 50%
## transp_dark_3 = N.12,75%

all_cols = list(plain = 'black', shadow='grey', transparency_light=my_red, transparency_dark_1=my_darkblue, transparency_dark_2=my_lightblue, transparency_dark_3=my_xlightblue)

data = read.table(paste('luminances_real_checkerboard.txt', sep=''),header=TRUE)

## empty matrix for fit output
params  = matrix(0, nr=6, nc=2, dimn=list( names(all_cols), c('slope','intercept')))
matches = tapply(data$luminance, data.frame(data$reflectance, data$medium), mean)

x_values = as.numeric(rownames(matches))
x_values = c(3.1,7.7,13.7,17.6,27.2,36.2,43.1,45.8,73.4)

svg('figs/atf_real.svg', w=6,h=6)

par(lwd=1.5, tck=.01, bty='n', cex=1.2, cex.axis=1.8, cex.lab=1.8, cex.main=1.8, font.lab=3, font=1, font.axis=3, font.main=1,  mar=c(3,3.5,2,1), mgp=c(2,0.5,0))

plot(x_values, matches[,'plain'],  pch=21, bg='white', col='black', xlab='reflectance', ylab='luminance', xlim=c(0,100), ylim=c(0,425), cex=.9)

for (curr_med in colnames(matches)[-5])
    {
    curr_lm = lm( matches[,curr_med] ~ x_values)
    abline(curr_lm, col=unlist(all_cols[curr_med]))
    points(x_values, matches[, curr_med], pch=21, bg=unlist(all_cols[curr_med]), col='black', cex=.9)
    params[curr_med, 'slope']     = summary(curr_lm)$coefficients[2]
    params[curr_med, 'intercept'] = summary(curr_lm)$coefficients[1]
    }

legend('topleft', c('plain', 'shadow', 'dark transp low', 'dark transp high', 'light transp'), pch=rep(21, 6), pt.bg=c('black', 'grey', my_darkblue, my_lightblue, my_red), bty='n')

dev.off()

write.table(format(params, digits=4), 'stat_out/atf_real_regress.txt', quote=FALSE)


