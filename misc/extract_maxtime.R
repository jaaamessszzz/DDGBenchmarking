library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('extract_ddgtime.png', height=4096, width=4096, bg="white", res=600)
txtalpha <- 0.25
redtxtalpha <- 0.5

par(mar=c(5, 5, 1, 1))
a <- read.csv('extract_ddgtime.txt', header=T)
head(a)
coefs <- coef(lm(ChainLength~DDGTime, data = a))
# Sanity check
coefs
fitcoefs = coef(lm(ChainLength~0 + DDGTime, data = a))
fitlmv_ChainLength <- as.numeric(fitcoefs[1])

# coefs contains two values: (Intercept) and DDGTime
lmv_intercept <- as.numeric(coefs[1])
lmv_ChainLength <- as.numeric(coefs[2])

lm(a$ChainLength~a$DDGTime)
fitcoefs

xlabel <- "DDG time in minutes"
ylabel <- "Monomer chain length"
rvalue <- cor(a$ChainLength, a$DDGTime)

paste('PYTHON_VALUE', 'float', 'correlation', rvalue)

# To change the font size of the axis labels (tick labels), use e.g.:
# 	p <- p + theme(axis.text.x=element_text(size=22))
# To change the font of the axis titles, use e.g.:
# 	p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),

# shape I(20) is a small dot, I(19) is a large dot, I(4) is a cross

p <- qplot(DDGTime, ChainLength, main="ddg_monomer runtime vs chain length", data=a, xlab=xlabel, ylab=ylabel, shape = I(19), alpha = I(txtalpha)) + # label=ProThermID
		geom_abline(size = 0.25, intercept = lmv_intercept, slope = lmv_ChainLength) +
		geom_abline(color="blue",size = 0.25, intercept = 0, slope = fitlmv_ChainLength  )

# Create labels for cor(y,x) and MAE
# Using hjust=0 in geom_text sets text to be left-aligned

minx <- min(a$DDGTime)
maxx <- max(a$DDGTime)
miny <- min(a$ChainLength)
maxy <- max(a$ChainLength)

# fontface can be plain, bold, italic
# sans, serif, mono does not work for my PostScript driver
# palatino, bookman, helvetica, times works for PostScript but both look the same

fface <- "sans"

xpos <- minx + ((maxx - minx) * 0.05)
ypos_cor <- maxy - ((maxy - miny) * 0.015)
ypos_mae <- maxy - ((maxy - miny) * 0.085) # different variable name as these seem to be evaluated later (if we use the same label, if affects the printing of cor(y,x) as well)
p <- p + geom_text(hjust=0, size=6, aes(xpos, ypos_cor, fontface="plain", family = fface, label=sprintf("cor(y,x) = %f", round(rvalue, digits = 4))))

p <- p + geom_text(hjust=0, size=4, aes(xpos, ypos_mae, fontface="plain", family = fface, label="LIMITNOTE"))

aexp = a$DDGTime
apre = a$ChainLength

# Plot graph
p

dev.off()