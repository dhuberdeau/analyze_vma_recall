library(aod)
library(R.matlab)

RELATIVE_PT_MINMAX = 0.5;

setwd(getwd())

pc_data <- read.csv('raw_data_mat_E1_pc',
                    header = FALSE, na.strings = 'NaN')
head(pc_data)
# V1=correct or not, V2=trial type, V3=PT, V4=subject

dat <- data.frame(
    pc = pc_data$V1,
    type = pc_data$V2,
    pt = pc_data$V3,
    subject = pc_data$V4)

dat <- na.exclude(dat)

### --- Fit all trial types --- ###
pt <- dat$pt
y <- dat$pc
ty <- dat$type
ty[ty == 4] <- 3
ty[ty == 3] <- 3
ty <- as.factor(ty)
fit <- glm(y ~ pt*ty, data = data.frame(pt,ty,y), family = binomial(link = logit))

summary(fit)

with(dat[dat$type == 0,], plot(pc ~ pt, xlim = c(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX), col="red"))
par(new = TRUE)
with(dat[dat$type == 1,], plot(pc ~ pt, xlim = c(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX), col="green"))
par(new = TRUE)
with(dat[dat$type == 2,], plot(pc ~ pt, xlim = c(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX), col="blue"))
par(new = TRUE)
with(dat[dat$type == 3 | dat$type == 4,], plot(pc ~ pt, xlim = c(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX), col= "black"))

t <- seq(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX, .05)
p <- seq(0,0,length.out = length(t))
p <- as.factor(p)
y_0 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_0, col = "red")

p <- seq(1,1,length.out = length(t))
p <- as.factor(p)
y_1 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_1, col = "green")

p <- seq(2,2,length.out = length(t))
p <- as.factor(p)
y_2 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_2, col = "blue")

p <- seq(3,3,length.out = length(t))
p <- as.factor(p)
y_3 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_3, col = "black")

mat_out = data.frame(t, y_0, y_1, y_2, y_3)
writeMat("sigmoid_fit_E1.mat", mat_out = mat_out)
