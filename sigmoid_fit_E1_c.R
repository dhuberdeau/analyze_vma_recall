library(aod)

pc_data <- read.csv('/Users/david/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/raw_data_mat_E1_pc',
                    header = FALSE, na.strings = 'NaN')
head(pc_data)
# V1=correct or not, V2=trial type, V3=PT, V4=subject

dat <- data.frame(
    pc = pc_data$V1,
    type = pc_data$V2,
    pt = pc_data$V3,
    subject = pc_data$V4)

dat <- na.exclude(dat)

### --- Fit No Pre-Cue trials --- ###
pt <- dat$pt[dat$type == 0 | dat$type == 3 | dat$type == 4]
y <- dat$pc[dat$type == 0 | dat$type == 3 | dat$type == 4]
ty <- dat$type[dat$type == 0 | dat$type == 3 | dat$type == 4]
ty[ty == 4] <- 1
ty[ty == 3] <- 1
ty <- as.factor(ty)
fit <- glm(y ~ pt*ty, data = data.frame(pt,ty,y), family = binomial(link = logit))

summary(fit)

with(dat[dat$type == 0,], plot(pc ~ pt, xlim = c(-.5, 1)))
par(new = TRUE)
with(dat[dat$type == 3 | dat$type == 4,], plot(pc ~ pt, xlim = c(-.5, 1), axes=FALSE, col= "red"))

t <- seq(-.5, 1, .05)
p <- seq(0,0,length.out = length(t))
p <- as.factor(p)
y_0 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_0)

p <- seq(1,1,length.out = length(t))
p <- as.factor(p)
y_1 <- predict(fit, newdata = data.frame(pt=t,ty=p), type="response")
lines(t, y_1, col = "red")
