# clear history
rm(list = ls(all = TRUE))
graphics.off()

# import data
stock = as.data.frame(read.csv("stock10.csv", header = T, sep = ","))

# get daily return of each stock
price = stock[, 2:11]
value = matrix(stock[, 12], nrow = T)
T = nrow(price)
n = ncol(price)
W = matrix(unlist(price[1, ]/value[1, ]), nrow = n, ncol = 1)
alpha95 = 1.65
alpha99 = 2.33

return = log(price[2, ]/price[1, ]) * 100
for (j in 2:(T - 1)) {
    return[j, ] = log(price[j + 1, ]/price[j, ])
}



# install.packages('ccgarch') install.packages('tseries')
library("tseries")
library("ccgarch")

# use AR(1) model to get the residual of stock use time interval t=200, there are
# 247-200+1=48 intervals(T-t)

t = 200

VaR95 = matrix(, nrow = T - t, ncol = 1)
VaR99 = matrix(, nrow = T - t, ncol = 1)

for (j in 1:(T - t)) {
    
    ts = return[j:(j + t - 1), ]
    
    residual = matrix(, nrow = t, ncol = n)
    
    for (i in 1:n) {
        residual[, i] = matrix(residuals(arma(ts[, i], order = c(1, 0))))
    }
    
    residual = residual[-1, ]
    
    coef = matrix(, nrow = n, ncol = 3)
    
    # initial GARCH model estimation
    for (i in 1:10) {
        coef[i, ] = matrix(coef(garch(residual[, i], order = c(1, 1), series = NULL)))
        
    }
    
    # DCC-GARCH model estimation
    a = coef[, 1]
    A = diag(coef[, 2])
    B = diag(coef[, 3])
    dcc.para = c(0.01, 0.97)
    results = dcc.estimation(inia = a, iniA = A, iniB = B, ini.dcc = dcc.para, dvar = residual, 
        model = "diagonal")
    h = results$h
    dcc = results$DCC
    v = sqrt(diag(h[t - 1, ]))
    R = matrix(data = dcc[1, ], nrow = n, ncol = n)
    H = v %*% R %*% v
    
    VaR95[j, ] = sqrt(t(W) %*% H %*% W) * alpha95 * value[j]
    VaR99[j, ] = sqrt(t(W) %*% H %*% W) * alpha99 * value[j]
    
}

VaR = cbind(VaR95, VaR99)

# backtest profit&loss
value = matrix(value, nrow = T)
PL = value[(t + 1):T, ] - value[t:(T - 1), ]

# install.packages('parallel') install.packages('rugarch')
library("parallel")
library("rugarch")

ind_test = function(V) {
    J = matrix(ncol = 4, nrow = length(V))
    for (i in 2:length(V)) {
        J[i, 1] = V[i - 1] == 0 & V[i] == 0
        J[i, 2] = V[i - 1] == 0 & V[i] == 1
        J[i, 3] = V[i - 1] == 1 & V[i] == 0
        J[i, 4] = V[i - 1] == 1 & V[i] == 1
    }
    V_00 = sum(J[, 1], na.rm = TRUE)
    V_01 = sum(J[, 2], na.rm = TRUE)
    V_10 = sum(J[, 3], na.rm = TRUE)
    V_11 = sum(J[, 4], na.rm = TRUE)
    p_00 = V_00/(V_00 + V_01)
    p_01 = V_01/(V_00 + V_01)
    p_10 = V_10/(V_10 + V_11)
    p_11 = V_11/(V_10 + V_11)
    hat_p = (V_01 + V_11)/(V_00 + V_01 + V_10 + V_11)
    a = (1 - hat_p)^(V_00 + V_10) * (hat_p)^(V_01 + V_11)
    b = (p_00)^(V_00) * (p_01)^(V_01) * (p_10)^(V_10) * p_11^(V_11)
    return(-2 * log(a/b))
}

V95 = matrix(, nrow = (T - t))
for (i in 1:(T - t)) {
    if (PL[i] > (-VaR95[i, 1])) {
        V95[i] = 0
    } else {
        V95[i] = 1
    }
}

V99 = matrix(, nrow = (T - t))
for (i in 1:(T - t)) {
    if (PL[i] > (-VaR99[i, 1])) {
        V99[i] = 0
    } else {
        V99[i] = 1
    }
}

inde95 = ind_test(V95)
inde99 = ind_test(V99)

vt95 = VaRTest(alpha = 0.05, as.numeric(PL), as.numeric(-VaR95), conf.level = 0.95)
vt99 = VaRTest(alpha = 0.01, as.numeric(PL), as.numeric(-VaR99), conf.level = 0.95)
vt95
vt99

cc.stat95 = vt95$uc.LRstat + inde95
cc.stat99 = vt99$uc.critical + inde99

if (cc.stat95 > vt95$cc.critical) {
    print("Reject H0")
} else {
    print("Fail to Reject H0")
}

if (cc.stat99 > vt99$cc.critical) {
    print("Reject H0")
} else {
    print("Fail to Reject H0")
}

bt = cbind(PL, -VaR)
colnames(bt) = c("PL", "95%VaR", "99%VaR")
matplot(c(1:(T - t)), bt[, 1:3], type = "l", xlab = "time", ylab = "P&L")
legend("topright", colnames(bt)[-1], lwd = 1, col = 2:3)
title("Portfolio P&L and estimated VaR")


