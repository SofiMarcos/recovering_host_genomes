######
### R script
######

module load intel/perflibs/
module load gcc/9.3.0
module load R/4.0.0

setwd()
a <- read.table("CA03_14_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
d <- read.table("CA19_13_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
e <- read.table("CA19_18_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
g <- read.table("CA24_01_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
k <- read.table("CC19_18_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
l <- read.table("CA06_05_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
m <- read.table("CA13_03_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
u <- read.table("CA21_09_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
p <- read.table("CB20_03_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
z <- read.table("CC09_13_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
v <- read.table("CC20_14_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)
w <- read.table("CB07_06_imp.diff.discordance_matrix", sep = "", header = T, row.names = 1)

or <- a+d+e+g+k+l+m+u+p+z
# concordance
sum(or[1,1],or[2,2],or[3,3])/sum(or[1,],or[2,],or[3,])
sum(a[1,1],a[2,2],a[3,3])/sum(a[1,],a[2,],a[3,])
sum(d[1,1],d[2,2],d[3,3])/sum(d[1,],d[2,],d[3,])
sum(e[1,1],e[2,2],e[3,3])/sum(e[1,],e[2,],e[3,])
sum(g[1,1],g[2,2],g[3,3])/sum(g[1,],g[2,],g[3,])
sum(k[1,1],k[2,2],k[3,3])/sum(k[1,],k[2,],k[3,])
sum(l[1,1],l[2,2],l[3,3])/sum(l[1,],l[2,],l[3,])
sum(m[1,1],m[2,2],m[3,3])/sum(m[1,],m[2,],m[3,])
sum(u[1,1],u[2,2],u[3,3])/sum(u[1,],u[2,],u[3,])
sum(p[1,1],p[2,2],p[3,3])/sum(p[1,],p[2,],p[3,])
sum(z[1,1],z[2,2],z[3,3])/sum(z[1,],z[2,],z[3,])
sum(v[1,1],v[2,2],v[3,3])/sum(v[1,],v[2,],v[3,])
sum(w[1,1],w[2,2],w[3,3])/sum(w[1,],w[2,],w[3,])

# precision
sum(or[2,2],or[3,3])/sum(or[1,2],or[2,2],or[3,2],or[1,3],or[2,3],or[3,3])
sum(a[2,2],a[3,3])/sum(a[1,2],a[2,2],a[3,2],a[1,3],a[2,3],a[3,3])
sum(d[2,2],d[3,3])/sum(d[1,2],d[2,2],d[3,2],d[1,3],d[2,3],d[3,3])
sum(e[2,2],e[3,3])/sum(e[1,2],e[2,2],e[3,2],e[1,3],e[2,3],e[3,3])
sum(g[2,2],g[3,3])/sum(g[1,2],g[2,2],g[3,2],g[1,3],g[2,3],g[3,3])
sum(k[2,2],k[3,3])/sum(k[1,2],k[2,2],k[3,2],k[1,3],k[2,3],k[3,3])
sum(l[2,2],l[3,3])/sum(l[1,2],l[2,2],l[3,2],l[1,3],l[2,3],l[3,3])
sum(m[2,2],m[3,3])/sum(m[1,2],m[2,2],m[3,2],m[1,3],m[2,3],m[3,3])
sum(u[2,2],u[3,3])/sum(u[1,2],u[2,2],u[3,2],u[1,3],u[2,3],u[3,3])
sum(p[2,2],p[3,3])/sum(p[1,2],p[2,2],p[3,2],p[1,3],p[2,3],p[3,3])
sum(z[2,2],z[3,3])/sum(z[1,2],z[2,2],z[3,2],z[1,3],z[2,3],z[3,3])
sum(v[2,2],v[3,3])/sum(v[1,2],v[2,2],v[3,2],v[1,3],v[2,3],v[3,3])

# missingness
sum(or[4,]/sum(or[1,],or[2,],or[3,],or[4,]))
1-sum(a[4,]/sum(a))
1-sum(d[4,]/sum(d))
1-sum(e[4,]/sum(e))
1-sum(g[4,]/sum(g))
1-sum(k[4,]/sum(k))
1-sum(l[4,]/sum(l))
1-sum(m[4,]/sum(m))
1-sum(u[4,]/sum(u))
1-sum(p[4,]/sum(p))
1-sum(z[4,]/sum(z))
1-sum(v[4,]/sum(v))
1-sum(w[4,]/sum(w))

# missingness of heterozygous sites
1-sum(a[4,2]/sum(a[,2]))
1-sum(d[4,2]/sum(d[,2]))
1-sum(e[4,2]/sum(e[,2]))
1-sum(g[4,2]/sum(g[,2]))
1-sum(k[4,2]/sum(k[,2]))
1-sum(l[4,2]/sum(l[,2]))
1-sum(m[4,2]/sum(m[,2]))
1-sum(u[4,2]/sum(u[,2]))
1-sum(p[4,2]/sum(p[,2]))
1-sum(z[4,2]/sum(z[,2]))
1-sum(v[4,2]/sum(v[,2]))
1-sum(w[4,2]/sum(w[,2]))
