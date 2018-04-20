# CS112 Final Project

This paper is a extension and replication of “Providing driver’s licenses to unauthorized immigrants in California improves traffic safety” paper by Hans Lueders, Jens Hainmueller and Duncan Lawrence published in 2017 that examines the short-term effects of California’s Assembly Bill 60 (AB60) implemented in 2015. The law legalizes driving for unauthorized immigrants in California.

## Running the analysis and code
`FARS_data_preparer.R` prepares and cleans the data from sources like FARS, US Census Bureau and RI Department of Labor and Training for use of the analysis
`replication_code.R` performs the replication of model created by Lueders et al using the accident data provided by FARS data set.
`synthetic_control.R` performs synthetic control on the treated unit California and creates a synthetic control unit from the donor pool of other US states who have not passed a similar law.
`vcovCluster.R` creates clusters. It has been copied from the replication data by Lueders et al.