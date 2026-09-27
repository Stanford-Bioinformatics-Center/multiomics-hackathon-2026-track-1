# Run the engine's tests: Rscript -e 'testthat::test_local("network/engine")'
library(testthat); library(exnet)
test_check("exnet")
