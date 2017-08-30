

library(nimble)
library(testthat)

test_that('unnecessary data do not break model building', {
    data <- list(a = 1, unwantedVariable = 2)

    toyModel <- nimbleModel(
        nimbleCode({
            a ~ dnorm(0,1)
        }),
        data = data
    )

    expect_true(inherits(toyModel, 'modelBaseClass'),
                'nimbleModel stopped due to unnecessary data.')
})

##c.toy.model <- compileNimble(toyModel)

## check if test works from package
library(testthat)
test_package('nimble', 'models', fixed = TRUE)

##

nimble:::modelBaseClass$trace('setData', browser)

ADCode1 <- nimbleCode({
    a[1] ~ dnorm(0, 1)
    a[2] ~ dnorm(0, 1)
    y[1] ~ dnorm(a[1], 1)
    y[2] ~ dnorm(a[2], 1)
})

ADMod1 <- nimbleModel(code = ADCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)))

cADMod1 <- compileNimble(ADMod1, showCompilerOutput = TRUE)
## derivs(calculate(model, nodes)) not fully implemented yet, this tests only that compilation is successful.

